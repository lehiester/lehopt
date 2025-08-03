mod bnb_node;

use bnb_node::{BnbNode, BnbNodeSortedByBound};

use crate::constr::lp::{
    LPModelStdform,
    LPStdformResult,
    LPStdformResultOptimal,
    solve_dual_simplex_warm,
    solve_primal_simplex_cold,
};
use crate::constr::milp::{
    MILPConstrSense,
    MILPModel,
    MILPResult,
    MILPResultUnbounded,
};
use crate::primitives::vector::{
    DenseVector,
    SparseVector,
};
use crate::util::{
    AvlTree,
    pos_int_is_power_of,
};

use std::time::{Duration, Instant};



/// State of the branch-and-bound procedure. Lifetime tied to the model.
struct BnbState<'mdl> {
    model: &'mdl MILPModel,
    relaxation: LPModelStdform,
    incmb_result: Option<LPStdformResultOptimal>,
    nodes: AvlTree<BnbNodeSortedByBound>,
    nodes_expl_count: usize,  // incremented before processing, =1 while processing root
    next_node_id_to_assign: usize,
}

impl<'mdl> BnbState<'mdl> {

    pub fn get_incmb_val(&self) -> Option<f64> {
        if self.incmb_result.is_some() {
            Some(self.incmb_result.as_ref().unwrap().obj_value() + self.model.obj_offset())
        }
        else {
            None
        }
    }

}



pub fn solve_milp_bnb(model: &MILPModel) -> MILPResult {
    // initialization
    let time_start = Instant::now();
    let mut time_next_statmsg = time_start + Duration::new(5, 0);
    let (relaxation, mappings) = model.relax_to_lp_stdform();
    let mut state = BnbState {
        model: &model,
        relaxation: relaxation,
        incmb_result: None,
        nodes: AvlTree::new(),
        nodes_expl_count: 0,
        next_node_id_to_assign: 1,
    };
    let is_var_int = |c: usize| {
        c < mappings.stdform_var_is_int().len()
        && mappings.stdform_var_is_int()[c]
    };

    state.nodes.insert(BnbNode::new(0, f64::INFINITY, None, vec![], vec![]).into());

    // main loop: process each node
    while state.nodes.len() > 0 {
        state.nodes_expl_count += 1;

        // "best-first", i.e. most loose bound first
        let node: BnbNode = state.nodes.pop_first().unwrap().into();

        // before solving, check if a new incumbent is better than the parent's bound
        if state.get_incmb_val().is_some_and(|z| z <= node.parent_obj_val()) {
            continue;
        }

        let node_result = solve_node(&mut state, &node);

        match node_result {
            // unbounded case: record unbounded ray and terminate
            LPStdformResult::Unbounded(unb_result) => {
                return MILPResultUnbounded {
                    unbounded_ray: unb_result.unbounded_ray().clone(),  // TODO: map to original space
                }.into();
            },

            // infeasible case: prune, skip to next node
            LPStdformResult::Infeasible(_) => {
                continue;
            },

            // lp-feasible case: check if bounded above incumbent, then for integrality
            LPStdformResult::Optimal(result) => {
                if state.get_incmb_val().is_some_and(|z| z <= result.obj_value()) {
                    // prune by bound
                    continue;
                }

                let basis_and_values = std::iter::zip(
                    result.basis().iter(),
                    result.basic_var_vals().values().iter()
                );

                // find largest integrality violation
                let mut max_frac = 0.0;
                let mut max_frac_idx = 0;
                let mut max_frac_val = 0.0;
                for (var_idx, value) in basis_and_values {
                    if is_var_int(*var_idx) {
                        let fractionality = (value.round() - value).abs();
                        if fractionality > max_frac {
                            max_frac = fractionality;
                            max_frac_idx = *var_idx;
                            max_frac_val = *value;
                        }
                    }
                }

                // branch on most-fractional variable if above integrality tolerance
                // TODO: parameterize integrality tolerance
                if max_frac > 1e-9 {
                    // branch on the chosen variable
                    let chosen_var = max_frac_idx;
                    let node_left = node.make_child(
                        state.next_node_id_to_assign,
                        result.obj_value(),
                        result.basis().clone(),
                        vec![(chosen_var, MILPConstrSense::Less, max_frac_val.floor())]
                    );
                    let node_right = node.make_child(
                        state.next_node_id_to_assign + 1,
                        result.obj_value(),
                        result.basis().clone(),
                        vec![(chosen_var, MILPConstrSense::Greater, max_frac_val.ceil())]
                    );
                    state.nodes.insert(node_left.into());
                    state.nodes.insert(node_right.into());
                    state.next_node_id_to_assign += 2;
                }
                else {
                    // no fractional variables: integer-feasible, check against incumbent
                    let node_obj_value = result.obj_value();
                    if !state.get_incmb_val().is_some_and(|z| z <= node_obj_value) {
                        state.incmb_result = Some(result);
                    }
                }
            }
        }

        // periodically print status message
        if Instant::now() > time_next_statmsg {
            time_next_statmsg += Duration::new(5, 0);

            let incmb_descr = if let Some(z) = state.get_incmb_val() {
                format!("{:10.6e}", z)
            } else {
                "    --    ".to_string()
            };

            let bound_descr = if state.nodes.len() > 0 {
                format!("{:10.6e}", state.nodes.peek_first().unwrap().node().parent_obj_val())
            } else {
                incmb_descr.clone()
            };

            println!("{} {} {:4}s", incmb_descr, bound_descr, time_start.elapsed().as_secs());
        }
    }

    // main loop finished; process result
    match state.incmb_result {
        // if no incumbent, model instance is infeasible
        None => MILPResult::Infeasible,

        // if an incumbent value, return it
        // TODO: need to map standard-form space back to original space
        Some(res) => { todo!("Finished branch and bound with optimal result of {}", res.obj_value()); }
    }
}



fn solve_node(state: &mut BnbState, node: &BnbNode) -> LPStdformResult {
    // IDEA: stdform: experiment with lazy-adding placeholder rows for bound changes to avoid reshuffling CSCMatrix
    // IDEA: stdform: keep track of associated bound constraints/slack and update in-place

    // -- end of notes --

    // TODO: if lower bound added for pos/neg component of a split variable, constrain other to zero (no slack)
    // TODO: ^ is this necessary?

    // apply bound changes for this node
    let num_rows_before = state.relaxation.coef_matrix().num_rows();
    let num_cols_before = state.relaxation.coef_matrix().num_cols();
    let mut next_row_idx = num_rows_before;
    let mut bound_rows: Vec<usize> = vec![];  // column in which the single "1" occurs (before slack added)
    let mut bound_rhs: Vec<f64> = vec![];  
    let mut bound_slack_cols: Vec<(usize,f64)> = vec![];  // (row, coef); one item per column
    let all_bounds = node.parent_bounds().iter().chain(node.new_bounds());
    for (var_idx, sense, rhs) in all_bounds {
        bound_rows.push(*var_idx);
        bound_rhs.push(*rhs);
        let slack_coef = match sense {
            MILPConstrSense::Greater => -1.0,
            MILPConstrSense::Less => 1.0,
            MILPConstrSense::Equal => unreachable!(),
        };
        bound_slack_cols.push((next_row_idx, slack_coef));
        next_row_idx += 1;
    }
    let bound_rows: Vec<SparseVector> = bound_rows.into_iter()
        .map(|c| SparseVector::new(num_cols_before, vec![c], vec![1.0])).collect();
    let bound_slack_cols: Vec<SparseVector> = bound_slack_cols.into_iter()
        .map(|(r, coef)| SparseVector::new(next_row_idx, vec![r], vec![coef])).collect();
    let num_bound_rows = bound_rows.len();
    let num_bound_cols = bound_slack_cols.len();
    state.relaxation.append_rows_inplace(bound_rows, bound_rhs);
    state.relaxation.append_cols_inplace(bound_slack_cols, DenseVector::new_zeros(num_bound_cols).to_vec());


    // solve this node
    let lp_result;
    if node.parent_basis().is_none() {
        // assuming that introduction of artificial variables will only happen
        // at root, or more precisely where there is no issue of them getting
        // entangled with bound rows or bound slack columns
        assert_eq!(num_bound_rows, 0);
        assert_eq!(num_bound_cols, 0);
        assert_eq!(state.nodes_expl_count, 1);
        lp_result = solve_primal_simplex_cold(&mut state.relaxation);
    }
    else {
        // add slack variables for new bounds to the basis
        let num_cols_after = state.relaxation.coef_matrix().num_cols();
        let mut basis = node.parent_basis().clone().unwrap();
        basis.extend((num_cols_after - node.new_bounds().len())..num_cols_after);
        assert_eq!(basis.len(), state.relaxation.coef_matrix().num_rows());

        lp_result = solve_dual_simplex_warm(&state.relaxation, basis);
    }

    // remove bound rows and (slack) columns
    // TODO: don't remove and re-add if visiting a child of this node next
    state.relaxation.pop_cols_inplace(num_bound_cols);
    state.relaxation.pop_rows_inplace(num_bound_rows);

    lp_result
}
