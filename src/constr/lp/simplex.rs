// TODO: testing with known-degenerate models



use crate::constr::lp::{
    LPModelStdform,
    LPStdformResult,
    LPStdformResultInfeasible,
    LPStdformResultOptimal,
    LPStdformResultUnbounded,
};
use crate::linalg::csparse::{
    CSparseLU,
    csparse_lu,
    csparse_solve_with_lu,
    csparse_solve_with_lu_transpose,
};
use crate::primitives::vector::{
    DenseVector,
    SparseVector,
};
use crate::util::indices_complement;



pub fn solve_dual_simplex_warm(model: &LPModelStdform, basis: Vec<usize>) -> LPStdformResult {
    solve_dual_simplex(model, basis)
}



/// May leave added columns and rows for degenerate artificial variables.
pub fn solve_primal_simplex_cold(model: &mut LPModelStdform) -> LPStdformResult {
    solve_primal_simplex_twophase_artif_vars(model)
}



fn solve_dual_simplex(model: &LPModelStdform, starting_basis: Vec<usize>) -> LPStdformResult {
    // TODO: refactor pivoting into some sort of "single pivot" function?
    // TODO: > and call with different arguments for a degenerate pivot? different function?

    let mut basis = starting_basis;
    let mut iter_count: usize = 0;

    loop {
        iter_count += 1;

        // extract the basis systems (same as primal simplex)
        let basis_matrix = model.coef_matrix().col_submatrix(&basis);
        let basis_obj_coefs = model.obj_coefs().subvector(&basis);
        let primal_rhs = model.constr_rhs().clone();
        let dual_rhs = basis_obj_coefs.clone();

        // solve the basis systems (same as primal simplex)
        let lu: CSparseLU = csparse_lu(&basis_matrix);
        let primal_vals = csparse_solve_with_lu(&lu, primal_rhs);
        let dual_vals = csparse_solve_with_lu_transpose(&lu, dual_rhs);
        let obj_val = basis_obj_coefs.dot(&primal_vals);

        // calculate dual slack (a.k.a. reduced costs)
        let nonbasis = indices_complement(&basis, model.coef_matrix().num_cols());
        let dual_slack: DenseVector = DenseVector::new(
            nonbasis.iter().map(|c|
                model.obj_coefs()[*c] - model.coef_matrix().col_dot(*c, &dual_vals)
            ).collect()
        );

        // validation: on first iteration, verify that basis is dual feasible
        // (dual slack should be nonnegative--satisfied by construction for basic columns)
        if iter_count == 1 {
            let neg_vals = dual_slack.values().iter().filter(|x| **x < 0.0).count();
            if neg_vals > 0 {
                panic!("Starting basis passed to solve_dual_simplex is not dual feasible!");
            }
        }

        // pick a negative primal variable to make zero by exiting the basis
        // (current implementation here: pick the most negative one)
        let (min_prim_bidx, min_prim_val) = primal_vals.min_elem();
        if min_prim_val >= 0.0 {
            // if all nonnegative, this is primal feasible and we're done
            return LPStdformResultOptimal::new(
                obj_val,
                primal_vals,
                basis,
                iter_count,
            ).into();
        }
        let bidx_exiting = min_prim_bidx;

        // solve for the dual-variable direction (which will increase the
        // dual slack for c_exiting but leave the dual slack at zero for
        // all the other columns remaining in the basis)
        let dual_dir_rhs: DenseVector = DenseVector::new_unit(basis.len(), bidx_exiting);
        let dual_dir = csparse_solve_with_lu_transpose(&lu, dual_dir_rhs);

        // calculate the dual slack velocity (change per unit of dual slack
        // added for the exiting variable)
        let dual_slack_velocity: DenseVector = DenseVector::new(
            nonbasis.iter().map(|c|
                -model.coef_matrix().col_dot(*c, &dual_dir)
            ).collect()
        );

        // calculate amount the exiting variable's dual slack can increase
        // before each candidate variable's dual slack hits zero
        // TODO: parameterize zero-tolerance? hard to find discussions of this
        let dist_nidxs: Vec<usize> = (0..nonbasis.len())
            .filter(|nidx| dual_slack_velocity[*nidx] > 1e-12).collect();
        if dist_nidxs.len() == 0 {
            // dual is unbounded, meaning primal is infeasible
            return LPStdformResultInfeasible {}.into();
        }
        let dist: DenseVector = DenseVector::new(
            dist_nidxs.iter().map(|nidx|
                dual_slack[*nidx] / dual_slack_velocity[*nidx]
            ).collect()
        );
        let (didx_entering, min_dist) = dist.min_elem();
        let nidx_entering = dist_nidxs[didx_entering];

        // TODO: degeneracy and cycle prevention (if min_dist is zero)
        basis[bidx_exiting] = nonbasis[nidx_entering];
    }
}



fn solve_primal_simplex(model: &LPModelStdform, starting_basis: Vec<usize>) -> LPStdformResult {
    // TODO: refactor pivoting into some sort of "single pivot" function?
    // TODO: > and call with different arguments for a degenerate pivot? different function?
    
    let mut basis = starting_basis;
    let mut iter_count: usize = 0;

    loop {
        iter_count += 1;

        // extract the basis systems
        let basis_matrix = model.coef_matrix().col_submatrix(&basis);
        let basis_obj_coefs = model.obj_coefs().subvector(&basis);
        let primal_rhs = model.constr_rhs().clone();
        let dual_rhs = basis_obj_coefs.clone();

        // solve the basis systems
        let lu: CSparseLU = csparse_lu(&basis_matrix);
        let primal_vals = csparse_solve_with_lu(&lu, primal_rhs);
        let dual_vals = csparse_solve_with_lu_transpose(&lu, dual_rhs);
        let obj_val = basis_obj_coefs.dot(&primal_vals);

        // validation: on first iteration, verify that basis is primal feasible
        if iter_count == 1 {
            let neg_vals = primal_vals.values().iter().filter(|x| **x < 0.0).count();
            if neg_vals > 0 {
                panic!("Starting basis passed to solve_primal_simplex is not primal feasible!");
            }
        }

        // calculate reduced costs
        let nonbasis = indices_complement(&basis, model.coef_matrix().num_cols());
        let rcosts: DenseVector = DenseVector::new(
            nonbasis.iter().map(|c|
                model.obj_coefs()[*c] - model.coef_matrix().col_dot(*c, &dual_vals)
            ).collect()
        );

        // find most-negative reduced cost, or terminate if all nonnegative
        // (most-negative: "Dantzig's rule")
        let (min_rcost_nidx, min_rcost_value) = rcosts.min_elem();
        if min_rcost_value >= 0.0 {
            return LPStdformResultOptimal::new(
                obj_val,
                primal_vals,
                basis,
                iter_count,
            ).into();
        }
        let c_entering = nonbasis[min_rcost_nidx];

        // solve for pivot direction
        // TODO: do this with cs_spsolve and a sparse rhs instead
        let dir_rhs = model.coef_matrix().col_dense(c_entering);
        let dir = csparse_solve_with_lu(&lu, dir_rhs);

        // identify exiting variable: which one hits zero first along the pivot direction
        // (if none ever hit zero in this direction, we found an unbounded ray!)
        let dist_bidxs: Vec<usize> = (0..basis.len()).filter(|bidx| dir[*bidx] > 0.0).collect();
        if dist_bidxs.len() == 0 {
            assert!(min_rcost_value < 0.0, "Problem in exiting variable logic!");  // refactor guardrail

            let mut unb_indices: Vec<usize> = basis.clone();
            let mut unb_values = dir.clone().to_vec();
            unb_indices.push(c_entering);
            unb_values.push(1.0);

            return LPStdformResultUnbounded::new(
                SparseVector::new_without_zeros(model.coef_matrix().num_cols(), unb_indices, unb_values)
            ).into();
        }
        let dist: DenseVector = DenseVector::new(
            dist_bidxs.iter().map(|bidx| primal_vals[*bidx] / dir[*bidx]).collect()
        );
        let (didx_exiting, min_dist) = dist.min_elem();
        let bidx_exiting = dist_bidxs[didx_exiting];

        if min_dist > 0.0 {
            // update basis and repeat
            basis[bidx_exiting] = c_entering;
        }
        else {
            // degenerate step; apply Bland's rule (this is not efficient!):
            // select first nonbasic variable with negative reduced cost as the entering variable,
            // select first basic variable eligible to drop as the exiting variable
            let nidx_entering = rcosts.values().iter().enumerate()
                .filter(|(_, x)| **x < 0.0).next().unwrap().0;
            let c_entering = nonbasis[nidx_entering];
            let dir_rhs = model.coef_matrix().col_dense(c_entering);
            let dir = csparse_solve_with_lu(&lu, dir_rhs);
            let dist_bidxs: Vec<usize> = (0..basis.len()).filter(|bidx| dir[*bidx] > 0.0).collect();
            let dist: DenseVector = DenseVector::new(
                dist_bidxs.iter().map(|bidx| primal_vals[*bidx] / dir[*bidx]).collect()
            );
            let (didx_exiting, _) = dist.min_elem();
            let bidx_exiting = dist_bidxs[didx_exiting];

            basis[bidx_exiting] = c_entering;

            // TODO: unit testing for known-degenerate models
        }
    }
}



fn solve_primal_simplex_twophase_artif_vars(model: &mut LPModelStdform) -> LPStdformResult {
    let num_vars_before = model.coef_matrix().num_cols();
    let (basis, obj_before) = model.phase1_add_artif_vars_inplace();
    
    let phase1_result = solve_primal_simplex(model, basis);

    match phase1_result {
        LPStdformResult::Infeasible(_) => panic!("Simplex phase 1 infeasible (should not be possible)"),
        LPStdformResult::Unbounded(_) => panic!("Simplex phase 1 unbounded (should not be possible)"),
        LPStdformResult::Optimal(res) => {
            // if phase 1 is feasible with positive objective value, original is infeasible
            if res.obj_value() > 0.0 {
                // TODO: feasibility tolerance?
                return LPStdformResultInfeasible{}.into();
            }

            // phase 1 successful, proceed to phase 2
            let phase2_basis = model.phase1_postproc_artif_vars_inplace(
                num_vars_before, obj_before, &res.basis()
            );
            return solve_primal_simplex(model, phase2_basis);
        }
    }
}
