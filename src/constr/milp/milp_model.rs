use crate::constr::lp::LPModelStdform;
use crate::primitives::matrix::CSCMatrix;
use crate::primitives::vector::DenseVector;
use crate::primitives::vector::SparseVector;



#[derive(Clone, Copy)]
pub enum MILPConstrSense {
    Equal,
    Less,
    Greater,
}



// TODO: preprocessing - round int bounds
pub struct MILPModel {
    model_name: String,
    obj_coefs: DenseVector,
    obj_offset: f64,
    coef_matrix: CSCMatrix,
    constr_rhs: DenseVector,
    constr_sense: Vec<MILPConstrSense>,
    var_lb: Vec<Option<f64>>,
    var_ub: Vec<Option<f64>>,
    var_is_int: Vec<bool>,
    tol_feas: f64,
}

impl MILPModel {
    
    pub fn new(
            model_name: String,
            obj_coefs: DenseVector,
            obj_offset: f64,
            coef_matrix: CSCMatrix,
            constr_rhs: DenseVector,
            constr_sense: Vec<MILPConstrSense>,
            var_lb: Vec<Option<f64>>,
            var_ub: Vec<Option<f64>>,
            var_is_int: Vec<bool>,
        ) -> Self {

        // TODO: other invariants to enforce?
        let num_rows = coef_matrix.num_rows();
        let num_cols = coef_matrix.num_cols();
        assert_eq!(obj_coefs.len(), num_cols);
        assert_eq!(constr_rhs.len(), num_rows);
        assert_eq!(constr_sense.len(), num_rows);
        assert_eq!(var_lb.len(), num_cols);
        assert_eq!(var_ub.len(), num_cols);
        assert_eq!(var_is_int.len(), num_cols);
        
        Self {
            model_name,
            obj_coefs,
            obj_offset,
            coef_matrix,
            constr_rhs,
            constr_sense,
            var_lb,
            var_ub,
            var_is_int,
            tol_feas: 1e-6,  // TODO: parameterize
        }
    }

    pub fn obj_offset(&self) -> f64 {
        self.obj_offset
    }

    pub fn relax_to_lp_stdform(&self) -> (LPModelStdform, MILPStdformMappings) {
        // note: "lower" and "upper" bounds on negative-component variables
        // here refer to lower/upper bounds on their absolute value, e.g.
        // a variable ranging from -10 to -1 will have a negative component
        // with a lower bound of 1 and an upper bound of 10

        let num_vars_orig = self.coef_matrix.num_cols();
        
        // determine which variables need positive and negative components,
        // and what bound constraints are needed (map things out, don't edit yet)
        // IDEA: memory-conserving mode where these mappings aren't stored densely?
        let mut var_pos_col: Vec<Option<usize>> = vec![];
        let mut var_neg_col: Vec<Option<usize>> = vec![];
        let mut var_pos_lb_row: Vec<Option<usize>> = vec![];
        let mut var_pos_ub_row: Vec<Option<usize>> = vec![];
        let mut var_neg_lb_row: Vec<Option<usize>> = vec![];
        let mut var_neg_ub_row: Vec<Option<usize>> = vec![];
        let mut next_neg_idx: usize = self.coef_matrix.num_cols();
        let mut next_bnd_idx: usize = self.coef_matrix.num_rows();
        let mut bound_rows_indices: Vec<usize> = vec![];
        let mut bound_rows_rhs: Vec<f64> = vec![];
        let mut stdform_var_is_int: Vec<bool> = vec![false; self.coef_matrix.num_cols()];
        let mut add_bound_if = |p: bool, bnd_map: &mut Vec<Option<usize>>, c: usize, rhs: Option<f64>, neg: bool| {
            if p {
                bnd_map.push(Some(next_bnd_idx));
                next_bnd_idx += 1;
                let rhs = rhs.unwrap();
                bound_rows_rhs.push(if neg {-rhs} else {rhs});
                bound_rows_indices.push(c);
            }
            else {
                bnd_map.push(None);
            }
        };
        for c in 0..num_vars_orig {
            let lb = self.var_lb[c];
            let ub = self.var_ub[c];

            stdform_var_is_int[c] = self.var_is_int[c];

            if lb.is_none() || lb.unwrap() < 0.0 {
                if ub.is_none() || ub.unwrap() > 0.0 {
                    // split this column into positive and negative
                    var_pos_col.push(Some(c));
                    var_neg_col.push(Some(next_neg_idx));

                    add_bound_if(ub.is_some(), &mut var_pos_ub_row, c, ub, false);
                    add_bound_if(lb.is_some(), &mut var_neg_ub_row, next_neg_idx, lb, true);
                    var_pos_lb_row.push(None);
                    var_neg_lb_row.push(None);
                    
                    next_neg_idx += 1;
                    stdform_var_is_int.push(self.var_is_int[c]);
                }
                else {
                    // negative only (flip sign of this column)
                    var_pos_col.push(None);
                    var_neg_col.push(Some(c));

                    add_bound_if(lb.is_some(), &mut var_neg_ub_row, c, lb, true);
                    add_bound_if(ub.is_some_and(|x| x < 0.0), &mut var_neg_lb_row, c, ub, true);
                    var_pos_lb_row.push(None);
                    var_pos_ub_row.push(None);
                }
            }
            else {
                // positive-only (leave this column as is)
                var_pos_col.push(Some(c));
                var_neg_col.push(None);

                add_bound_if(ub.is_some(), &mut var_pos_ub_row, c, ub, false);
                add_bound_if(lb.is_some_and(|x| x > 0.0), &mut var_pos_lb_row, c, lb, false);
                var_neg_lb_row.push(None);
                var_neg_ub_row.push(None);
            }
        };

        // okay now actually modify (a clone of) the data structures
        // according to the mappings constructed above
        let mut coef_matrix = self.coef_matrix.clone();
        let mut obj_coefs = self.obj_coefs.clone();
        let mut constr_rhs = self.constr_rhs.clone();

        // add columns for negative components for positive-or-negative variables
        let mut cols_to_split: Vec<usize> = vec![];
        for c in 0..num_vars_orig {
            if var_pos_col[c].is_some() && var_neg_col[c].is_some() {
                cols_to_split.push(c);
                obj_coefs.append_inplace(-obj_coefs[c]);
            }
        }
        let num_vars_split = cols_to_split.len();
        coef_matrix.append_col_multiples_inplace(cols_to_split, -1.0);

        // flip sign of columns for negative-only variables
        for c in 0..num_vars_orig {
            if var_pos_col[c].is_none() && var_neg_col[c].is_some() {
                coef_matrix.scale_col_inplace(c, -1.0);
                obj_coefs.mul_value_inplace(c, -1.0);
            }
        }

        // add rows for variable bounds
        let bound_row_length: usize = num_vars_orig + num_vars_split;
        let bound_rows: Vec<SparseVector> = bound_rows_indices.iter().map(|c| {
            SparseVector::new(bound_row_length, vec![*c], vec![1.0])
        }).collect();
        coef_matrix.append_rows_inplace(bound_rows);
        for b in bound_rows_rhs {
            constr_rhs.append_inplace(b);
        }

        // add columns for slack/surplus variables, including added bound rows
        let mut slack_cols: Vec<SparseVector> = vec![];
        let mut add_slack_if_some = |r: Option<usize>, coef: f64| {
            if let Some(r) = r {
                slack_cols.push(SparseVector::new(coef_matrix.num_rows(), vec![r], vec![coef]));
                stdform_var_is_int.push(false);
            }
        };
        for r in 0..self.constr_sense.len() {
            match self.constr_sense[r] {
                MILPConstrSense::Equal => (),
                MILPConstrSense::Less => add_slack_if_some(Some(r), 1.0),
                MILPConstrSense::Greater => add_slack_if_some(Some(r), -1.0)
            }
        }
        for c in 0..num_vars_orig {
            add_slack_if_some(var_pos_ub_row[c], 1.0);
            add_slack_if_some(var_pos_lb_row[c], -1.0);
            add_slack_if_some(var_neg_ub_row[c], 1.0);
            add_slack_if_some(var_neg_lb_row[c], -1.0);
        }
        obj_coefs.append_zeros_inplace(slack_cols.len());
        coef_matrix.append_cols_inplace(slack_cols);

        // done; bundle everything up
        let stdform = LPModelStdform::new(
            obj_coefs,
            coef_matrix,
            constr_rhs,
        );
        let mappings = MILPStdformMappings {
            stdform_var_is_int,
        };

        (stdform, mappings)
    }

}



pub struct MILPStdformMappings {
    stdform_var_is_int: Vec<bool>,
}

impl MILPStdformMappings {

    pub fn stdform_var_is_int(&self) -> &Vec<bool> {
        &self.stdform_var_is_int
    }

}
