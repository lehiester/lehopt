use crate::primitives::matrix::CSCMatrix;
use crate::primitives::vector::DenseVector;
use crate::primitives::vector::SparseVector;



pub struct LPModelStdform {
    obj_coefs: DenseVector,  // no offset!
    coef_matrix: CSCMatrix,
    constr_rhs: DenseVector,
}

impl LPModelStdform {
    
    pub fn new(
            obj_coefs: DenseVector,
            coef_matrix: CSCMatrix,
            constr_rhs: DenseVector,
        ) -> Self {

        assert_eq!(obj_coefs.len(), coef_matrix.num_cols(), "LPModelStdform::new - objective vector does not match number of columns!");
        assert_eq!(constr_rhs.len(), coef_matrix.num_rows(), "LPModelStdform::new - rhs vector does not match number of rows!");

        Self {
            obj_coefs,
            coef_matrix,
            constr_rhs,
        }
    }

    // Getters

    #[inline]
    pub fn coef_matrix(&self) -> &CSCMatrix {
        &self.coef_matrix
    }

    #[inline]
    pub fn constr_rhs(&self) -> &DenseVector {
        &self.constr_rhs
    }

    #[inline]
    pub fn obj_coefs(&self) -> &DenseVector {
        &self.obj_coefs
    }

    // Instance methods

    pub fn append_cols_inplace(&mut self, cols: Vec<SparseVector>, obj: Vec<f64>) {
        assert_eq!(cols.len(), obj.len());

        self.coef_matrix.append_cols_inplace(cols);
        self.obj_coefs.extend_inplace(obj);
    }

    pub fn append_rows_inplace(&mut self, rows: Vec<SparseVector>, rhs: Vec<f64>) {
        assert_eq!(rows.len(), rhs.len());

        self.coef_matrix.append_rows_inplace(rows);
        self.constr_rhs.extend_inplace(rhs);
    }

    /// Adds artificial variables to this model and returns an initial
    /// primal-feasible basis, which may include non-artificial variables if a
    /// full set was not added due to the presence of other suitable ones.
    /// Also switches to the phase 1 objective and returns the original one.
    /// 
    /// Artificial variables are appended to the end; indices for all existing
    /// variables will remain unchanged.
    pub fn phase1_add_artif_vars_inplace(&mut self) -> (Vec<usize>, DenseVector) {
        // NOTE: current implementation just adds a full set
        // TODO: don't add artificial variables for rows with suitable slack/other variables

        let m = self.coef_matrix.num_rows();
        let n = self.coef_matrix.num_cols();

        // add columns for artificial variables
        let coef = |r| {
            if self.constr_rhs.values()[r] < 0.0 {-1.0} else {1.0}
        };
        let new_cols: Vec<SparseVector> = (0..m)
            .map(|r| SparseVector::new(m, vec![r], vec![coef(r)])).collect();
        self.coef_matrix.append_cols_inplace(new_cols);

        // swap in phase 1 objective
        let mut obj_phase1: DenseVector = DenseVector::new_zeros(n);
        obj_phase1.append_rep_inplace(1.0, m);
        let obj_before = std::mem::replace(&mut self.obj_coefs, obj_phase1);

        // sanity checks on matrix/vector dimensions
        assert_eq!(self.coef_matrix.num_cols(), self.obj_coefs.len());
        assert_eq!(self.coef_matrix.num_rows(), self.constr_rhs.len());

        // return original objective and starting basis
        ((n..(n + m)).collect(), obj_before)
    }

    /// Removes artificial variables from this model and returns the starting
    /// basis for phase 2, after reindexing if needed. Degenerate artificial
    /// variables that were in the ending basis for phase 1 will be kept and
    /// rows added to constrain them to zero; they will typically be reindexed,
    /// as any nonbasic artificial columns to their left will be removed.
    /// 
    /// Any remaining artificial variables will be at the end; indices for all
    /// original variables will remain unchanged.
    pub fn phase1_postproc_artif_vars_inplace(
            &mut self,
            num_orig_vars: usize,
            orig_obj: DenseVector,
            phase1_ending_basis: &Vec<usize>
        ) -> Vec<usize> {
        
        // remove nonbasic artificial variables
        let num_artif_vars = self.coef_matrix.num_cols() - num_orig_vars;
        let mut artif_var_in_basis = vec![false; num_artif_vars];
        for c in phase1_ending_basis {
            if *c >= num_orig_vars {
                artif_var_in_basis[*c - num_orig_vars] = true;
            }
        }
        let vars_to_remove: Vec<usize> = (0..artif_var_in_basis.len())
            .filter(|i| !artif_var_in_basis[*i])
            .map(|i| i + num_orig_vars)
            .collect();
        if vars_to_remove.len() > 0 {
            self.coef_matrix.remove_cols_inplace(&vars_to_remove);
        }
        self.obj_coefs = orig_obj;

        // constrain basic (degenerate) artificial variables to zero
        let num_artif_basic = self.coef_matrix.num_cols() - num_orig_vars;
        if num_artif_basic > 0 {
            self.obj_coefs.append_zeros_inplace(num_artif_basic);
            let artif_fix_rows: Vec<SparseVector> = (num_orig_vars..self.coef_matrix.num_cols())
                .map(|c| SparseVector::new(self.coef_matrix.num_cols(), vec![c], vec![1.0]))
                .collect();
            self.coef_matrix.append_rows_inplace(artif_fix_rows);
            self.constr_rhs.append_zeros_inplace(num_artif_basic);
        }

        // sanity checks on matrix/vector dimensions
        assert_eq!(self.coef_matrix.num_cols(), self.obj_coefs.len());
        assert_eq!(self.coef_matrix.num_rows(), self.constr_rhs.len());

        // construct phase 2 basis
        let mut phase2_basis: Vec<usize> = phase1_ending_basis.iter().copied()
            .filter(|c| *c < num_orig_vars).collect();
        for i in 0..num_artif_basic {
            phase2_basis.push(num_orig_vars + i);
        }

        phase2_basis
    }

    pub fn pop_cols_inplace(&mut self, num_to_pop: usize) {
        assert!(num_to_pop < self.coef_matrix.num_cols());

        self.coef_matrix.pop_cols_inplace(num_to_pop);
        self.obj_coefs.trunc_inplace(self.coef_matrix.num_cols());
    }

    pub fn pop_rows_inplace(&mut self, num_to_pop: usize) {
        assert!(num_to_pop < self.coef_matrix.num_rows());

        self.coef_matrix.pop_rows_inplace(num_to_pop);
        self.constr_rhs.trunc_inplace(self.coef_matrix.num_rows());
    }

}
