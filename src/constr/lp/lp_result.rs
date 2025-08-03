use crate::primitives::vector::DenseVector;
use crate::primitives::vector::SparseVector;



pub enum LPStdformResult {
    Infeasible(LPStdformResultInfeasible),
    Optimal(LPStdformResultOptimal),
    Unbounded(LPStdformResultUnbounded),
}

impl From<LPStdformResultOptimal> for LPStdformResult {
    fn from(value: LPStdformResultOptimal) -> Self {
        LPStdformResult::Optimal(value)
    }
}

impl From<LPStdformResultInfeasible> for LPStdformResult {
    fn from(value: LPStdformResultInfeasible) -> Self {
        LPStdformResult::Infeasible(value)
    }
}

impl From<LPStdformResultUnbounded> for LPStdformResult {
    fn from(value: LPStdformResultUnbounded) -> Self {
        LPStdformResult::Unbounded(value)
    }
}



pub struct LPStdformResultInfeasible {
    // ... TODO
}



pub struct LPStdformResultOptimal {
    obj_value: f64,
    basic_var_vals: DenseVector,
    basis: Vec<usize>,
    iterations: usize,
    // ... TODO
}

impl LPStdformResultOptimal {

    pub fn new(
            obj_value: f64,
            basic_var_vals: DenseVector,
            basis: Vec<usize>,
            iterations: usize,
        ) -> Self {
        
        Self {
            obj_value,
            basic_var_vals,
            basis,
            iterations,
        }
    }

    pub fn obj_value(&self) -> f64 {
        self.obj_value
    }

    pub fn basic_var_vals(&self) -> &DenseVector {
        &self.basic_var_vals
    }

    pub fn basis(&self) -> &Vec<usize> {
        &self.basis
    }

    pub fn iterations(&self) -> usize {
        self.iterations
    }

}



pub struct LPStdformResultUnbounded {
    unbounded_ray: SparseVector,
}

impl LPStdformResultUnbounded {

    pub fn new(unbounded_ray: SparseVector) -> Self {
        Self {
            unbounded_ray,
        }
    }

    pub fn unbounded_ray(&self) -> &SparseVector {
        &self.unbounded_ray
    }

}
