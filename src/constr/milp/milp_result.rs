use crate::primitives::vector::SparseVector;



pub enum MILPResult where {
    Infeasible,
    Optimal,
    Unbounded(MILPResultUnbounded),
}

impl From<MILPResultUnbounded> for MILPResult {
    fn from(value: MILPResultUnbounded) -> Self {
        MILPResult::Unbounded(value)
    }
}



pub struct MILPResultUnbounded {
    pub unbounded_ray: SparseVector<f64>,
}
