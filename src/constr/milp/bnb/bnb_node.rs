use core::f64;

use crate::constr::milp::MILPConstrSense;
use crate::primitives::scalar::{NonNan, Scalar};
use crate::util::AvlTreeItem;



pub struct BnbNode {
    id: usize,
    parent_obj_val: NonNan<f64>,
    parent_basis: Option<Vec<usize>>,
    parent_bounds: Vec<(usize,MILPConstrSense,f64)>,
    new_bounds: Vec<(usize,MILPConstrSense,f64)>,
}

impl BnbNode {

    pub fn new(
            id: usize,
            parent_obj_val: f64,
            parent_basis: Option<Vec<usize>>,
            parent_bounds: Vec<(usize,MILPConstrSense,f64)>,
            new_bounds: Vec<(usize,MILPConstrSense,f64)>,
        ) -> Self {

        let obj_nn = parent_obj_val.non_nan()
            .unwrap_or_else(|| panic!("Node {}: parent objective value was NaN", id));
        
        Self {
            id,
            parent_obj_val: obj_nn,
            parent_basis,
            parent_bounds,
            new_bounds,
        }
    }

    pub fn parent_obj_val(&self) -> f64 {
        self.parent_obj_val.value()
    }

    pub fn parent_basis(&self) -> &Option<Vec<usize>> {
        &self.parent_basis
    }

    pub fn parent_bounds(&self) -> &Vec<(usize,MILPConstrSense,f64)> {
        &self.parent_bounds
    }

    pub fn new_bounds(&self) -> &Vec<(usize,MILPConstrSense,f64)> {
        &self.new_bounds
    }

    pub fn make_child(
            &self,
            child_id: usize,
            obj_val: f64,
            basis: Vec<usize>,
            child_bounds: Vec<(usize,MILPConstrSense,f64)>,
        ) -> Self {

        let mut self_all_bounds = self.parent_bounds.clone();
        self_all_bounds.extend(&self.new_bounds);

        let obj_nn = obj_val.non_nan()
            .unwrap_or_else(|| panic!("Node {}: parent objective value was NaN", child_id));

        Self {
            id: child_id,
            parent_obj_val: obj_nn,
            parent_basis: Some(basis),
            parent_bounds: self_all_bounds,
            new_bounds: child_bounds,
        }
    }

}

impl From<BnbNodeSortedByBound> for BnbNode {
    fn from(bnsbb: BnbNodeSortedByBound) -> Self {
        bnsbb.node
    }
}



pub struct BnbNodeSortedByBound {
    node: BnbNode,
}

impl BnbNodeSortedByBound {
    pub fn node(&self) -> &BnbNode {
        &self.node
    }
}

impl AvlTreeItem for BnbNodeSortedByBound {
    type Key = (NonNan<f64>, usize);

    fn key(&self) -> Self::Key {
        (self.node.parent_obj_val, self.node.id)
    }
}

/// Temporary(?), pending AvlTree rework
impl Default for BnbNodeSortedByBound {

    fn default() -> Self {
        Self {
            node: BnbNode {
                id: usize::MAX,
                parent_obj_val: f64::INFINITY.non_nan().unwrap(),
                parent_basis: None,
                parent_bounds: vec![],
                new_bounds: vec![],
            }
        }
    }

}

impl From<BnbNode> for BnbNodeSortedByBound {
    fn from(node: BnbNode) -> Self {
        Self { node }
    }
}
