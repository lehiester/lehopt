use crate::primitives::int::IntNonnegFitsInUsize;
use crate::primitives::scalar::Scalar;
use crate::primitives::vector::DenseVector;
use crate::util::cosort_two_vecs_inplace;



/// A sparse vector with zero-based indices.
/// 
/// The invariants must be very carefully maintained because CSCMatrix relies
/// on them, and CSCMatrix is used with unsafe code.
/// 
/// Invariants:
/// - 01: positive length (zero length not allowed!)
/// - 02: values and indices are same length
/// - 03: value indices fit within max for the given int type
/// - 04: indices are sorted ascending
/// - 05: indices have no duplicates
/// - 06: indices are between 0 and (length - 1) inclusive
pub struct SparseVector<S=f64,I=usize> where S: Scalar, I: IntNonnegFitsInUsize {
    length: I,
    indices: Vec<I>,
    values: Vec<S>,
}

impl<S,I> SparseVector<S,I> where S: Scalar, I: IntNonnegFitsInUsize {

    // Constructors

    pub fn new(length: I, mut indices: Vec<I>, mut values: Vec<S>) -> Self {
        // TODO: unit tests for invariants

        // validate invariant 01: positive length
        assert!(length > I::ZERO);

        // validate invariant 02: values and indices are same length
        assert_eq!(indices.len(), values.len());

        // assuming at least one element for other invariant checks (satisfied automatically otherwise)
        if indices.len() > 0 {
            // validate invariant 03: value indices fit within max for the given int type
            assert!(values.len() <= I::MAX.to_usize_unchecked());
    
            // coerce to invariant 04: indices are sorted ascending
            cosort_two_vecs_inplace(&mut indices, &mut values, |i, _| *i);
    
            // validate invariant 05: indices have no duplicates
            for i in 1..indices.len() {
                assert!(indices[i] > indices[i - 1]);
            }
    
            // validate invariant 06: indices are between 0 and (length - 1) inclusive
            assert!(indices[0] >= I::ZERO);
            assert!(indices[indices.len() - 1] < length);
        }

        SparseVector {
            length,
            indices,
            values,
        }
    }

    pub fn new_from_dense(dv: &DenseVector<S,I>) -> Self {
        // TODO: invariant checking
        let length: I = dv.len();

        let mut indices: Vec<I> = vec![];
        let mut values: Vec<S> = vec![];

        for i in I::ZERO.up_to(length) {
            if dv[i] != S::ZERO {
                indices.push(i);
                values.push(dv[i]);
            }
        }

        Self {
            length,
            indices,
            values,
        }
    }

    pub fn new_without_zeros(length: I, indices: Vec<I>, values: Vec<S>) -> Self {
        // TODO: invariant checking

        // validate invariant 01: positive length
        assert!(length > I::ZERO);

        assert_eq!(indices.len(), values.len());
        
        let num_zeros = values.iter().filter(|v| **v == S::ZERO).count();

        if num_zeros == 0 {
            return Self::new(length, indices, values);
        }

        let mut new_indices: Vec<I> = vec![];
        let mut new_values: Vec<S> = vec![];

        for i in 0..values.len() {
            if values[i] != S::ZERO {
                new_indices.push(indices[i]);
                new_values.push(values[i]);
            }
        }

        Self {
            length,
            indices: new_indices,
            values: new_values,
        }
    }

    pub fn new_zeros(length: I) -> Self {
        // validate invariant 01: positive length
        assert!(length > I::ZERO);

        // satisfies invariants 02 through 06 trivially
        SparseVector {
            length,
            indices: vec![],
            values: vec![],
        }
    }

    // Getters

    pub fn length(&self) -> I {
        self.length
    }

    pub fn indices(&self) -> &Vec<I> {
        &self.indices
    }

    pub fn values(&self) -> &Vec<S> {
        &self.values
    }

    // Instance methods

    pub fn clear(&mut self) {
        // satisfies invariant 01 (positive length) due to existing `self`;
        // satisfies invariants 02 through 06 trivially
        self.indices = vec![];
        self.values = vec![];
    }

    pub fn element(&self, index: I) -> S {
        // binary search over indices

        if index < I::ZERO {
            panic!("Attempted to use negative index: {}", index);
        }
        else if index >= self.length {
            panic!("Attempted to use (zero-based) index of {} for vector with length {}", index, self.length);
        }
        else if index < self.indices[0] || index  > self.indices[self.indices.len() - 1] {
            return S::ZERO;
        }

        let mut a: usize = 0;
        let mut b: usize = self.indices.len();

        loop {
            let k: usize = (a + b) / 2;
            let i = self.indices[k];
            if i == index {
                return self.values[k];
            }
            else if i < index {
                a = k + 1;
            }
            else {
                b = k;
            }

            if (a == b) && self.indices[a] != index {
                // finished searching, did not find this index; return zero
                return S::ZERO;
            }
        }
    }
    
    pub fn set_value_in_order(&mut self, index: I, value: S) {
        // maintains invariant 01: positive length
        // maintains invariant 02: values and indices are same length

        // validate invariant 03: value indices fit within max for the given int type
        assert!(self.values.len() < I::MAX.to_usize_unchecked());

        // validate invariants 04 and 05: indices are ascending with no duplicates
        if self.indices.len() > 0 {
            let prev_idx = self.indices[self.indices.len() - 1];
            if prev_idx == index {
                panic!("Attempted to add duplicate index to sparse vector");
            }
            else if prev_idx > index {
                panic!("Attempted to call SparseVector.append_in_order with out-of-order indices");
            }
        }

        // validate invariant 06: indices are between 0 and (length - 1) inclusive
        assert!(index >= I::ZERO);
        assert!(index < self.length);

        // append new element
        self.indices.push(index);
        self.values.push(value);
    }

    pub fn to_vecs(self) -> (Vec<I>, Vec<S>) {
        (self.indices, self.values)
    }

}



impl<S,I> Clone for SparseVector<S,I> where S: Scalar, I: IntNonnegFitsInUsize {

    fn clone(&self) -> Self {
        Self {
            length: self.length,
            indices: self.indices.clone(),
            values: self.values.clone(),
        }
    }

}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_element_000() {
        let indices1 = vec![1, 2, 999];
        let values1 = vec![7.0, 8.0, 9.0];
        let v1: SparseVector = SparseVector::new(1000, indices1, values1);
        assert_eq!(v1.length, 1000);
        assert_eq!(v1.element(2), 8.0);
        assert_eq!(v1.element(3), 0.0);
        assert_eq!(v1.element(1), 7.0);
        assert_eq!(v1.element(0), 0.0);
        assert_eq!(v1.element(998), 0.0);
        assert_eq!(v1.element(999), 9.0);

        let indices2 = vec![5, 8, 9, 10, 13, 16, 17, 20, 22, 23];
        let values2 = vec![11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0];
        let v2: SparseVector = SparseVector::new(24, indices2, values2);
        assert_eq!(v2.length, 24);
        let dense_values = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 11.0, 0.0, 0.0, 12.0, 13.0, 14.0, 0.0,
            0.0, 15.0, 0.0, 0.0, 16.0, 17.0, 0.0, 0.0, 18.0, 0.0, 19.0, 20.0
        ];
        for (i, v) in dense_values.iter().enumerate() {
            assert_eq!(v2.element(i), *v);
        }
    }

    #[test]
    #[should_panic]
    fn test_element_001() {
        let indices = vec![5, 8, 9, 10, 13, 16, 17, 21, 23, 24];
        let values = vec![11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0];
        let v: SparseVector = SparseVector::new(25, indices, values);
        let _ = v.element(25);
    }

    #[test]
    fn test_new_000() {
        let v0: SparseVector = SparseVector::new(10, vec![], vec![]);
        assert_eq!(v0.length, 10);
        assert_eq!(v0.indices, vec![]);
        assert_eq!(v0.values, vec![]);

        let v1: SparseVector = SparseVector::new(10, vec![0], vec![1.0]);
        assert_eq!(v1.length, 10);
        assert_eq!(v1.indices, vec![0]);
        assert_eq!(v1.values, vec![1.0]);

        let v2: SparseVector = SparseVector::new(10, vec![3, 9], vec![2.0, 3.0]);
        assert_eq!(v2.length, 10);
        assert_eq!(v2.indices, vec![3, 9]);
        assert_eq!(v2.values, vec![2.0, 3.0]);
    }

}
