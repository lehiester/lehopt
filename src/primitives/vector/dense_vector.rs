use crate::primitives::{int::IntNonnegFitsInUsize, scalar::Scalar};

use std::iter::zip;
use std::marker::PhantomData;
use std::ops::Index;



pub struct DenseVector<S=f64,I=usize> where S: Scalar, I: IntNonnegFitsInUsize {
    values: Vec<S>,
    _i: PhantomData<I>,
}

impl<S,I> DenseVector<S,I> where S: Scalar, I: IntNonnegFitsInUsize {

    // Constructors

    pub fn new(values: Vec<S>) -> Self {
        // validate arguments: `values` do not exceed capacity of int type I
        assert!(values.len() <= I::MAX.to_usize_unchecked());

        // validate arguments: nonzero length
        assert!(values.len() > I::ZERO.to_usize_unchecked());

        Self {
            values,
            _i: PhantomData
        }
    }

    pub fn new_from_sparse_lists(length: I, indices: Vec<I>, values: Vec<S>) -> Self {
        // TODO: validate indices in bounds

        // validate arguments: `indices` and `values` same length
        assert_eq!(indices.len(), values.len());

        let mut dense_values = vec![S::ZERO; length.to_usize_unchecked()];
        let mut assigned = vec![false; length.to_usize_unchecked()];

        for (i, v) in zip(indices, values) {
            let i_usize = i.to_usize_unchecked();
            // TODO FUTURE: graceful error handling?
            if assigned[i_usize] {
                panic!("Duplicate indices passed to DenseVector::new_from_sparse_lists");
            }
            assigned[i_usize] = true;
            dense_values[i_usize] = v;
        }

        Self {
            values: dense_values,
            _i: PhantomData,
        }
    }

    pub fn new_unit(length: I, index: I) -> Self {
        assert!(index >= I::ZERO);
        assert!(index < length);

        let mut values = vec![S::ZERO; length.to_usize_unchecked()];
        values[index.to_usize_unchecked()] = S::ONE;

        Self {
            values,
            _i: PhantomData
        }
    }

    pub fn new_zeros(length: I) -> Self {
        Self {
            values: vec![S::ZERO; length.to_usize_unchecked()],
            _i: PhantomData,
        }
    }

    // Getters

    #[inline]
    pub fn values(&self) -> &Vec<S> {
        &self.values
    }

    // Instance methods

    pub fn append_inplace(&mut self, value: S) {
        assert!(self.values.len() < I::MAX.to_usize_unchecked());
        self.values.push(value);
    }

    pub fn append_rep_inplace<I2>(&mut self, value: S, num_reps: I2)
        where I2: IntNonnegFitsInUsize {

        assert!(num_reps > I2::ZERO);
        let num_reps_usize = num_reps.to_usize_unchecked();

        assert!(self.values.len().checked_add(num_reps_usize).unwrap() < I::MAX.to_usize_unchecked());
        self.values.extend((0..num_reps_usize).map(|_| value));
    }

    pub fn append_zeros_inplace<I2>(&mut self, num_zeros: I2)
        where I2: IntNonnegFitsInUsize {
        
        self.append_rep_inplace(S::ZERO, num_zeros);
    }

    pub fn dot(&self, other: &Self) -> S {
        assert_eq!(self.len(), other.len());
        let mut total = S::ZERO;
        for i in 0..self.len().to_usize_unchecked() {
            total += self.values[i] * other.values[i];
        }
        total
    }

    pub fn extend_inplace(&mut self, mut new_values: Vec<S>) {
        assert!(self.values.len().checked_add(new_values.len()).unwrap() < I::MAX.to_usize_unchecked());
        self.values.append(&mut new_values);
    }

    pub fn len(&self) -> I {
        self.values.len().try_into().unwrap()
    }

    /// Returns the smallest index for ties.
    pub fn min_elem(&self) -> (I, S) {
        let mut min_idx: usize = 0;
        let mut min_value = self.values[0];
        for i in 1..self.values.len() {
            if self.values[i] < min_value {
                min_idx = i;
                min_value = self.values[i];
            }
        }
        (min_idx.try_into().unwrap(), min_value)
    }

    pub fn mul_value_inplace<I2>(&mut self, index: I2, multiplier: S)
        where I2: IntNonnegFitsInUsize {

        self.values[index.to_usize_unchecked()] *= multiplier;
    }

    pub fn subvector(&self, indices: &Vec<I>) -> Self {
        assert!(indices.len() <= I::MAX.to_usize_unchecked());
        
        Self {
            values: indices.iter().map(|i| self.values[i.to_usize_unchecked()]).collect(),
            _i: PhantomData
        }
    }

    pub fn to_vec(self) -> Vec<S> {
        self.values
    }

    pub fn trunc_inplace<I2>(&mut self, new_length: I2)
        where I2: IntNonnegFitsInUsize {
        
        assert!(new_length > I2::ZERO);
        let new_length_usize = new_length.to_usize_unchecked();
        assert!(new_length_usize <= self.values.len());

        if new_length_usize < self.values.len() {
            self.values.truncate(new_length_usize);
        }
    }

}

impl<S,I> Clone for DenseVector<S,I> where S: Scalar, I: IntNonnegFitsInUsize {

    fn clone(&self) -> Self {
        Self {
            values: self.values.clone(),
            _i: PhantomData,
        }
    }

}

impl<S,I,I2> Index<I2> for DenseVector<S,I>
    where S: Scalar, I: IntNonnegFitsInUsize, I2: IntNonnegFitsInUsize {
    
    type Output = S;

    fn index(&self, index: I2) -> &S {
        &self.values[index.to_usize_unchecked()]
    }

}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_dot_000() {
        let v1: DenseVector = DenseVector::new(vec![1.0, 2.0, 3.0, 4.0]);
        let v2: DenseVector = DenseVector::new(vec![2.0, 3.0, 5.0, 7.0]);

        assert_eq!(v1.dot(&v2), 51.0);
        assert_eq!(v2.dot(&v1), 51.0);
    }

    #[test]
    fn test_len_000() {
        let v1 = DenseVector::<f64>::new_zeros(13);
        assert_eq!(v1.len(), 13);

        let v2 = DenseVector::<f64>::new_from_sparse_lists(5, vec![1, 2, 3], vec![3.0, 2.0, 1.0]);
        assert_eq!(v2.len(), 5);
    }

    #[test]
    fn test_new_from_sparse_lists_000() {
        let v1 = DenseVector::<f64>::new_from_sparse_lists(7, vec![], vec![]);
        assert_eq!(v1.values, vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);

        let v2 = DenseVector::<f64>::new_from_sparse_lists(
            7,
            vec![4, 2, 6, 3, 0],
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
        );
        assert_eq!(v2.values, vec![5.0, 0.0, 2.0, 4.0, 1.0, 0.0, 3.0]);
    }

    #[test]
    #[should_panic(expected="Duplicate indices passed to DenseVector::new_from_sparse_lists")]
    fn test_new_from_sparse_lists_001() {
        // duplicate element checking

        let _ = DenseVector::<f64>::new_from_sparse_lists(
            100,
            vec![0, 13, 29, 37, 43, 13, 71],
            vec![0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        );
    }

    #[test]
    fn test_subvector_000() {
        let v: DenseVector = DenseVector::new(vec![2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0]);

        let sv1 = v.subvector(&vec![0, 2, 1, 5]);
        assert_eq!(sv1.values, vec![2.0, 5.0, 3.0, 13.0]);

        let sv2 = v.subvector(&vec![6, 3, 7, 4]);
        assert_eq!(sv2.values, vec![17.0, 7.0, 19.0, 11.0]);
    }

}
