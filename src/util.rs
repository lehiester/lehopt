mod avl_tree;
mod sort;
mod sorted_list_ops;

pub use avl_tree::{
    AvlTree,
    AvlTreeItem,
};
pub use sort::cosort_two_vecs_inplace;
pub use sorted_list_ops::sorted_list_setminus_keyed;

use crate::primitives::int::IntNonnegFitsInUsize;



pub fn indices_complement<I,I2>(indices: &Vec<I>, length: I2) -> Vec<I>
    where I: IntNonnegFitsInUsize, I2: IntNonnegFitsInUsize {

    let mut is_complement = vec![true; length.to_usize_unchecked()];
    for i in indices {
        assert!(*i >= I::ZERO);
        let i_usize = i.to_usize_unchecked();
        assert!(i_usize <= I2::MAX.to_usize_unchecked());
        is_complement[i_usize] = false;
    }

    let mut complement = vec![];
    for i in 0..is_complement.len(){
        if is_complement[i] {
            complement.push(i.try_into().unwrap());
        }
    }

    complement
}



pub fn pos_int_is_power_of<I,I2>(int: I, base: I2) -> bool
    where I: IntNonnegFitsInUsize, I2: IntNonnegFitsInUsize {
    
    assert!(base > I2::ZERO);
    assert!(int > I::ZERO);
    if int == I::ONE {
        return true;
    }

    let int_usize = int.to_usize_unchecked();
    let base_usize = base.to_usize_unchecked();

    if int_usize < base_usize {
        return false;
    }

    let mut x = int_usize;
    while x > base_usize {
        if x % base_usize != 0 {
            return false;
        }
        x /= base_usize;
    }

    x == base_usize
}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_indices_complement_000() {
        let v1 = vec![1, 2, 3];
        assert_eq!(indices_complement(&v1, 4), vec![0]);
        assert_eq!(indices_complement(&v1, 5), vec![0, 4]);
        assert_eq!(indices_complement(&v1, 6), vec![0, 4, 5]);

        let v2 = vec![0, 2, 3];
        assert_eq!(indices_complement(&v2, 4), vec![1]);
        assert_eq!(indices_complement(&v2, 5), vec![1, 4]);
        assert_eq!(indices_complement(&v2, 6), vec![1, 4, 5]);

        let v3 = vec![1, 7, 3, 0, 4];
        assert_eq!(indices_complement(&v3, 8), vec![2, 5, 6]);
        assert_eq!(indices_complement(&v3, 9), vec![2, 5, 6, 8]);
    }

    #[test]
    #[should_panic]
    fn test_indices_complement_001() {
        let v = vec![3, -1, 5];
        let _ = indices_complement(&v, 8);
    }

    #[test]
    #[should_panic]
    fn test_indices_complement_002() {
        let v = vec![3, 8, 5];
        let _ = indices_complement(&v, 8);
    }

}
