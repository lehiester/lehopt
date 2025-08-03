/// Simultaneously sorts two Vecs, applying the same permutation to both.
/// The Vecs must be of the same length.
/// 
/// Uses heap sort, which is not a "stable" sort, i.e. it does not preserve the
/// order of elements that compare equal, a property which is often used for
/// multi-level sorting. Multi-level sorting can easily be implemented here by
/// having `key_func` return a tuple.
pub fn cosort_two_vecs_inplace<T1,T2,K,F>(vec1: &mut Vec<T1>, vec2: &mut Vec<T2>, key_func: F)
    where K: Ord, F: Fn(&T1,&T2) -> K {
    
    // TODO FUTURE: explore whether other inplace sorting methods have better cache performance

    if vec1.len() != vec2.len() {
        panic!("Function `cosort_two_vecs_inplace` requires same length, got {} and {}", vec1.len(), vec2.len());
    }
    let length = vec1.len();

    // turn into max heap in place
    for i in 1..length {
        // insert target element at index `i` into the heap, sift up
        let mut t = i;  // current index
        let t_key = key_func(&vec1[t], &vec2[t]);
        while t != 0 {
            let p = (t + 1) / 2 - 1;
            let p_key = key_func(&vec1[p], &vec2[p]);

            if t_key > p_key {
                // swap current element with parent
                vec1.swap(t, p);
                vec2.swap(t, p);
                t = p;
            }
            else {
                // done sifting up
                break;
            }
        }
    }

    // pop max element off and place after heap one by one
    for i in (1..length).rev() {
        // swap root (max) with target element at the end of the heap range
        vec1.swap(0, i);
        vec2.swap(0, i);

        // indices `i` and beyond are now outside of the heap
        let heap_length = i;

        // sift target element down
        let mut t = 0;
        let t_key = key_func(&vec1[t], &vec2[t]);
        loop {
            let cl = 2*t + 1;
            let cr = 2*t + 2;
            let s;
            
            // case: has both children
            if cr < heap_length {
                let cl_key = key_func(&vec1[cl], &vec2[cl]);
                let cr_key = key_func(&vec1[cr], &vec2[cr]);
                if cl_key > t_key {
                    if cl_key >= cr_key {
                        s = cl;  // swap with left child
                    }
                    else {
                        s = cr;  // swap with right child (greater transitively)
                    }
                }
                else if cr_key > t_key {
                    s = cr;  // swap with right child
                }
                else {
                    break;  // target is greater than both children, done sifting
                }
            }

            // case: has only left child
            else if cl < heap_length {
                let cl_key = key_func(&vec1[cl], &vec2[cl]);
                if cl_key > t_key {
                    s = cl;  // swap with left child
                }
                else {
                    break;  // target is greater than only child, done sifting
                }
            }

            // case: no children
            else {
                break;  // definitely done sifting down once there is no more "down"
            }

            // swap target element with the chosen child
            vec1.swap(t, s);
            vec2.swap(t, s);
            t = s;
        }
    }
}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_cosort_two_vecs_inplace_000() {
        let mut vec1 = vec![7, 5, 3, 2, 19, 17, 13, 11];
        let mut vec2 = vec![1, 2, 3, 4, 5, 6, 7, 8];

        cosort_two_vecs_inplace(&mut vec1, &mut vec2, |i1, _| *i1);
        assert_eq!(vec1, vec![2, 3, 5, 7, 11, 13, 17, 19]);
        assert_eq!(vec2, vec![4, 3, 2, 1, 8, 7, 6, 5]);

        cosort_two_vecs_inplace(&mut vec1, &mut vec2, |_, i2| -(*i2));
        assert_eq!(vec1, vec![11, 13, 17, 19, 2, 3, 5, 7]);
        assert_eq!(vec2, vec![8, 7, 6, 5, 4, 3, 2, 1]);
    }

    #[test]
    fn test_cosort_two_vecs_inplace_001() {
        // test for proper handling of boundary conditions at smallest sizes

        // length 0
        let mut vec1_0: Vec<i32> = vec![];
        let mut vec2_0: Vec<usize> = vec![];
        cosort_two_vecs_inplace(&mut vec1_0, &mut vec2_0, |x1, _| *x1);
        cosort_two_vecs_inplace(&mut vec1_0, &mut vec2_0, |_, x2| *x2);
        assert_eq!(vec1_0.len(), 0);
        assert_eq!(vec2_0.len(), 0);

        // length 1
        let mut vec1_1: Vec<u16> = vec![65535];
        let mut vec2_1: Vec<i64> = vec![-13];
        cosort_two_vecs_inplace(&mut vec1_1, &mut vec2_1, |x1, _| *x1);
        cosort_two_vecs_inplace(&mut vec1_1, &mut vec2_1, |_, x2| *x2);
        assert_eq!(vec1_1, vec![65535]);
        assert_eq!(vec2_1, vec![-13]);

        // length 2
        let mut vec1_2: Vec<i32> = vec![54321, -654321];
        let mut vec2_2: Vec<u8> = vec![255, 0];
        cosort_two_vecs_inplace(&mut vec1_2, &mut vec2_2, |x1, _| *x1);
        assert_eq!(vec1_2, vec![-654321, 54321]);
        assert_eq!(vec2_2, vec![0, 255]);
        cosort_two_vecs_inplace(&mut vec1_2, &mut vec2_2, |_, x2| -(*x2 as i16));
        assert_eq!(vec1_2, vec![54321, -654321]);
        assert_eq!(vec2_2, vec![255, 0]);
    }

}
