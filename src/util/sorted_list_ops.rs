pub fn sorted_list_setminus_keyed<T,K,R,FK,FR>(list1: &Vec<T>, list2: &Vec<T>, key_func: FK, result_func: FR) -> Vec<R> 
    where K: PartialOrd, FK: Fn(&T) -> K, FR: Fn(&T) -> R {

    if list2.len() == 0 {
        return list1.iter().map(result_func).collect();
    }
    
    let mut result = vec![];

    if list1.len() > 0 {
        let mut i1 = 0;
        let mut i2 = 0;
        let mut list1_end = false;
        let mut list2_end = false;
        let mut i1_incr = false;
        let mut i2_incr = false;
        loop {
            let k1 = key_func(&list1[i1]);
            let k2 = key_func(&list2[i2]);
            
            if k1 < k2 || (k1 > k2 && list2_end) {
                // check if k1 is between previous and current values of k2 (or past end of list2, or before start of list2)
                if list2_end || (i2 == 0 || key_func(&list2[i2 - 1]) < k1) {
                    // found an item in list1 not in list2
                    result.push(result_func(&list1[i1]));
                }
                
                i1_incr = true;
            }
            else if k1 == k2 {
                i1_incr = true;
                i2_incr = true;
            }
            else {  // k1 > k2 (and not list2_end)
                i2_incr = true;
            }

            // check if i1 should be (and can be) incremented
            if i1_incr {
                if i1 < list1.len() - 1 {
                    i1 += 1;
                }
                else {
                    list1_end = true;
                }
            }

            // check if i2 should be (and can be) incremented
            if i2_incr {
                if i2 < list2.len() - 1 {
                    i2 += 1;
                }
                else {
                    list2_end = true;
                }
            }

            if list1_end {
                break;
            }
        }
    }

    result
}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_sorted_list_setminus_keyed_000() {
        let list1a = vec![];
        let list1b = vec![1, 2, 3];
        let result1 = sorted_list_setminus_keyed(&list1a, &list1b, |x| *x, |x| *x);
        assert_eq!(result1.len(), 0);

        let list2a = vec![2, 3, 5, 7];
        let list2b = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let result2 = sorted_list_setminus_keyed(&list2a, &list2b, |x| *x, |x| *x);
        assert_eq!(result2.len(), 0);

        let list3a = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let list3b = vec![2, 3, 5, 7];
        let result3 = sorted_list_setminus_keyed(&list3a, &list3b, |x| *x, |x| *x);
        assert_eq!(result3, vec![1, 4, 6, 8, 9]);

        let list4a = vec![1, 2, 3, 4];
        let list4b = vec![3, 5, 9];
        let result4 = sorted_list_setminus_keyed(&list4a, &list4b, |x| *x, |x| *x);
        assert_eq!(result4, vec![1, 2, 4]);

        let list5a = vec![1, 2, 3];
        let list5b = vec![];
        let result5 = sorted_list_setminus_keyed(&list5a, &list5b, |x| *x, |x| *x);
        assert_eq!(result5, list5a);

        let list6a = vec![1, 2, 3];
        let list6b = vec![-3, -2, -1];
        let result6 = sorted_list_setminus_keyed(&list6a, &list6b, |x| *x, |x| *x);
        assert_eq!(result6, list6a);
    }

}
