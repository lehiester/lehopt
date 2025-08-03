pub mod constr;
pub mod linalg;
pub mod primitives;
pub mod util;



// TEMPORARY EXPERIMENTATION
/*
pub type CSInt = i64;

// `#[repr(C)]` - see https://doc.rust-lang.org/nomicon/other-reprs.html
#[repr(C)]
pub struct CSMatrix {
    pub nzmax: CSInt,
    pub m: CSInt,
    pub n: CSInt,
    pub p: *mut CSInt,
    pub i: *mut CSInt,
    pub x: *mut f64,
    pub nz: CSInt,
}

#[link(name="CSparse")]
extern {
    pub fn cs_lusol(order: CSInt, A: *const CSMatrix, b: *mut f64, tol: f64) -> CSInt;
}



#[cfg(test)]
mod tests {
    
    use super::*;

    #[test]
    fn test_foo() {
        let mut mat_a_p = vec![0, 1, 2, 3];  // IMPORTANT: this have must n + 1 entries!
        let mut mat_a_i = vec![2, 0, 1];
        let mut mat_a_x = vec![1.0, 1.0, 1.0];
        let mat_a = CSMatrix {
            nzmax: 3,
            m: 3,
            n: 3,
            p: mat_a_p.as_mut_ptr(),
            i: mat_a_i.as_mut_ptr(),
            x: mat_a_x.as_mut_ptr(),
            nz: -1,
        };
    
        let mut vec_b = vec![11.0, 13.0, 17.0];
    
        unsafe {
            // TODO: how to handle checking for/reporting singularity?
            // order = 2 for LU factorization (see comment on implementation of `cs_amd`)
            // overwrites the contents of `vec_b`!
            let return_code = cs_lusol(2, &mat_a, vec_b.as_mut_ptr(), 1e-12);
            assert_eq!(return_code, 1);
        }
    
        assert_eq!(vec_b, vec![17.0, 11.0, 13.0]);
    }

}
*/