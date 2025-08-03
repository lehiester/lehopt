use crate::primitives::matrix::CSCMatrix;
use crate::primitives::vector::DenseVector;

use std::ffi::c_void;



type CSInt = i64;

// `#[repr(C)]` - see https://doc.rust-lang.org/nomicon/other-reprs.html
#[repr(C)]
#[allow(non_snake_case)]
struct CSNumFactor {
    L: *const CSMatrix,
    U: *const CSMatrix,
    pinv: *const CSInt,
    B: *const f64,
}

#[repr(C)]
struct CSMatrix {
    nzmax: CSInt,
    m: CSInt,
    n: CSInt,
    p: *const CSInt,
    i: *const CSInt,
    x: *const f64,
    nz: CSInt,
}

#[repr(C)]
struct CSSymAnalysis {
    pinv: *const CSInt,
    q: *const CSInt,
    parent: *const CSInt,
    cp: *const CSInt,
    leftmost: *const CSInt,
    m2: CSInt,
    lnz: f64,
    unz: f64,
}

#[link(name="CSparse")]
extern {
    fn cs_free(p: *mut c_void) -> *mut c_void;
    fn cs_ipvec(p: *const CSInt, b: *const f64, x: *mut f64, n: CSInt) -> CSInt;
    fn cs_lsolve(L: *const CSMatrix, x: *mut f64) -> CSInt;
    fn cs_ltsolve(L: *const CSMatrix, x: *mut f64) -> CSInt;
    fn cs_lu(A: *const CSMatrix, S: *const CSSymAnalysis, tol: f64) -> *const CSNumFactor;
    fn cs_lusol(order: CSInt, A: *const CSMatrix, b: *mut f64, tol: f64) -> CSInt;
    fn cs_malloc(n: CSInt, size: usize) -> *mut c_void;
    fn cs_nfree(N: *const CSNumFactor) -> *const CSNumFactor;
    fn cs_pvec(p: *const CSInt, b: *const f64, x: *mut f64, n: CSInt) -> CSInt;
    fn cs_sfree(S: *const CSSymAnalysis) -> *const CSSymAnalysis;
    fn cs_sqr(order: CSInt, A: *const CSMatrix, qr: CSInt) -> *const CSSymAnalysis;
    fn cs_usolve(U: *const CSMatrix, x: *mut f64) -> CSInt;
    fn cs_utsolve(U: *const CSMatrix, x: *mut f64) -> CSInt;

    static CS_LU_RANK_DEFICIENT: CSNumFactor;
}



// "Safe" wrapper layer
// ----------------------------------------------------------------------------

#[allow(non_snake_case)]
pub struct CSparseLU {
    N: *const CSNumFactor,    // has the LU factors and row permutation
    S: *const CSSymAnalysis,  // has the column permutation, s.q
}

impl Drop for CSparseLU {
    fn drop(&mut self) {
        unsafe {
            cs_nfree(self.N);
            cs_sfree(self.S);
        }
    }
}

trait CSIntExt {
    fn check_retcode(self, origin: &str);
}

impl CSIntExt for CSInt {
    fn check_retcode(self, origin: &str) {
        if self == 0 {
            panic!("CSparse returned error from {}", origin);
        }
    }
}

pub fn csparse_lu(mat_a: &CSCMatrix<f64,usize>) -> CSparseLU {
    assert!(mat_a.num_rows() <= CSInt::MAX as usize);
    assert!(mat_a.num_cols() <= CSInt::MAX as usize);
    assert!(mat_a.col_boundaries().len() <= CSInt::MAX as usize);
    assert!(mat_a.row_indices().len() <= CSInt::MAX as usize);
    assert!(mat_a.values().len() <= CSInt::MAX as usize);

    let mat_a_cs = CSMatrix {
        nzmax: mat_a.values().len().try_into().unwrap(),
        m: mat_a.num_rows() as CSInt,
        n: mat_a.num_cols() as CSInt,
        p: mat_a.col_boundaries().as_ptr() as *const CSInt,
        i: mat_a.row_indices().as_ptr() as *const CSInt,
        x: mat_a.values().as_ptr(),
        nz: -1,
    };

    let sym: *const CSSymAnalysis;
    let num: *const CSNumFactor;
    unsafe {
        // (order = 2 for LU factorization--see comment on implementation of `cs_amd`)
        // (qr = 0 for LU factorization--see body of `cs_sqr`)
        sym = cs_sqr(2, &mat_a_cs, 0);
        if sym.is_null() {
            panic!("CSparse: error in cs_sqr");
        }

        // (tol = 1.0 for the standard "partial pivoting" rule, i.e. no
        // preference for existing order when considering row permutations)
        num = cs_lu(&mat_a_cs, sym, 1.0);
        if num == &CS_LU_RANK_DEFICIENT {
            panic!("CSparse: rank deficiency in cs_lu");
        }
        if num.is_null() {
            // note the `return cs_ndone(..., 0)` results in body of `cs_lu`
            // will also return a null pointer: calling with `ok` = 0
            panic!("CSparseL error in cs_lu");
        }
    }

    CSparseLU {
        N: num,
        S: sym,
    }
}

pub fn csparse_lusol(mat_a: &CSCMatrix<f64,usize>, vec_b: &Vec<f64>) -> Vec<f64> {
    let mut result = vec_b.clone();
    csparse_lusol_inplace(mat_a, &mut result);
    result
}

pub fn csparse_lusol_inplace(mat_a: &CSCMatrix<f64,usize>, vec_b: &mut Vec<f64>) {
    // TODO: thorough unit tests for CSCMatrix for any invariants that this relies on
    // TODO: unit tests for panics for all invalid stuff

    // ensure RHS is same height as coefficient matrix (necessary for memory safety!)
    assert_eq!(mat_a.num_rows(), vec_b.len().try_into().unwrap());

    // invariants of CSCMatrix relied on (not checked here):
    // - poitive row and column counts (m and n)
    // - number of columns is at least 1 less than the max forCSInt
    // - column boundaries array is (num_cols + 1) long
    // - values and row indices are same length
    // - last "column boundary" is length of values
    // - value indices fit within CSInt range
    // - all indices in `col_boundaries` are between 0 and values-length, inclusive
    // - all indices in `row_indices` are between 0 and (num_rows - 1), inclusive
    // - column boundaries are nondecreasing

    assert!(mat_a.num_rows() <= CSInt::MAX as usize);
    assert!(mat_a.num_cols() <= CSInt::MAX as usize);
    assert!(mat_a.col_boundaries().len() <= CSInt::MAX as usize);
    assert!(mat_a.row_indices().len() <= CSInt::MAX as usize);
    assert!(mat_a.values().len() <= CSInt::MAX as usize);

    let mat_a_cs = CSMatrix {
        nzmax: mat_a.values().len().try_into().unwrap(),
        m: mat_a.num_rows() as CSInt,
        n: mat_a.num_cols() as CSInt,
        p: mat_a.col_boundaries().as_ptr() as *const CSInt,
        i: mat_a.row_indices().as_ptr() as *const CSInt,
        x: mat_a.values().as_ptr(),
        nz: -1,
    };

    unsafe {
        // TODO: how to handle checking for/reporting singularity?
        // overwrites the contents of `vec_b`!
        // (order = 2 for LU factorization--see comment on implementation of `cs_amd`)
        cs_lusol(2, &mat_a_cs, vec_b.as_mut_ptr(), 1.0).check_retcode("cs_lusol");
    }
}

pub fn csparse_solve_with_lu(lu: &CSparseLU, rhs: DenseVector<f64>) -> DenseVector<f64> {
    let mut buffer = rhs.to_vec();
    csparse_solve_with_lu_buf(lu, &mut buffer);
    DenseVector::new(buffer)
}

/// Solves (Ax = b) for x, using the LU factorization `lu` already computed
/// for matrix A via `csparse_lu`. Overwrites `rhs` with the solution.
/// 
/// More precisely, solves (P'LUQ'x = b), where P and Q are the row and column
/// permutation matrices applied to A during the factorization, and ' denotes
/// the transpose.
pub fn csparse_solve_with_lu_buf(lu: &CSparseLU, rhs: &mut Vec<f64>) {
    unsafe {
        let m: CSInt = (*(*lu.N).L).m;
        let n: CSInt = (*(*lu.N).L).n;
        assert_eq!(m, n);

        // validate arguments
        assert_eq!(n as usize, rhs.len());

        let x: *mut f64 = cs_malloc(n, std::mem::size_of::<f64>()) as *mut f64;
        if x.is_null() {
            panic!("Out of memory!");
        }

        // cs_ipvec return code unchecked; only error conditions are if a null
        // pointer is passed for the source or destination array
        cs_ipvec((*lu.N).pinv, rhs.as_ptr(), x, n);
        cs_lsolve((*lu.N).L, x).check_retcode("cs_lsolve");
        cs_usolve((*lu.N).U, x).check_retcode("cs_usolve");
        cs_ipvec((*lu.S).q, x, rhs.as_mut_ptr(), n);

        cs_free(x as *mut c_void);
    }
}

pub fn csparse_solve_with_lu_transpose(lu: &CSparseLU, rhs: DenseVector<f64>) -> DenseVector<f64> {
    let mut buffer = rhs.to_vec();
    csparse_solve_with_lu_transpose_buf(lu, &mut buffer);
    DenseVector::new(buffer)
}

/// Solves the transposed system (A'x = b) for x, using the LU factorization
/// `lu` already computed for original matrix A via `csparse_lu`.
/// Overwrites `rhs` with the solution.
pub fn csparse_solve_with_lu_transpose_buf(lu: &CSparseLU, rhs: &mut Vec<f64>) {
    unsafe {
        let m: CSInt = (*(*lu.N).L).m;
        let n: CSInt = (*(*lu.N).L).n;
        assert_eq!(m, n);

        // validate arguments
        assert_eq!(n as usize, rhs.len());

        let x: *mut f64 = cs_malloc(n, std::mem::size_of::<f64>()) as *mut f64;
        if x.is_null() {
            panic!("Out of memory!");
        }

        // cs_pvec return code unchecked; only error conditions are if a null
        // pointer is passed for the source or destination array
        cs_pvec((*lu.S).q, rhs.as_ptr(), x, n);
        cs_utsolve((*lu.N).U, x).check_retcode("cs_utsolve");
        cs_ltsolve((*lu.N).L, x).check_retcode("cs_ltsolve");
        cs_pvec((*lu.N).pinv, x, rhs.as_mut_ptr(), n);

        cs_free(x as *mut c_void);
    }
}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_solve_csparse_lusol_000() {
        let vec_eq_within_tol = |v1, v2| std::iter::zip(v1, v2).all(|t| f64::abs(t.0 - t.1) < 1e-9);

        let mat_a1: CSCMatrix<f64> = CSCMatrix::new(
            3,
            3,
            vec![0, 1, 2, 3],
            vec![0, 1, 2],
            vec![1.0, 1.0, 1.0]
        );

        let mut vec_b1 = vec![2.0, 3.0, 5.0];
        let mut vec_b2 = vec![17.0, 13.0, 11.0];

        let tgt1 = vec_b1.clone();
        let tgt2 = vec_b2.clone();

        csparse_lusol_inplace(&mat_a1, &mut vec_b1);
        csparse_lusol_inplace(&mat_a1, &mut vec_b2);

        assert!(vec_eq_within_tol(vec_b1, tgt1));
        assert!(vec_eq_within_tol(vec_b2, tgt2));

        let mat_a2: CSCMatrix<f64> = CSCMatrix::new(
            3,
            3,
            vec![0, 1, 2, 3],
            vec![1, 2, 0],
            vec![1.0, 1.0, 1.0]
        );

        let vec_b3 = vec![23.0, 29.0, 31.0];
        let vec_b4 = vec![37.0, 41.0, 43.0];

        let vec_b3_clone = vec_b3.clone();
        let vec_b4_clone = vec_b4.clone();

        let result3 = csparse_lusol(&mat_a2, &vec_b3);
        let result4 = csparse_lusol(&mat_a2, &vec_b4);

        assert_eq!(vec_b3, vec_b3_clone);
        assert_eq!(vec_b4, vec_b4_clone);
        assert!(vec_eq_within_tol(result3, vec![29.0, 31.0, 23.0]));
        assert!(vec_eq_within_tol(result4, vec![41.0, 43.0, 37.0]));
    }

}
