use super::scalar_trait::Scalar;



impl Scalar for f64 {

    const ONE: Self = 1.0f64;
    const ZERO: Self = 0.0f64;

    #[inline]
    fn is_nan(self) -> bool {
        self.is_nan()
    }

    #[inline]
    fn is_nonzero(self) -> bool {
        // positive zero and negative zero must compare equal per IEE 754:
        // "Comparisons shall ignore the sign of zero (so +0 = -0)." (2019)
        // (also ensured via `test_is_nonzero` below)
        self != 0.0
    }

    #[inline]
    fn is_zero(self) -> bool {
        // see comment in `is_nonzero` regarding negative zero
        self == 0.0
    }

}



#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_is_nonzero_000() {
        assert_eq!(0.0f64, -0.0f64);
        assert!(!(-0.0f64).is_nonzero());
    }

}
