use std::cmp::Ordering;
use std::ops::{
    Add, AddAssign,
    Div, DivAssign,
    Mul, MulAssign,
    Sub, SubAssign,
    Neg,
};



/// Scalar types are allowed (but not required) to have NaN values; if they do,
/// they should follow the same comparison semantics as floats.
/// 
/// Non-NaN values (as determined by `is_nan()`) MUST be totally ordered.
pub trait Scalar:
    Add<Output=Self> + Div<Output=Self> + Mul<Output=Self> + Sub<Output=Self>
    + AddAssign + DivAssign + MulAssign + SubAssign + Neg<Output=Self>
    + PartialEq + PartialOrd  // note f32 and f64 do not implement Eq or Ord due to NaN values
    + Copy + Sized + 'static {

    // if a struct type, note Rust supports immutable global const struct values
    const ONE: Self;  
    const ZERO: Self;

    /// Tests whether this is a NaN value.
    fn is_nan(self) -> bool;

    /// Tests (exactly) for a nonzero value.
    fn is_nonzero(self) -> bool;

    /// Tests (exactly) for a zero value.
    fn is_zero(self) -> bool;

    /// Returns a wrapper signifying that this value is non-NaN, or returns
    /// None if it is a NaN value.
    fn non_nan(self) -> Option<NonNan<Self>> {
        if self.is_nan() {
            None
        }
        else {
            Some(NonNan { value: self })
        }
    }

}



/// Wrapper that signifies this Scalar value is non-NaN.
#[derive(Clone, Copy)]
pub struct NonNan<S: Scalar> {
    value: S,
}

impl<S: Scalar> NonNan<S> {

    #[inline]
    pub fn value(&self) -> S {
        self.value
    }

}

impl<S: Scalar> Eq for NonNan<S> {}

impl<S: Scalar> Ord for NonNan<S> {

    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.value.partial_cmp(&other.value)
            .unwrap_or_else(|| panic!("Comparison failed for NonNan<{}>", std::any::type_name::<S>()))
    }

}

impl<S: Scalar> PartialEq for NonNan<S> {

    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }

}

impl<S: Scalar> PartialOrd for NonNan<S> {

    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.value.partial_cmp(&other.value)
    }
    
}
