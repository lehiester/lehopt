// currently only implemented for 32-bit and 64-bit architectures, should be
// easy to extend if needed (do_not_compile is not a real config attribute but
// it's lazy-evaluated so it only throws an error when intended)
#[cfg_attr(not(any(target_pointer_width = "32", target_pointer_width = "64")), do_not_compile)]



use std::cmp::{Eq, Ord};
use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Sub, SubAssign};



/// Int trait: for finite, two's-complement (or unsigned), pass-by-value integer types.
pub trait Int: Copy + Display + Debug + Eq + Ord + Sized + 'static
    + Add<Output=Self> + Sub<Output=Self>
    + AddAssign + SubAssign
    + TryFrom<usize,Error:Debug> + TryInto<usize,Error:Debug> {
    
    const BITS: u32;
    const MAX: Self;
    const ONE: Self;
    const SIGNED: bool;
    const ZERO: Self;

    /// Whether the cardinality should be checked for sets of elements that are
    /// intended to *all* be processed; primarily intended for array storage.
    /// This is determined primarily with single-machine use in mind, or
    /// perhaps small-scale clusters (say, up to single-digit-teraflops scale).
    /// This choice may not be appropriate for petaflops-scale systems, but
    /// probably nothing else in this package will be, either.
    const CHECK_CARDINALITY: bool;

    fn from_usize_unchecked(value: usize) -> Self;

    fn up_to(self, end_excl: Self) -> IntIterator<Self> {
        IntIterator {
            next: self,
            end_excl: end_excl,
        }
    }
}



// IntIterator for `for i in i1.up_to(i2)` in lieu of `for i in i1..i2` syntax;
// at time of writing (Dec. 2024) the necessary internal traits that would need
// to be added as bounds to `Int` in order for the compiler to allow `i1..i2`
// syntax are not all part of the "stable" API so this syntax cannot be used
// with Int in stable versions of Rust

pub struct IntIterator<I> where I: Int {
    next: I,
    end_excl: I,
}

impl<I> Iterator for IntIterator<I> where I: Int {
    type Item = I;

    fn next(&mut self) -> Option<I> {
        if self.next < self.end_excl {
            let val = self.next;
            self.next += I::ONE;
            Some(val)
        }
        else {
            None
        }
    }
}



// Int trait "implementation" for allowed types

impl Int for i32 {
    const BITS: u32 = i32::BITS;
    const MAX: i32 = i32::MAX;
    const ONE: i32 = 1;
    const SIGNED: bool = true;
    const ZERO: i32 = 0;

    const CHECK_CARDINALITY: bool = true;

    #[inline]
    fn from_usize_unchecked(value: usize) -> Self {
        value as i32
    }
}

impl Int for i64 {
    const BITS: u32 = i64::BITS;
    const MAX: i64 = i64::MAX;
    const ONE: i64 = 1;
    const SIGNED: bool = true;
    const ZERO: i64 = 0;

    // see comments for `usize`
    const CHECK_CARDINALITY: bool = false;

    #[inline]
    fn from_usize_unchecked(value: usize) -> Self {
        // hopefully the compiler makes this a no-op
        value as i64
    }
}

impl Int for usize {
    const BITS: u32 = usize::BITS;
    const MAX: usize = usize::MAX;
    const ONE: usize = 1;
    const SIGNED: bool = false;
    const ZERO: usize = 0;

    // 18 quintillion is too large to have any practical use cases on
    // gigaflops-scale systems; would need petaflops-scale for one FLOP per
    // element to take less than a second (gigaflops-scale would be in the
    // order of CPU-days per operation, and probably much longer actually due
    // to such systems needing to rely on much slower storage for such a large
    // amount of data)
    #[cfg(target_pointer_width = "64")]
    const CHECK_CARDINALITY: bool = false;
    
    #[cfg(target_pointer_width = "32")]
    const CHECK_CARDINALITY: bool = true;

    #[inline]
    fn from_usize_unchecked(value: usize) -> Self {
        // hopefully the compiler makes this a no-op
        value
    }
}



/// Trait for integer types whose nonnegative values fit in `usize`.
pub trait IntNonnegFitsInUsize: Int {
    /// Casts this integer to usize.
    /// Does *not* check for nonnegativity, will wrap around if negative!
    fn to_usize_unchecked(self) -> usize;
}

impl IntNonnegFitsInUsize for usize {
    #[inline]
    fn to_usize_unchecked(self) -> usize {
        // hopefully the compiler makes this a no-op
        self
    }
}

impl IntNonnegFitsInUsize for i32 {
    #[inline]
    fn to_usize_unchecked(self) -> usize {
        self as usize
    }
}

#[cfg(target_pointer_width = "64")]
impl IntNonnegFitsInUsize for i64 {
    #[inline]
    fn to_usize_unchecked(self) -> usize {
        // hopefully the compiler makes this a no-op
        self as usize
    }
}
