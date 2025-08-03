// Julia doesn't maintain a factorization of the numerator and denominator,
// just does a GCD calculation on the fly; is this more performant?
// Start with that because it's definitely easier at least



use super::scalar_trait::Scalar;

use std::cmp::Ordering;
use std::ops::{
    Add, AddAssign,
    Div, DivAssign,
    Mul, MulAssign,
    Sub, SubAssign,
    Neg,
};



pub struct Rational {
    // ... TODO
}

impl Rational {
    // ... TODO
}

impl Scalar for Rational {

    const ONE: Self = todo!();
    const ZERO: Self = todo!();

    fn is_nan(self) -> bool {
        todo!();
    }

    fn is_nonzero(self) -> bool {
        todo!();
    }

    fn is_zero(self) -> bool {
        todo!();
    }

}

impl Add for Rational {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        todo!();
    }
}

impl AddAssign for Rational {

    fn add_assign(&mut self, other: Self) {
        todo!();
    }
    
}

impl Copy for Rational {}

impl Clone for Rational {

    fn clone(&self) -> Self {
        todo!();
    }

}

impl Div for Rational {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        todo!();
    }
}

impl DivAssign for Rational {

    fn div_assign(&mut self, rhs: Self) {
        todo!();
    }

}

impl Mul for Rational {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        todo!();
    }
}

impl MulAssign for Rational {
    
    fn mul_assign(&mut self, other: Self) {
        todo!();
    }

}

impl Neg for Rational {
    type Output = Self;

    fn neg(self) -> Self {
        todo!();
    }
}

impl PartialEq for Rational {
    fn eq(&self, other: &Rational) -> bool {
        todo!();
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Rational) -> Option<Ordering> {
        todo!();
    }
}

impl Sub for Rational {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        todo!();
    }
}

impl SubAssign for Rational {

    fn sub_assign(&mut self, rhs: Self) {
        todo!();
    }

}
