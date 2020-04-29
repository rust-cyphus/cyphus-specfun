use log::warn;
use num::{Complex, Float};

/// Codes supplying additional information about special function evaluations.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum SpecFunCode {
    /// Evaluation was successful
    Success = 0,
    /// Input domain error, e.g. sqrt(-1)
    DomainErr = 1,
    /// Output range error, e.g. exp(1e100)
    RangeErr = 2,
    /// Tried to divide by zero
    ZeroDivErr = 3,
    /// Floating point underflow
    UnderflowErr = 4,
    /// Floating point overflow
    OverflowErr = 5,
    /// Loss of accuracy
    AccLossErr = 6,
    /// Failed because of roundoff error
    RoundoffErr = 7,
}

/// Result structure from a special function evaluation, contianing the value,
/// error and code
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct SpecFunResult<T> {
    pub val: T,
    pub err: T,
    pub code: SpecFunCode,
}

impl<T: Float> SpecFunResult<T> {
    /// Generate a new SpecFunResult from the error code. This is helpful
    /// if there is an error.
    pub fn new_from_code(code: SpecFunCode) -> SpecFunResult<T> {
        match code {
            Success => SpecFunResult {
                val: T::zero(),
                err: T::zero(),
                code,
            },
            DomainErr => SpecFunResult {
                val: T::nan(),
                err: T::nan(),
                code,
            },
            RangeErr => SpecFunResult {
                val: T::nan(),
                err: T::nan(),
                code,
            },
            ZeroDivErr => SpecFunResult {
                val: T::nan(),
                err: T::nan(),
                code,
            },
            UnderflowErr => SpecFunResult {
                val: T::zero(),
                err: T::zero(),
                code,
            },
            OverflowErr => SpecFunResult {
                val: T::infinity(),
                err: T::infinity(),
                code,
            },
            AccLossErr => SpecFunResult {
                val: T::zero(),
                err: T::zero(),
                code,
            },
            RoundoffErr => SpecFunResult {
                val: T::zero(),
                err: T::zero(),
                code,
            },
        }
    }
}

impl<T: Float> SpecFunResult<Complex<T>> {
    /// Generate a new SpecFunResult from the error code. This is helpful
    /// if there is an error.
    pub fn new_from_code(code: SpecFunCode) -> SpecFunResult<Complex<T>> {
        match code {
            Success => SpecFunResult {
                val: Complex::new(T::zero(), T::zero()),
                err: Complex::new(T::zero(), T::zero()),
                code,
            },
            DomainErr => SpecFunResult {
                val: Complex::new(T::nan(), T::nan()),
                err: Complex::new(T::nan(), T::nan()),
                code,
            },
            RangeErr => SpecFunResult {
                val: Complex::new(T::nan(), T::nan()),
                err: Complex::new(T::nan(), T::nan()),
                code,
            },
            ZeroDivErr => SpecFunResult {
                val: Complex::new(T::nan(), T::nan()),
                err: Complex::new(T::nan(), T::nan()),
                code,
            },
            UnderflowErr => SpecFunResult {
                val: Complex::new(T::zero(), T::zero()),
                err: Complex::new(T::zero(), T::zero()),
                code,
            },
            OverflowErr => SpecFunResult {
                val: Complex::new(T::infinity(), T::infinity()),
                err: Complex::new(T::infinity(), T::infinity()),
                code,
            },
            AccLossErr => SpecFunResult {
                val: Complex::new(T::zero(), T::zero()),
                err: Complex::new(T::zero(), T::zero()),
                code,
            },
            RoundoffErr => SpecFunResult {
                val: Complex::new(T::zero(), T::zero()),
                err: Complex::new(T::zero(), T::zero()),
                code,
            },
        }
    }
}

impl<T: std::fmt::Debug> SpecFunResult<T> {
    /// Generate a warning message for a function call with a given set of
    /// arguments and a SpecFunResult
    pub fn issue_warning(&self, fname: &str, vars: &[T]) {
        match self.code {
            Success => (),
            DomainErr => warn!(
                "SpecFunDomainErr: Domain error occured in {:?} with args {:?}",
                fname, vars
            ),
            RangeErr => warn!(
                "SpecFunRangeErr: Range error occured in {:?} with args {:?}",
                fname, vars
            ),
            ZeroDivErr => warn!(
                "SpecFunZeroDivErr: Division by zero in {:?} with args {:?}",
                fname, vars
            ),
            UnderflowErr => warn!(
                "SpecFunUnderflow: Underflow occured in {:?} with args {:?}",
                fname, vars
            ),
            OverflowErr => warn!(
                "SpecFunOverflow: Overflow occured in {:?} with args {:?}",
                fname, vars
            ),
            AccLossErr => warn!(
                "SpecFunAccLossErr: Loss of accuracy in {:?} with args {:?}",
                fname, vars
            ),
            RoundoffErr => warn!(
                "SpecFunRoundoffErr: Roundoff encountered in {:?} with args {:?}",
                fname, vars
            ),
        }
    }
}
