use log::warn;

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
    /// Failed because too many iterations
    MaxIterErr = 8,
    /// Sanity check failed... should never happen
    SanityCheckErr = 9,
}

/// Result structure from a special function evaluation, contianing the value,
/// error and code
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct SpecFunResult<T> {
    pub val: T,
    pub err: T,
    pub code: SpecFunCode,
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct SpecFunResultE10<T> {
    pub val: T,
    pub err: T,
    pub code: SpecFunCode,
    pub e10: i32,
}

impl<T: std::fmt::Debug> SpecFunResult<T> {
    /// Generate a warning message for a function call with a given set of
    /// arguments and a SpecFunResult
    pub fn issue_warning(&self, fname: &str, vars: &[T]) {
        match &self.code {
            SpecFunCode::Success => (),
            SpecFunCode::DomainErr => warn!(
                "SpecFunDomainErr: Domain error occured in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::RangeErr => warn!(
                "SpecFunRangeErr: Range error occured in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::ZeroDivErr => warn!(
                "SpecFunZeroDivErr: Division by zero in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::UnderflowErr => warn!(
                "SpecFunUnderflow: Underflow occured in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::OverflowErr => warn!(
                "SpecFunOverflow: Overflow occured in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::AccLossErr => warn!(
                "SpecFunAccLossErr: Loss of accuracy in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::RoundoffErr => warn!(
                "SpecFunRoundoffErr: Roundoff encountered in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::MaxIterErr => warn!(
                "SpecFunMaxIterErr: Maximum number of iterations excceded in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::SanityCheckErr => warn!(
                "SanityCheckErr: Sanity check failed in {:?} with args {:?}",
                fname, vars
            ),
        }
    }
}

impl<T: std::fmt::Debug> SpecFunResultE10<T> {
    /// Generate a warning message for a function call with a given set of
    /// arguments and a SpecFunResult
    pub fn issue_warning(&self, fname: &str, vars: &[T]) {
        match &self.code {
            SpecFunCode::Success => (),
            SpecFunCode::DomainErr => warn!(
                "SpecFunDomainErr: Domain error occured in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::RangeErr => warn!(
                "SpecFunRangeErr: Range error occured in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::ZeroDivErr => warn!(
                "SpecFunZeroDivErr: Division by zero in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::UnderflowErr => warn!(
                "SpecFunUnderflow: Underflow occured in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::OverflowErr => warn!(
                "SpecFunOverflow: Overflow occured in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::AccLossErr => warn!(
                "SpecFunAccLossErr: Loss of accuracy in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::RoundoffErr => warn!(
                "SpecFunRoundoffErr: Roundoff encountered in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::MaxIterErr => warn!(
                "SpecFunMaxIterErr: Maximum number of iterations excceded in {:?} with args {:?}",
                fname, vars
            ),
            SpecFunCode::SanityCheckErr => warn!(
                "SanityCheckErr: Sanity check failed in {:?} with args {:?}",
                fname, vars
            ),
        }
    }
}
