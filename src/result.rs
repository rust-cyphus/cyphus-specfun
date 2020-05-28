use log::warn;
use num::Num;

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
    /// Some sort of Failure occured
    Failure = 10,
}

impl Default for SpecFunCode {
    fn default() -> Self {
        SpecFunCode::Success
    }
}

/// Result structure from a special function evaluation, contianing the value,
/// error and code
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Default)]
pub struct SpecFunResult<T>
where
    T: Num,
{
    pub val: T,
    pub err: T,
    pub code: SpecFunCode,
}

#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Default)]
pub(crate) struct SpecFunResultE10<T>
where
    T: Num,
{
    pub val: T,
    pub err: T,
    pub code: SpecFunCode,
    pub e10: i32,
}

impl<T: std::fmt::Debug + Num> SpecFunResult<T> {
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
            SpecFunCode::Failure => warn!(
                "Failure: Unknown failure occured in {:?} with args {:?}",
                fname, vars
            ),
        }
    }
}

impl<T: std::fmt::Debug + Num> SpecFunResultE10<T> {
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
            SpecFunCode::Failure => warn!(
                "Failure: Unknown failure occured in {:?} with args {:?}",
                fname, vars
            ),
        }
    }
}

/// Convert SpecFunResultE10 object into a SpecFunResult object.
pub(crate) fn result_smash_e(re: &mut SpecFunResultE10<f64>) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    result.code = re.code.clone();
    if re.e10 == 0 {
        result.val = re.val;
        result.err = re.err;
        result
    } else {
        let av = re.val;
        let ae = re.err;

        if (crate::consts::SQRT_DBL_MIN < av)
            && (av < crate::consts::SQRT_DBL_MAX)
            && (crate::consts::SQRT_DBL_MIN < ae)
            && (ae < crate::consts::SQRT_DBL_MAX)
            && (0.49 * crate::consts::LN_DBL_MIN < re.e10 as f64)
            && ((re.e10 as f64) < 0.49 * crate::consts::LN_DBL_MAX)
        {
            let scale = (re.e10 as f64 * std::f64::consts::LN_10).exp();
            result.val = re.val * scale;
            result.err = re.err * scale;
            result
        } else {
            crate::exp::core::exp_mult_err_e(
                re.e10 as f64 * std::f64::consts::LN_10,
                0.0,
                re.val,
                re.err,
            )
        }
    }
}
