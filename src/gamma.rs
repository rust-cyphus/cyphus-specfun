pub(crate) mod data;
pub(crate) mod utils;




use crate::{
    exp::core::exp_mult_err_e,
    result::{SpecFunCode, SpecFunResult},
    zeta::Zeta,
};
use data::*;
use utils::*;

/// Implementation for the digamma function
pub trait DiGamma {
    /// Compute the digamma function of a number along with the error.
    ///
    /// ```
    /// use stercus::specfun::gamma::DiGamma;
    ///
    /// let res = 2.0_f64.digamma_e();
    /// let abserr = (res.val + 0.5772156649).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn digamma_e(&self) -> SpecFunResult<f64>;
    /// Compute the digamma function of a number.
    ///
    /// ```
    /// use stercus::specfun::gamma::DiGamma;
    ///
    /// let res = 2.0_f64.digamma_e();
    /// let abserr = (res.val + 0.5772156649).abs();
    /// assert!(abserr < 1e-10);
    /// assert!(abserr < res.err);
    /// ```
    fn digamma(&self) -> f64;
}

/// Implementation for the trigamma function
pub trait TriGamma {
    /// Compute the trigamma function of a number along with the error.
    ///
    /// ```
    /// use stercus::specfun::gamma::TriGamma;
    ///
    /// let res = 2.0_f64.trigamma_e();
    /// let abserr = (res.val - 1.644934067).abs();
    /// assert!(abserr < 1e-10);
    /// assert!(abserr < res.err);
    /// ```
    fn trigamma_e(&self) -> SpecFunResult<f64>;
    /// Compute the trigamma function of a number.
    ///
    /// ```
    /// use stercus::specfun::gamma::TriGamma;
    ///
    /// let res = 2.0_f64.trigamma();
    /// let abserr = (res - 1.644934067).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn trigamma(&self) -> f64;
}

/// Implementation for the polygamma function
pub trait PolyGamma {
    /// Compute the polygamma function of a number along with the error.
    ///
    /// ```
    /// use stercus::specfun::gamma::PolyGamma;
    ///
    /// let res = 2.0_f64.polygamma_e(5);
    /// let abserr = (res.val - 122.0811674).abs();
    /// assert!(abserr < 1e-10);
    /// assert!(abserr < res.err);
    /// ```
    fn polygamma_e(&self, n: usize) -> SpecFunResult<f64>;
    /// Compute the polygamma function of a number.
    ///
    /// ```
    /// use stercus::specfun::gamma::PolyGamma;
    ///
    /// let res = 2.0_f64.polygamma_e(5);
    /// let abserr = (res.val - 122.0811674).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn polygamma(&self, n: usize) -> f64;
}

/// Implementation for the gamma function
pub trait Gamma {
    /// Gamma-function with error estimate
    fn gamma_e(&self) -> SpecFunResult<f64>;
    /// Gamma-function
    fn gamma(&self) -> f64;
}

pub trait LnFactorial {
    fn lnfact_e(&self) -> SpecFunResult<f64>;
    fn lnfact(&self) -> f64;
}

macro_rules! impl_digamma_int {
    ($T:ty) => {
        impl DiGamma for $T {
            fn digamma_e(&self) -> SpecFunResult<f64> {
                let n = *self;
                if n <= 0 {
                    SpecFunResult {
                        val: f64::NAN,
                        err: f64::NAN,
                        code: SpecFunCode::DomainErr,
                    }
                } else if (n as usize) < PSI_TABLE.len() {
                    let val = PSI_TABLE[n as usize];
                    let err = std::f64::EPSILON * val.abs();
                    SpecFunResult {
                        val,
                        err,
                        code: SpecFunCode::Success,
                    }
                } else {
                    // Abramowitz+Stegun 6.3.18
                    let nf = n as f64;
                    let c2 = -1.0 / 12.0;
                    let c3 = 1.0 / 120.0;
                    let c4 = -1.0 / 252.0;
                    let c5 = 1.0 / 240.0;
                    let ni2 = (1.0 / nf) * (1.0 / nf);
                    let ser = ni2 * (c2 + ni2 * (c3 + ni2 * (c4 + ni2 * c5)));
                    let val = nf.ln() - 0.5 / nf + ser;
                    let mut err =
                        std::f64::EPSILON * (nf.ln().abs() + (0.5 / nf).abs() + ser.abs());
                    err += std::f64::EPSILON * val.abs();
                    SpecFunResult {
                        val,
                        err,
                        code: SpecFunCode::Success,
                    }
                }
            }
            fn digamma(&self) -> f64 {
                self.digamma_e().val
            }
        }
    };
}

impl_digamma_int!(i8);
impl_digamma_int!(i16);
impl_digamma_int!(i32);
impl_digamma_int!(i64);
impl_digamma_int!(i128);
impl_digamma_int!(u8);
impl_digamma_int!(u16);
impl_digamma_int!(u32);
impl_digamma_int!(usize);
impl_digamma_int!(u128);

impl DiGamma for f64 {
    fn digamma_e(&self) -> SpecFunResult<f64> {
        psi_x_e(*self)
    }
    fn digamma(&self) -> f64 {
        self.digamma_e().val
    }
}

macro_rules! impl_trigamma_int {
    ($T:ty) => {
        impl TriGamma for $T {
            fn trigamma_e(&self) -> SpecFunResult<f64> {
                let n = *self;
                if n <= 0 {
                    SpecFunResult {
                        val: f64::EPSILON,
                        err: f64::EPSILON,
                        code: SpecFunCode::DomainErr,
                    }
                } else if (n as usize) < PSI_1_TABLE.len() {
                    let val = PSI_1_TABLE[n as usize];
                    let err = std::f64::EPSILON * val;
                    SpecFunResult {
                        val,
                        err,
                        code: SpecFunCode::Success,
                    }
                } else {
                    // Abramowitz+Stegun 6.4.12
                    // double-precision for n > 100
                    let nf = n as f64;
                    let c0 = -1.0 / 30.0;
                    let c1 = 1.0 / 42.0;
                    let c2 = -1.0 / 30.0;
                    let ni2 = (1.0 / nf) * (1.0 / nf);
                    let ser = ni2 * ni2 * (c0 + ni2 * (c1 + c2 * ni2));
                    let val = (1.0 + 0.5 / nf + 1.0 / (6.0 * nf * nf) + ser) / nf;
                    let err = std::f64::EPSILON * val;
                    SpecFunResult {
                        val,
                        err,
                        code: SpecFunCode::Success,
                    }
                }
            }
            fn trigamma(&self) -> f64 {
                self.digamma_e().val
            }
        }
    };
}

impl_trigamma_int!(i8);
impl_trigamma_int!(i16);
impl_trigamma_int!(i32);
impl_trigamma_int!(i64);
impl_trigamma_int!(i128);
impl_trigamma_int!(u8);
impl_trigamma_int!(u16);
impl_trigamma_int!(u32);
impl_trigamma_int!(usize);
impl_trigamma_int!(u128);

impl TriGamma for f64 {
    fn trigamma_e(&self) -> SpecFunResult<f64> {
        let x = *self;

        if x.abs() < std::f64::EPSILON
            || (x + 1.0).abs() <= std::f64::EPSILON
            || (x + 2.0).abs() <= std::f64::EPSILON
        {
            let result = SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("trigamma_e", &[x]);
            result
        } else if x > 0.0 {
            psi_n_xg0(1, x)
        } else if x > -5.0 {
            // Abramowitz + Stegun 6.4.6
            let m = -x.floor() as i32;
            let fx = x + m as f64;

            if fx.abs() < std::f64::EPSILON {
                let result = SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                    code: SpecFunCode::DomainErr,
                };
                result.issue_warning("trigamma_e", &[x]);
                return result;
            }

            let sum = (0..m)
                .map(|mm| (x + mm as f64).powi(2).recip())
                .sum::<f64>();

            let mut result = psi_n_xg0(1, fx);
            result.val += sum;
            result.err += (m as f64) * std::f64::EPSILON * sum;
            result
        } else {
            let sin_px = (std::f64::consts::PI * x).sin();
            let d = std::f64::consts::PI.powi(2) / (sin_px * sin_px);
            let mut r = psi_n_xg0(1, 1.0 - x);

            r.val = d - r.val;
            r.err += 2.0 * std::f64::EPSILON * d;

            r
        }
    }
    fn trigamma(&self) -> f64 {
        self.trigamma_e().val
    }
}

impl PolyGamma for f64 {
    fn polygamma_e(&self, n: usize) -> SpecFunResult<f64> {
        if n == 0 {
            self.digamma_e()
        } else if n == 1 {
            self.trigamma_e()
        } else if *self <= 0.0 {
            let result = SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("polygamma_e", &[n as f64]);
            result
        } else {
            let hz = (n as f64 + 1.0).hzeta_e(*self);
            let lnnf = n.lnfact_e();
            let mut result = exp_mult_err_e(lnnf.val, lnnf.err, hz.val, hz.err);
            if n % 2 == 0 {
                result.val *= -1.0;
            }
            result
        }
    }
    fn polygamma(&self, n: usize) -> f64 {
        self.polygamma_e(n).val
    }
}

impl LnFactorial for usize {
    fn lnfact_e(&self) -> SpecFunResult<f64> {
        lnfact_int_e(*self)
    }
    fn lnfact(&self) -> f64 {
        lnfact_int_e(*self).val
    }
}

impl Gamma for f64 {
    /// Gamma-function with error estimate
    fn gamma_e(&self) -> SpecFunResult<f64> {
        gamma_e(*self)
    }
    /// Gamma-function
    fn gamma(&self) -> f64 {
        self.gamma_e().val
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{SQRT_DLB_EPS};
    use crate::test_utils::*;
    

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131072.0 * f64::EPSILON;

    #[test]
    fn test_lngamma_sgn_e() {
        test_sf_check_result_and_code(
            utils::lngamma_e(0.7),
            0.260_867_246_531_666_54,
            TOL1,
            SpecFunCode::Success,
        );
    }

    #[test]
    fn test_gamma_e() {
        test_sf_check_result_and_code(
            (1.0 + 1.0 / 4096.0).gamma_e(),
            0.999_859_137_145_940_3,
            TOL0,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            (1.0 + 1.0 / 32.0).gamma_e(),
            0.982_901_099_283_626_9,
            TOL0,
            SpecFunCode::Success,
        );
    }
}
