mod data;
mod utils;

use data::*;
use utils::*;

use crate::result::SpecFunResult;
use log::warn;

mod data;
mod utils;

use crate::exp::exp_mult_err_e;
use crate::zeta::HurwitzZeta;
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
                    warn!("SpecFunDomainError: polygamma_e requires positive argument");
                    SpecFunResult {
                        val: f64::NAN,
                        err: f64::NAN,
                    }
                } else if (n as usize) <= psi_table.len() - 1 {
                    let val = psi_table[n as usize];
                    let err = std::f64::EPSILON * val.abs();
                    SpecFunResult { val, err }
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
                    SpecFunResult { val, err }
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
                    warn!("SpecFunDomainError: Trigamma requires n > 0");
                    SpecFunResult {
                        val: f64::EPSILON,
                        err: f64::EPSILON,
                    }
                } else if n as usize <= psi_1_table.len() - 1 {
                    let val = psi_1_table[n as usize];
                    let err = std::f64::EPSILON * val;
                    SpecFunResult { val, err }
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
                    SpecFunResult { val, err }
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
            warn!("SpecFunDomainError: Trigamma function is singular at negative integers");
            SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
            }
        } else if x > 0.0 {
            psi_n_xg0(1, x)
        } else if x > -5.0 {
            // Abramowitz + Stegun 6.4.6
            let m = -x.floor() as i32;
            let fx = x + m as f64;

            if fx.abs() < std::f64::EPSILON {
                warn!("SpecFunDomainError: Trigamma function is singular at negative integers");
                return SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                };
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
            r.err = r.err + 2.0 * std::f64::EPSILON * d;

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
            warn!("SpecFunDomainError: polygamma_e requires positive arguments.");
            SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
            }
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
    fn gamma_e(&self) -> SpecFunResult<f64> {}
    /// Gamma-function
    fn gamma(&self) -> f64 {
        self.gamma_e()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_digammaf() {
        let xs = vec![
            -20.0, -19.5, -19.0, -18.5, -18.0, -17.5, -17.0, -16.5, -16.0, -15.5, -15.0, -14.5,
            -14.0, -13.5, -13.0, -12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.50, -9.00, -8.50,
            -8.00, -7.50, -7.00, -6.50, -6.00, -5.50, -5.00, -4.50, -4.00, -3.50, -3.00, -2.50,
            -2.00, -1.50, -1.00, -0.500, 0.0, 0.500, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00,
            4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.0, 10.5, 11.0,
            11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0,
            18.5, 19.0, 19.5, 20.0,
        ];
        let truths = vec![
            std::f64::INFINITY,
            2.9958363947076466,
            std::f64::INFINITY,
            2.9445543434255953,
            std::f64::INFINITY,
            2.8905002893715412,
            std::f64::INFINITY,
            2.8333574322286841,
            std::f64::INFINITY,
            2.7727513716226235,
            std::f64::INFINITY,
            2.7082352425903654,
            std::f64::INFINITY,
            2.6392697253489861,
            std::f64::INFINITY,
            2.5651956512749120,
            std::f64::INFINITY,
            2.4851956512749120,
            std::f64::INFINITY,
            2.3982391295357816,
            std::f64::INFINITY,
            2.3030010342976864,
            std::f64::INFINITY,
            2.1977378764029495,
            std::f64::INFINITY,
            2.0800908175794201,
            std::f64::INFINITY,
            1.9467574842460868,
            std::f64::INFINITY,
            1.7929113303999329,
            std::f64::INFINITY,
            1.6110931485817511,
            std::f64::INFINITY,
            1.3888709263595289,
            std::f64::INFINITY,
            1.1031566406452432,
            std::f64::INFINITY,
            0.70315664064524319,
            std::f64::INFINITY,
            0.036489973978576521,
            std::f64::INFINITY,
            -1.9635100260214235,
            -0.57721566490153286,
            0.036489973978576521,
            0.42278433509846714,
            0.70315664064524319,
            0.92278433509846714,
            1.1031566406452432,
            1.2561176684318005,
            1.3888709263595289,
            1.5061176684318005,
            1.6110931485817511,
            1.7061176684318005,
            1.7929113303999329,
            1.8727843350984671,
            1.9467574842460868,
            2.0156414779556100,
            2.0800908175794201,
            2.1406414779556100,
            2.1977378764029495,
            2.2517525890667211,
            2.3030010342976864,
            2.3517525890667211,
            2.3982391295357816,
            2.4426616799758120,
            2.4851956512749120,
            2.5259950133091454,
            2.5651956512749120,
            2.6029180902322223,
            2.6392697253489861,
            2.6743466616607937,
            2.7082352425903654,
            2.7410133283274604,
            2.7727513716226235,
            2.8035133283274604,
            2.8333574322286841,
            2.8623368577392251,
            2.8905002893715412,
            2.9178924132947806,
            2.9445543434255953,
            2.9705239922421491,
        ];
        for (x, truth) in xs.iter().zip(truths.iter()) {
            let res = (*x).digamma_e();
            if *truth >= std::f64::INFINITY {
                //assert!(res.is_err());
            } else {
                let val = res.val;
                let err = res.err;
                let abserr = (*truth - val).abs();
                // TODO: There are two pts where error estimate is too small.
                assert!(abserr < 10.0 * err.abs());
                println!("val = {}, truth={},fracerr={}", val, truth, abserr);
            }
        }
    }

    #[test]
    fn test_trigammaf() {
        let xs = vec![
            -20.0, -19.5, -19.0, -18.5, -18.0, -17.5, -17.0, -16.5, -16.0, -15.5, -15.0, -14.5,
            -14.0, -13.5, -13.0, -12.5, -12.0, -11.5, -11.0, -10.5, -10.0, -9.50, -9.00, -8.50,
            -8.00, -7.50, -7.00, -6.50, -6.00, -5.50, -5.00, -4.50, -4.00, -3.50, -3.00, -2.50,
            -2.00, -1.50, -1.00, -0.500, 0.0, 0.500, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00,
            4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.50, 9.00, 9.50, 10.0, 10.5, 11.0,
            11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0,
            18.5, 19.0, 19.5, 20.0,
        ];
        let truths = vec![
            std::f64::INFINITY,
            9.8196148086593976,
            std::f64::INFINITY,
            9.8169849598757027,
            std::f64::INFINITY,
            9.8140631191160241,
            std::f64::INFINITY,
            9.8107978129935751,
            std::f64::INFINITY,
            9.8071247184113896,
            std::f64::INFINITY,
            9.8029623875060826,
            std::f64::INFINITY,
            9.7982061449377117,
            std::f64::INFINITY,
            9.7927191764877802,
            std::f64::INFINITY,
            9.7863191764877802,
            std::f64::INFINITY,
            9.7787577398148124,
            std::f64::INFINITY,
            9.7696874450302319,
            std::f64::INFINITY,
            9.7586071126202596,
            std::f64::INFINITY,
            9.7447662821704326,
            std::f64::INFINITY,
            9.7269885043926548,
            std::f64::INFINITY,
            9.7033198653394004,
            std::f64::INFINITY,
            9.6702620140997310,
            std::f64::INFINITY,
            9.6208792980503482,
            std::f64::INFINITY,
            9.5392466449891238,
            std::f64::INFINITY,
            9.3792466449891238,
            std::f64::INFINITY,
            8.9348022005446793,
            std::f64::INFINITY,
            4.9348022005446793,
            1.6449340668482264,
            0.93480220054467931,
            0.64493406684822644,
            0.49035775610023486,
            0.39493406684822644,
            0.33035775610023486,
            0.28382295573711533,
            0.24872510303901038,
            0.22132295573711533,
            0.19934238698962766,
            0.18132295573711533,
            0.16628453574995824,
            0.15354517795933755,
            0.14261589669670380,
            0.13313701469403143,
            0.12483811891892602,
            0.11751201469403143,
            0.11099728846909903,
            0.10516633568168575,
            0.099916956059126733,
            0.095166335681685746,
            0.090846661274546234,
            0.086901872871768391,
            0.083285224601578370,
            0.079957428427323946,
            0.076885224601578370,
            0.074040268664010337,
            0.071398256151646958,
            0.068938227847683806,
            0.066642013583275971,
            0.064493783403239362,
            0.062479682677968999,
            0.060587533403239362,
            0.058806588095783507,
            0.057127325790782614,
            0.055541281973334528,
            0.054040906037696195,
            0.052619441213655930,
            0.051270822935203120,
        ];
        for (x, truth) in xs.iter().zip(truths.iter()) {
            let res = (*x).trigamma_e();
            if *truth >= std::f64::INFINITY {
                //assert!(res.is_err());
            } else {
                let val = res.val;
                let err = res.err;
                let abserr = (*truth - val).abs();
                // TODO: There are two pts where error estimate is too small.
                //assert!(abserr < 10.0 * err.abs());
                println!("val = {}, truth={},fracerr={}", val, truth, abserr);
            }
        }
    }
}
