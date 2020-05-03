use super::psi::digamma_int_e;
use super::psi::polygamma_e;
use super::psi::trigamma_int_e;

use crate::trig::angle_restrict_symm_e;
use crate::trig::complex_lnsin_e;
use crate::{
    cheb::cheb_eval_e,
    consts::{LNPI, ROOT4_DBL_EPS},
    exp::core::{exp_err_e, exp_mult_err_e},
    result::{SpecFunCode, SpecFunResult},
};

use num::{complex::Complex, Float};

use super::data::*;

/// Compute Log(Gamma(z)) using Lanczos method for complex numbers
fn lngamma_lanczos_complex_e(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut y = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(f64::EPSILON, f64::EPSILON),
        code: SpecFunCode::Success,
    };

    let zz = Complex::new(z.re - 1.0, z.im);
    let mut ag = Complex::new(LANCZOS_7_C[0], 0.0);
    for (k, lnzc) in LANCZOS_7_C.iter().enumerate().skip(1) {
        let r = zz.re + (k as f64);
        let i = zz.im;
        let a = lnzc / (r * r + i * i);
        ag.re += a * r;
        ag.im -= a * i;
    }

    let log1 = (zz + 7.5).ln();
    let logag = ag.ln();

    y.val.re = (zz.re + 0.5) * log1.re - zz.im * log1.im - (zz.re + 7.5)
        + 0.918_938_533_204_672_8_f64
        + logag.re;
    y.val.im = zz.im * log1.re + (zz.re + 0.5) * log1.im - zz.im + logag.im;

    y.err.re = 4.0 * f64::EPSILON * y.val.re.abs();
    y.err.im = 4.0 * f64::EPSILON * y.val.im.abs();

    let res = angle_restrict_symm_e(y.val.im);
    y.val.im = res.val;
    y.err.im += res.err;
    y
}

/// Compute Log(Gamma(x)) using Lanczos method
fn lngamma_lanczos_e(x: f64) -> SpecFunResult<f64> {
    let xx = x - 1.0;
    let ag = LANCZOS_7_C
        .iter()
        .enumerate()
        .skip(1)
        .fold(LANCZOS_7_C[0], |acc, (k, val)| {
            acc + val / (xx + (k as f64))
        });
    let term1 = (xx + 0.5) * ((xx + 7.5) * (-1.0).exp()).ln();
    let term2 = 0.918_938_533_204_672_8_f64 + ag.ln();

    let val = term1 + (term2 - 7.0);
    let err = f64::EPSILON * (2.0 * (term1.abs() + term2.abs() + 7.0) + val.abs());

    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

/// Calculate series for g(eps) = Gamma(eps) eps - 1/(1+eps) - eps / 2, as well as its sign
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
fn lngamma_sgn_0_e(eps: f64) -> (SpecFunResult<f64>, f64) {
    let c1 = -0.07721566490153286061_f64;
    let c2 = -0.01094400467202744461_f64;
    let c3 = 0.09252092391911371098_f64;
    let c4 = -0.01827191316559981266_f64;
    let c5 = 0.01800493109685479790_f64;
    let c6 = -0.00685088537872380685_f64;
    let c7 = 0.00399823955756846603_f64;
    let c8 = -0.00189430621687107802_f64;
    let c9 = 0.00097473237804513221_f64;
    let c10 = -0.00048434392722255893_f64;

    let g = eps
        * c10
            .mul_add(eps, c9)
            .mul_add(eps, c8)
            .mul_add(eps, c7)
            .mul_add(eps, c6)
            .mul_add(eps, c5)
            .mul_add(eps, c4)
            .mul_add(eps, c3)
            .mul_add(eps, c2)
            .mul_add(eps, c1);

    let gee = g + 1.0 / (1.0 + eps) + 0.5 * eps;
    let val = (gee / eps.abs()).ln();
    let err = 4.0 * f64::EPSILON * val.abs();
    (
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        },
        eps.signum(),
    )
}

#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
fn lngamma_sgn_sing_e(n: usize, eps: f64) -> (SpecFunResult<f64>, f64) {
    if eps == 0.0 {
        (
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            },
            0.0,
        )
    } else if n == 1 {
        let c0 = 0.07721566490153286061_f64;
        let c1 = 0.08815966957356030521_f64;
        let c2 = -0.00436125434555340577_f64;
        let c3 = 0.01391065882004640689_f64;
        let c4 = -0.00409427227680839100_f64;
        let c5 = 0.00275661310191541584_f64;
        let c6 = -0.00124162645565305019_f64;
        let c7 = 0.00065267976121802783_f64;
        let c8 = -0.00032205261682710437_f64;
        let c9 = 0.00016229131039545456_f64;
        let g5 = c5 + eps * (c6 + eps * (c7 + eps * (c8 + eps * c9)));
        let g = eps * (c0 + eps * (c1 + eps * (c2 + eps * (c3 + eps * (c4 + eps * g5)))));

        // calculate eps gamma(-1+eps), a negative quantity
        let gam_e = g - 1.0 - 0.5 * eps * (1.0 + 3.0 * eps) / (1.0 - eps * eps);
        let val = (gam_e.abs() * eps.abs().recip()).ln();
        let err = 2.0 * f64::EPSILON * val.abs();
        (
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            },
            if eps > 0.0 { -1.0 } else { 1.0 },
        )
    } else {
        let cs1 = -1.6449340668482264365_f64;
        let cs2 = 0.8117424252833536436_f64;
        let cs3 = -0.1907518241220842137_f64;
        let cs4 = 0.0261478478176548005_f64;
        let cs5 = -0.0023460810354558236_f64;
        let e2 = eps * eps;
        let sin_ser = 1.0 + e2 * (cs1 + e2 * (cs2 + e2 * (cs3 + e2 * (cs4 + e2 * cs5))));

        // Calculate series for ln(gamma(1+n-eps))
        let aeps = eps.abs();

        let c0 = lnfact_e(n);
        let psi0 = digamma_int_e(n + 1);
        let psi1 = trigamma_int_e(n + 1);
        let psi2 = if aeps > 0.00001 {
            polygamma_e(2, n as f64 + 1.0)
        } else {
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            }
        };
        let psi3 = if aeps > 0.0002 {
            polygamma_e(3, n as f64 + 1.0)
        } else {
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            }
        };
        let psi4 = if aeps > 0.001 {
            polygamma_e(4, n as f64 + 1.0)
        } else {
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            }
        };
        let psi5 = if aeps > 0.005 {
            polygamma_e(5, n as f64 + 1.0)
        } else {
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            }
        };
        let psi6 = if aeps > 0.01 {
            polygamma_e(6, n as f64 + 1.0)
        } else {
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            }
        };

        let c1 = psi0.val;
        let c2 = psi1.val / 2.0;
        let c3 = psi2.val / 6.0;
        let c4 = psi3.val / 24.0;
        let c5 = psi4.val / 120.0;
        let c6 = psi5.val / 720.0;
        let c7 = psi6.val / 5040.0;

        let lng_ser = c7
            .mul_add(-eps, c6)
            .mul_add(-eps, c5)
            .mul_add(-eps, c4)
            .mul_add(-eps, c3)
            .mul_add(-eps, c2)
            .mul_add(-eps, c1)
            .mul_add(-eps, c0.val);

        let g = -lng_ser - sin_ser.ln();
        let val = g - eps.abs().ln();
        let err = c0.err + 2.0 * f64::EPSILON * (g.abs() + val.abs());

        let sgn = (if n % 2 == 1 { -1.0 } else { 1.0 }) * (if eps > 0 as f64 { 1.0 } else { -1.0 });

        (
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            },
            sgn,
        )
    }
}

/// Compute log(Gamma(1+eps))/eps using the (2,2) pade
/// approximate plus a correction series
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
#[inline]
fn lngamma_1_pade_e(eps: f64) -> SpecFunResult<f64> {
    let n1: f64 = -1.0017419282349508699871138440;
    let n2: f64 = 1.7364839209922879823280541733;
    let d1: f64 = 1.2433006018858751556055436011;
    let d2: f64 = 5.0456274100274010152489597514;
    let num: f64 = (eps + n1) * (eps + n2);
    let den: f64 = (eps + d1) * (eps + d2);
    let pade: f64 = 2.0816265188662692474880210318 * num * den.recip();

    let c0: f64 = 0.004785324257581753;
    let c1: f64 = -0.01192457083645441;
    let c2: f64 = 0.01931961413960498;
    let c3: f64 = -0.02594027398725020;
    let c4: f64 = 0.03141928755021455;
    let eps5 = eps * eps * eps * eps * eps;
    let corr = eps5 * (c0 + eps * (c1 + eps * (c2 + eps * (c3 + c4 * eps))));

    // Error in case wanted for future use
    // let err = 2.0 * T::epsilon() * eps * (pade + corr);
    let val = eps * (pade + corr);
    let err = 2.0 * f64::EPSILON * val.abs();
    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

/// Compute log(Gamma(2+eps))/eps using the (2,2) pade
/// approximate plus a correction series
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
#[inline]
fn lngamma_2_pade_e(eps: f64) -> SpecFunResult<f64> {
    let n1: f64 = 1.000895834786669227164446568;
    let n2: f64 = 4.209376735287755081642901277;
    let d1: f64 = 2.618851904903217274682578255;
    let d2: f64 = 10.85766559900983515322922936;
    let num: f64 = (eps + n1) * (eps + n2);
    let den: f64 = (eps + d1) * (eps + d2);
    let pade: f64 = 2.85337998765781918463568869 * num / den;
    let c0: f64 = 0.0001139406357036744;
    let c1: f64 = -0.0001365435269792533;
    let c2: f64 = 0.0001067287169183665;
    let c3: f64 = -0.0000693271800931282;
    let c4: f64 = 0.0000407220927867950;
    let corr = c4
        .mul_add(eps, c3)
        .mul_add(eps, c2)
        .mul_add(eps, c1)
        .mul_add(eps, c0)
        .mul_add(eps.powi(5), 0.0);

    let val = eps * (pade + corr);
    let err = 2.0 * f64::EPSILON * val.abs();

    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

fn gammastar_ser_e(x: f64) -> SpecFunResult<f64> {
    let y = 1.0 / (x * x);
    let c0 = 1.0 / 12.0;
    let c1 = -1.0 / 360.0;
    let c2 = 1.0 / 1_260.0;
    let c3 = -1.0 / 1_680.0;
    let c4 = 1.0 / 1_188.0;
    let c5 = -691.0 / 360_360.0;
    let c6 = 1.0 / 156.0;
    let c7 = -3617.0 / 122_400.0;
    let ser = c0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * c7))))));

    let val = (ser / x).exp();
    let err = 2.0 * f64::EPSILON * val * (ser / x).max(1.0);

    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

/// Compute Gamma(x) for x >= 1/2
fn gamma_x_gt_half_e(x: f64) -> SpecFunResult<f64> {
    if x == 0.5 {
        // Error term
        let val = 1.77245385090551602729817_f64;
        let err = f64::EPSILON * val;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x <= (FACT_TABLE.len() as f64) && x == x.floor() {
        let n = x.floor() as usize;
        let val = FACT_TABLE[n - 1];
        let err = val * f64::EPSILON;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if (x - 1.0).abs() < 0.01 {
        let eps = x - 1.0;
        let c1: f64 = 0.4227843350984671394;
        let c2: f64 = -0.01094400467202744461;
        let c3: f64 = 0.09252092391911371098;
        let c4: f64 = -0.018271913165599812664;
        let c5: f64 = 0.018004931096854797895;
        let c6: f64 = -0.006850885378723806846;
        let c7: f64 = 0.003998239557568466030;
        let err = f64::EPSILON;
        let val = c7
            .mul_add(eps, c6)
            .mul_add(eps, c5)
            .mul_add(eps, c4)
            .mul_add(eps, c3)
            .mul_add(eps, c2)
            .mul_add(eps, c1)
            .mul_add(eps, x.recip());

        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if (x - 2.0).abs() < 0.01 {
        let eps = x - 2.0;
        let c1 = 0.4227843350984671394_f64;
        let c2 = 0.4118403304264396948_f64;
        let c3 = 0.08157691924708626638_f64;
        let c4 = 0.07424901075351389832_f64;
        let c5 = -0.0002669820687450147683_f64;
        let c6 = 0.011154045718130991049_f64;
        let c7 = -0.002852645821155340816_f64;
        let c8 = 0.0021039333406973880085_f64;

        let err = f64::EPSILON;
        let val = c8
            .mul_add(eps, c7)
            .mul_add(eps, c6)
            .mul_add(eps, c5)
            .mul_add(eps, c4)
            .mul_add(eps, c3)
            .mul_add(eps, c2)
            .mul_add(eps, c1)
            .mul_add(eps, 1.0);

        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 5.0 {
        // Exponentiating the logarithm is fine, as
        // long as the exponential is not so large
        // that it greatly amplifies the error.
        let mut lg = lngamma_lanczos_e(x);
        lg.val = lg.val.exp();
        lg.err = lg.val * (lg.err + 2.0 * f64::EPSILON);
        lg
    } else if x < 10.0 {
        // This is a sticky area. The logarithm
        // is too large and the gammastar series
        // is not good.
        let gamma_8 = 5040.0;
        let t = (2.0 * x - 15.0) / 5.0;
        let c = cheb_eval_e(t, &GAMMA_5_10_DATA, -1.0, 1.0);
        let val = c.val.exp() * gamma_8;
        let err = val * c.err + 2.0 * f64::EPSILON * val;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 171.0 {
        // 171 = max x s.t. gamma(x) is not overflow
        // We do not want to exponentiate the logarithm
        // if x is large because of the inevitable
        // inflation of the error. So we carefully
        // use pow() and exp() with exact quantities.
        let p = x.powf(x * 0.5);
        let e = (-x).exp();
        let q = (p * e) * p;
        let pre = std::f64::consts::SQRT_2 * std::f64::consts::PI.sqrt() * q / x.sqrt();
        let mut gstar = gammastar_ser_e(x);

        gstar.val *= pre;
        gstar.err = (x + 2.5) * f64::EPSILON * gstar.val;

        gstar
    } else {
        let result = SpecFunResult {
            val: std::f64::NAN,
            err: std::f64::NAN,
            code: SpecFunCode::OverflowErr,
        };
        result.issue_warning("psi_x_e", &[x]);
        result
    }
}

// --------------------------------
// ------ External Functions ------
// --------------------------------

/// Compute the factorial of a number, returning an f64
pub fn factorial(n: usize) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };
    if n < 18 {
        result.val = FACT_TABLE[n];
        result.err = 0.0;
    } else if n < FACT_TABLE.len() {
        result.val = FACT_TABLE[n];
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("factorial", &[n as f64]);
    }
    result
}

/// Compute the double factorial of a number, returning an f64
pub fn double_factorial(n: usize) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };
    if n < 26 {
        result.val = DOUBLE_FACT_TABLE[n];
        result.err = 0.0;
    } else if n < FACT_TABLE.len() {
        result.val = DOUBLE_FACT_TABLE[n];
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("double_factorial", &[n as f64]);
    }
    result
}

/// Compute the natural log of n!
pub(crate) fn lnfact_e(n: usize) -> SpecFunResult<f64> {
    if n < FACT_TABLE.len() {
        let val = FACT_TABLE[n].ln();
        SpecFunResult {
            val,
            err: 2.0 * f64::EPSILON * val.abs(),
            code: SpecFunCode::Success,
        }
    } else {
        lngamma_e((n as f64) + 1.0)
    }
}

/// ln(gamma(x)) where x is not a negative integer
///
// Uses real Lanczos method.
// Returns the real part of ln(gamma(x)) when x < 0,
// i.e. ln(|gamma(x)|).
pub(crate) fn lngamma_e(x: f64) -> SpecFunResult<f64> {
    if (x - 1.0).abs() < 0.01 {
        // Note that we must amplify the errors
        // from the Pade evaluations because of
        // the way we must pass the argument, i.e.
        // writing (1-x) is a loss of precision
        // when x is near 1.
        let mut result = lngamma_1_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 1.0).abs());
        result
    } else if (x - 2.0).abs() < 0.01 {
        let mut result = lngamma_2_pade_e(x - 2.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 2.0).abs());
        result
    } else if x >= 0.5 {
        lngamma_lanczos_e(x)
    } else if x == 0.0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("lngamma_e", &[x]);
        result
    } else if x.abs() < 0.02 {
        lngamma_sgn_0_e(x).0
    } else if x > -0.5 / (f64::EPSILON * std::f64::consts::PI) {
        // Try tp extract a fractional part from x
        let z = 1.0 - x;
        let s = (std::f64::consts::PI * z).sin();
        let abss = s.abs();
        if s == 0.0 {
            let result = SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("lngamma_e", &[x]);
            result
        } else if abss < std::f64::consts::PI * 0.015 {
            // x is near a negative integer
            if x < std::i32::MIN as f64 + 2.0 {
                let result = SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                    code: SpecFunCode::RoundoffErr,
                };
                result.issue_warning("lngamma_e", &[x]);
                result
            } else {
                let n = -(x - 0.5) as usize;
                let eps = x + n as f64;
                lngamma_sgn_sing_e(n, eps).0
            }
        } else {
            let mut result = lngamma_lanczos_e(z);
            result.val = LNPI - (abss.ln() + result.val);
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            result
        }
    } else {
        // |x| was too large to extract any fraction part
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::RoundoffErr,
        };
        result.issue_warning("lngamma_e", &[x]);
        result
    }
}

/// Compute ln(gamma(x)) where x is not a negative integer.
///
// Uses real Lanczos method. Determines
// the sign of gamma(x) as well as ln(|gamma(x)|) for x < 0.
// So gamma(x) = sgn * exp(lngamma_sgn_e(x).res).
pub(crate) fn lngamma_sgn_e(x: f64) -> (SpecFunResult<f64>, f64) {
    if (x - 1.0).abs() < 0.01 {
        // Note that we must amplify the errors
        // from the Pade evaluations because of
        // the way we must pass the argument, i.e.
        // writing (1-x) is a loss of precision
        // when x is near 1.
        let mut result = lngamma_1_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 1.0).abs());
        (result, 1.0)
    } else if (x - 2.0).abs() < 0.01 {
        let mut result = lngamma_2_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 2.0).abs());
        (result, 1.0)
    } else if x >= 0.5 {
        (lngamma_lanczos_e(x), 1.0)
    } else if x == 0.0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("lngamma_sgn_e", &[x]);
        (result, 0.0)
    } else if x.abs() < 0.02 {
        lngamma_sgn_0_e(x)
    } else if x > -0.5 / (f64::EPSILON * std::f64::consts::PI) {
        // Try to extract a fractional part from x
        let z = 1.0 - x;
        let s = (std::f64::consts::PI * z).sin();
        let abss = s.abs();
        if s == 0.0 {
            let result = SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("lngamma_sgn_e", &[x]);
            (result, 0.0)
        } else if abss < std::f64::consts::PI * 0.015 {
            // x is near a negative integer
            if x < std::i32::MIN as f64 + 2.0 {
                let result = SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                    code: SpecFunCode::RoundoffErr,
                };
                result.issue_warning("lngamma_sgn_e", &[x]);
                (result, 0.0)
            } else {
                let n = -(x - 0.5) as usize;
                let eps = x + n as f64;
                lngamma_sgn_sing_e(n, eps)
            }
        } else {
            let mut result = lngamma_lanczos_e(z);
            let sgn = if s > 0.0 { 1.0 } else { -1.0 };
            result.val = LNPI - (abss.ln() + result.val);
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            (result, sgn)
        }
    } else {
        // |x| was too large to extract any fraction part
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("lngamma_sgn_e", &[x]);
        (result, 0.0)
    }
}

/// Gamma(x), where x is not a negative integer
///
/// Uses real Lanczos method.
pub(crate) fn gamma_e(x: f64) -> SpecFunResult<f64> {
    if x < 0.5 {
        let rint_x = (x + 0.5).floor() as i32;
        let f_x = x - rint_x as f64;
        let sgn_gamma = if rint_x % 2 == 0 { 1.0 } else { -1.0 };
        let sin_term = sgn_gamma * (std::f64::consts::PI * f_x).sin() / std::f64::consts::PI;

        if sin_term == 0.0 {
            let result = SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("gamma_e", &[x]);
            result
        } else if x > -169.0 {
            let g = gamma_x_gt_half_e(1.0 - x);
            if sin_term.abs() * g.val * f64::MIN_POSITIVE < 1.0 {
                let val = 1.0 / (sin_term * g.val);
                let mut err = (g.err / g.val).abs() * val.abs();
                err += 2.0 * f64::EPSILON * val.abs();
                SpecFunResult {
                    val,
                    err,
                    code: SpecFunCode::Success,
                }
            } else {
                let result = SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                    code: SpecFunCode::UnderflowErr,
                };
                result.issue_warning("gamma_e", &[x]);
                result
            }
        } else {
            // It is hard to control it here.
            // We can only exponentiate the
            // logarithm and eat the loss of
            // precision.
            let (lng, sgn) = lngamma_sgn_e(x);
            exp_mult_err_e(lng.val, lng.err, sgn, 0.0)
        }
    } else {
        gamma_x_gt_half_e(x)
    }
}

/// Regulated Gamma Function, x > 0
///
/// Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
///            = (1 + 1/(12x) + ...), x->Inf
pub(crate) fn gammastar_e(x: f64) -> SpecFunResult<f64> {
    if x <= 0.0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("gammastar_e", &[x]);
        result
    } else if x < 0.5 {
        let lg = lngamma_e(x);
        let lx = x.ln();
        let c = 0.5 * (std::f64::consts::LN_2 + LNPI);
        let lnr_val = lg.val - (x - 0.5) * lx + x - c;
        let lnr_err = lg.err + 2.0 * f64::EPSILON * ((x + 0.5) * lx.abs() + c);
        exp_err_e(lnr_val, lnr_err)
    } else if x < 2.0 {
        let t = 4.0 / 3.0 * (x - 0.5) - 1.0;
        cheb_eval_e(t, &GSTAR_A_DATA, -1.0, 1.0)
    } else if x < 10.0 {
        let t = 0.25 * (x - 2.0) - 1.0;
        let c = cheb_eval_e(t, &GSTAR_B_DATA, -1.0, 1.0);
        let val = c.val / (x * x) + 1.0 + 1.0 / (12.0 * x);
        let mut err = c.err / (x * x);
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 1.0 / ROOT4_DBL_EPS {
        gammastar_ser_e(x)
    } else if x < 1.0 / f64::EPSILON {
        let xi = 1.0 / x;
        let val = 1.0
            + xi / 12.0 * (1.0 + xi / 24.0 * (1.0 - xi * (139.0 / 180.0 + 571.0 / 8640.0 * xi)));
        let err = 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        SpecFunResult {
            val: 1.0,
            err: x.recip(),
            code: SpecFunCode::Success,
        }
    }
}

/// Compute 1 / gamma(x)
///
/// Uses real Lanczos method.
pub(crate) fn gammainv_e(x: f64) -> SpecFunResult<f64> {
    if x <= 0.0 && x == x.floor() {
        SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        }
    } else if x < 0.5 {
        let (lng, sgn) = lngamma_sgn_e(x);
        if lng.val.is_nan() {
            SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::Success,
            }
        } else {
            exp_mult_err_e(-lng.val, lng.err, sgn, 0.0)
        }
    } else {
        let g = gamma_x_gt_half_e(x);
        let val = 1.0 / g.val;
        let mut err = (g.err / g.val).abs() * val.abs();
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

/// Compute ln(gamma(z)) for a complex number z
pub(crate) fn lngamma_complex_e(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut result = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
        code: SpecFunCode::Success,
    };

    if z.re <= 0.5 {
        // Transform to right half plane using reflection;
        // in fact we do a little better by stopping at 1/2.
        let zz = Complex::new(1.0 - z.re, -z.im);

        let lg = lngamma_lanczos_complex_e(zz);
        let lnsin = complex_lnsin_e(std::f64::consts::PI * zz);

        if lnsin.code == SpecFunCode::Success {
            result.val.re = LNPI - lnsin.val.re - lg.val.re;
            result.val.im = -lnsin.val.im - lg.val.im;

            result.err.re = lnsin.err.re + lg.err.re + 2.0 * f64::EPSILON * result.val.re.abs();
            result.err.im = lnsin.err.im + lg.err.im + 2.0 * f64::EPSILON * result.val.im.abs();

            result.val.im = angle_restrict_symm_e(result.val.im).val;
            result
        } else {
            result.val = Complex::new(f64::NAN, f64::NAN);
            result.err = Complex::new(f64::NAN, f64::NAN);
            result.code = SpecFunCode::DomainErr;
            result
        }
    } else {
        lngamma_lanczos_complex_e(z)
    }
}
