use super::data::*;
use super::mono::lnfact_e;

use crate::cheb::cheb_eval_e;
use crate::exp::core::exp_mult_err_e;
use crate::result::{SpecFunCode, SpecFunResult};
use crate::zeta::hurwitz::hzeta_e;

use num::complex::Complex;

/// Compute the digamma function for either positive or negative `x`.
fn psi_x_e(x: f64) -> SpecFunResult<f64> {
    let y = x.abs();

    if x == 0.0 || x == -1.0 || x == -2.0 {
        let result = SpecFunResult {
            val: std::f64::NAN,
            err: std::f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("psi_x_e", &[x]);
        result
    } else if y >= 2.0 {
        let t = 8.0 / (y * y) - 1.0;
        let result_c = cheb_eval_e(t, &APSICS_DATA, -1.0, 1.0);
        if x < 0.0 {
            let s = (x * std::f64::consts::PI).sin();
            let c = (x * std::f64::consts::PI).cos();
            if s.abs() < 2.0 * f64::MIN_POSITIVE.sqrt() {
                let result = SpecFunResult {
                    val: std::f64::NAN,
                    err: std::f64::NAN,
                    code: SpecFunCode::DomainErr,
                };
                result.issue_warning("psi_x_e", &[x]);
                result
            } else {
                let val = y.ln() - 0.5 / x + result_c.val - std::f64::consts::PI * c / s;
                let mut err = std::f64::consts::PI * x.abs() * f64::EPSILON / (s * s);
                err += result_c.err;
                err += f64::EPSILON * val.abs();
                SpecFunResult {
                    val,
                    err,
                    code: SpecFunCode::Success,
                }
            }
        } else {
            let val = y.ln() - 0.5 / x + result_c.val;
            let mut err = result_c.err;
            err += f64::EPSILON * val;
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            }
        }
    } else if x < -1.0 {
        let v = x + 2.0;
        let t1 = 1.0 / x;
        let t2 = 1.0 / (x + 1.0);
        let t3 = 1.0 / v;
        let result_c = cheb_eval_e(2.0 * v - 1.0, &PSICS_DATA, -1.0, 1.0);

        let val = -(t1 + t2 + t3) + result_c.val;
        let mut err = f64::EPSILON * (t1.abs() + (x / (t2 * t2)).abs() + (x / (t3 * t3)).abs());
        err += result_c.err;
        err += f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 0.0 {
        let v = x + 1.0;
        let t1 = 1.0 / x;
        let t2 = 1.0 / v;
        let result_c = cheb_eval_e(2.0 * v - 1.0, &PSICS_DATA, -1.0, 1.0);

        let val = -(t1 + t2) + result_c.val;
        let mut err = f64::EPSILON * (t1.abs() + (x / (t2 * t2)).abs());
        err += result_c.err;
        err += f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else if x < 1.0 {
        let t1 = 1.0 / x;
        let result_c = cheb_eval_e(2.0 * x - 1.0, &PSICS_DATA, -1.0, 1.0);

        let val = -t1 + result_c.val;
        let mut err = f64::EPSILON * t1;
        err += result_c.err;
        err += f64::EPSILON * val;
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        let v = x - 1.0;
        cheb_eval_e(2.0 * v - 1.0, &PSICS_DATA, -1.0, 1.0)
    }
}

/// Compute psi(z) for large |z| in the right-half plane
/// ref: [Abramowitz + Stegun, 6.3.18]
fn psi_complex_asymp(z: Complex<f64>) -> Complex<f64> {
    // coefficients in the asymptotic expansion for large z;
    // let w = z^(-2) and write the expression in the form
    //
    //   ln(z) - 1/(2z) - 1/12 w (1 + c1 w + c2 w + c3 w + ... )
    let c1 = -0.1;
    let c2 = 1.0 / 21.0;
    let c3 = -0.05;

    let zinv = z.inv();
    let w = zinv.powi(2);

    // Horner method evaluation of term in parentheses
    let mut sum: Complex<f64>;
    sum = w * (c3 / c2);
    sum += 1.0;
    sum *= c2 / c1;
    sum *= w;
    sum += 1.0;
    sum *= c1;
    sum *= w;
    sum += 1.0;

    // Correction added to log(z)
    let mut cs = sum * w;
    cs *= -1.0 / 12.0;
    cs += zinv * (-0.5);

    cs + z.ln()
}

/// Compute Psi(z) for complex z in the right-half plane
fn psi_complex_rhp(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut n_recurse: usize = 0;
    let mut res = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
        code: SpecFunCode::Success,
    };

    if z.re.abs() < std::f64::EPSILON && z.im.abs() < std::f64::EPSILON {
        res.val = Complex::new(std::f64::NAN, std::f64::NAN);
        res.err = Complex::new(std::f64::NAN, std::f64::NAN);
        res.code = SpecFunCode::DomainErr;
        res.issue_warning("psi_complex_rhp", &[z]);
        return res;
    }

    // Compute the number of recurrences to apply
    if z.re < 20.0 && z.im.abs() < 20.0 {
        let sp = (20.0 + z.im).sqrt();
        let sn = (20.0 - z.im).sqrt();
        let rhs = sp * sn - z.re;
        if rhs > 0.0 {
            n_recurse = rhs.ceil() as usize;
        }
    }

    // Compute asymptotic at the large value z + n_recurse
    let mut a = psi_complex_asymp(z + n_recurse as f64);
    res.err = 2.0 * std::f64::EPSILON * Complex::new(a.re.abs(), a.im.abs());

    // Descend recursively, if necessary
    for i in (1..(n_recurse + 1)).rev() {
        let zn = z + (i as f64 - 1.0);
        let zn_inv = zn.inv();
        a = z - zn_inv;

        // Accumulate the error, to catch cancellations
        res.err.re += 2.0 * std::f64::EPSILON * zn_inv.re.abs();
        res.err.im += 2.0 * std::f64::EPSILON * zn_inv.im.abs();
    }

    res.val = a;

    res.err.re += 2.0 * std::f64::EPSILON * res.val.re.abs();
    res.err.im += 2.0 * std::f64::EPSILON * res.val.im.abs();

    res
}

/// Compute Psi^{(n)}(x) for x > 0
fn psi_n_xg0(n: usize, x: f64) -> SpecFunResult<f64> {
    if n == 0 {
        psi_x_e(x)
    } else {
        // Abramowitz + Stegun 6.4.10
        let hz = hzeta_e(n as f64 + 1.0, x);
        let lnnf = lnfact_e(n);
        let mut result = exp_mult_err_e(lnnf.val, lnnf.err, hz.val, hz.err);

        if n % 2 == 0 {
            result.val *= -1.0;
        }

        result
    }
}

// --------------------------------
// ------ External Functions ------
// --------------------------------

/// Compute the digamma function
pub(crate) fn digamma_int_e(n: usize) -> SpecFunResult<f64> {
    if n == 0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("digamma_int_e", &[n as f64]);
        result
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
        let mut err = std::f64::EPSILON * (nf.ln().abs() + (0.5 / nf).abs() + ser.abs());
        err += std::f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

/// Compute the digamma function for float
pub(crate) fn digamma_e(x: f64) -> SpecFunResult<f64> {
    psi_x_e(x)
}

/// Compute the trigamma function for integer
pub(crate) fn trigamma_int_e(n: usize) -> SpecFunResult<f64> {
    if n == 0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("trigamma_int_e", &[n as f64]);
        result
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

/// Compute the trigamma function for float
pub(crate) fn trigamma_e(x: f64) -> SpecFunResult<f64> {
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

/// Compute the polygamma function for float
pub(crate) fn polygamma_e(n: usize, x: f64) -> SpecFunResult<f64> {
    if n == 0 {
        digamma_e(x)
    } else if n == 1 {
        trigamma_e(x)
    } else if x <= 0.0 {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        result.issue_warning("polygamma_e", &[n as f64]);
        result
    } else {
        let hz = hzeta_e(n as f64 + 1.0, x);
        let lnnf = lnfact_e(n);
        let mut result = exp_mult_err_e(lnnf.val, lnnf.err, hz.val, hz.err);
        if n % 2 == 0 {
            result.val *= -1.0;
        }
        result
    }
}
