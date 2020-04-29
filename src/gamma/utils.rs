use crate::cheb::cheb_eval_e;
use crate::consts::{LNPI, ROOT4_DBL_EPS};
use crate::exp::{exp_err_e, exp_mult_err_e};
use crate::result::{ErrorCode, SpecFunResult};
use crate::zeta::HurwitzZeta;
use log::warn;
use num::complex::Complex;
use num::Float;

use super::data::*;

// ------------------
// Internal Functions
// ------------------

/// Compute the factorial of a number, returning an f64
pub fn factorial(n: usize) -> f64 {
    if n < fact_table.len() {
        fact_table[n]
    } else {
        warn!(
            "SpecFunOverflowError: Overflow occured in factorial at n = {}",
            n
        );
        std::f64::NAN
    }
}

/// Compute the double factorial of a number, returning an f64
pub fn double_factorial(n: usize) -> f64 {
    if n < double_fact_table.len() {
        double_fact_table[n]
    } else {
        warn!(
            "SpecFunOverflowError: Overflow occured in double_factorial at n = {}",
            n
        );
        std::f64::NAN
    }
}

/// Compute Log(Gamma(z)) using Lanczos method for complex numbers
pub fn lngamma_lanczos_complex_e(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut y = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(f64::EPSILON, f64::EPSILON),
        code: ErrorCode::CONTINUE,
    };

    let zz = Complex::new(z.re - 1.0, z.im);
    let mut ag = Complex::new(lanczos_7_c[0], 0.0);
    for k in 1..9 {
        let r = zz.re + (k as f64);
        let i = zz.im;
        let a = lanczos_7_c[k] / (r * r + i * i);
        ag.re += a * r;
        ag.im -= a * i;
    }

    let log1 = (zz + 7.5).ln();
    let logag = ag.ln();

    y.val.re = (zz.re + 0.5) * log1.re - zz.im * log1.im - (zz.re + 7.5)
        + 0.9189385332046727418
        + logag.re;
    y.val.im = zz.im * log1.re + (zz.re + 0.5) * log1.im - zz.im + logag.im;

    y.err.re = 4.0 * f64::EPSILON * y.val.re.abs();
    y.err.im = 4.0 * f64::EPSILON * y.val.im.abs();

    y
}

/// Compute Log(Gamma(x)) using Lanczos method
pub fn lngamma_lanczos_e(x: f64) -> SpecFunResult<f64> {
    let xx = x - 1.0;
    let ag = lanczos_7_c
        .iter()
        .skip(1)
        .enumerate()
        .fold(lanczos_7_c[0], |acc, (k, val)| {
            acc + val / (xx + (k as f64))
        });
    let term1 = (xx + 0.5) * ((xx + 7.5) * (-1.0).exp()).ln();
    let term2 = 0.9189385332046727418 + ag.ln();

    let val = term1 + (term2 - 7.0);
    let err = f64::EPSILON * (2.0 * (term1.abs() + term2.abs() + 7.0) + val.abs());

    SpecFunResult { val, err }
}

/// Calculate series for g(eps) = Gamma(eps) eps - 1/(1+eps) - eps / 2, as well as its sign
pub fn lngamma_sgn_0_e(eps: f64) -> (SpecFunResult<f64>, f64) {
    let c1 = -0.07721566490153286061;
    let c2 = -0.01094400467202744461;
    let c3 = 0.09252092391911371098;
    let c4 = -0.01827191316559981266;
    let c5 = 0.01800493109685479790;
    let c6 = -0.00685088537872380685;
    let c7 = 0.00399823955756846603;
    let c8 = -0.00189430621687107802;
    let c9 = 0.00097473237804513221;
    let c10 = -0.00048434392722255893;
    let g6 = c6 + eps * (c7 + eps * (c8 + eps * (c9 + eps * c10)));
    let g = eps * (c1 + eps * (c2 + eps * (c3 + eps * (c4 + eps * (c5 + eps * g6)))));

    let gee = g + 1.0 / (1.0 + eps) + 0.5 * eps;
    let val = gee * eps.abs().recip();
    let err = 4.0 * f64::EPSILON * val.abs();
    (SpecFunResult { val, err }, val.signum())
}

pub fn lngamma_sgn_sing_e(n: i32, eps: f64) -> (SpecFunResult<f64>, f64) {
    if eps.abs() < f64::EPSILON {
        (SpecFunResult { val: 0.0, err: 0.0 }, 0.0)
    } else if n == 1 {
        let c0 = 0.07721566490153286061;
        let c1 = 0.08815966957356030521;
        let c2 = -0.00436125434555340577;
        let c3 = 0.01391065882004640689;
        let c4 = -0.00409427227680839100;
        let c5 = 0.00275661310191541584;
        let c6 = -0.00124162645565305019;
        let c7 = 0.00065267976121802783;
        let c8 = -0.00032205261682710437;
        let c9 = 0.00016229131039545456;
        let g5 = c5 + eps * (c6 + eps * (c7 + eps * (c8 + eps * c9)));
        let g = eps * (c0 + eps * (c1 + eps * (c2 + eps * (c3 + eps * (c4 + eps * g5)))));

        // calculate eps gamma(-1+eps), a negative quantity
        let gam_e = g - 1.0 - 0.5 * eps * (1.0 + 3.0 * eps) / (1.0 - eps * eps);
        let val = (gam_e.abs() * eps.abs().recip()).ln();
        let err = 2.0 * f64::EPSILON * val.abs();
        (
            SpecFunResult { val, err },
            if eps > 0.0 { -1.0 } else { 1.0 },
        )
    } else {
        let cs1 = -1.6449340668482264365;
        let cs2 = 0.8117424252833536436;
        let cs3 = -0.1907518241220842137;
        let cs4 = 0.0261478478176548005;
        let cs5 = -0.0023460810354558236;
        let e2 = eps * eps;
        let sin_ser = 1.0 + e2 * (cs1 + e2 * (cs2 + e2 * (cs3 + e2 * (cs4 + e2 * cs5))));

        // Calculate series for ln(gamma(1+n-eps))
        let aeps = eps.abs();

        let mut psi0 = SpecFunResult { val: 0.0, err: 0.0 };
        let mut psi1 = SpecFunResult { val: 0.0, err: 0.0 };
        let mut psi2 = SpecFunResult { val: 0.0, err: 0.0 };
        let mut psi3 = SpecFunResult { val: 0.0, err: 0.0 };
        let mut psi4 = SpecFunResult { val: 0.0, err: 0.0 };
        let mut psi5 = SpecFunResult { val: 0.0, err: 0.0 };
        let mut psi6 = SpecFunResult { val: 0.0, err: 0.0 };

        let c0 = lnfact_int_e(n as usize);
        psi0 = digamma_int_e((n + 1) as usize);
        psi1 = trigamma_int_e((n + 1) as usize);

        (SpecFunResult { val: 0.0, err: 0.0 }, 0.0)
    }
}

pub fn lnfact_int_e(n: usize) -> SpecFunResult<f64> {
    if n < fact_table.len() {
        let val = fact_table[n].ln();
        SpecFunResult {
            val: val,
            err: 2.0 * f64::EPSILON * val.abs(),
        }
    } else {
        lngamma_e((n as f64) + 1.0)
    }
}

/// Compute log(Gamma(1+eps))/eps using the (2,2) pade
/// approximate plus a correction series
pub fn lngamma_1_pade_e(eps: f64) -> SpecFunResult<f64> {
    let n1 = -1.0017419282349508699871138440;
    let n2 = 1.7364839209922879823280541733;
    let d1 = 1.2433006018858751556055436011;
    let d2 = 5.0456274100274010152489597514;
    let num = (eps + n1) * (eps + n2);
    let den = (eps + d1) * (eps + d2);
    let pade = 2.0816265188662692474880210318 * num * den.recip();

    let c0 = 0.004785324257581753;
    let c1 = -0.01192457083645441;
    let c2 = 0.01931961413960498;
    let c3 = -0.02594027398725020;
    let c4 = 0.03141928755021455;
    let eps5 = eps * eps * eps * eps * eps;
    let corr = eps5 * (c0 + eps * (c1 + eps * (c2 + eps * (c3 + c4 * eps))));

    // Error in case wanted for future use
    // let err = 2.0 * T::epsilon() * eps * (pade + corr);
    let val = eps * (pade + corr);
    let err = 2.0 * f64::EPSILON * val.abs();
    SpecFunResult { val, err }
}

/// Compute log(Gamma(2+eps))/eps using the (2,2) pade
/// approximate plus a correction series
pub fn lngamma_2_pade_e(eps: f64) -> SpecFunResult<f64> {
    let n1 = 1.000895834786669227164446568;
    let n2 = 4.209376735287755081642901277;
    let d1 = 2.618851904903217274682578255;
    let d2 = 10.85766559900983515322922936;
    let num = (eps + n1) * (eps + n2);
    let den = (eps + d1) * (eps + d2);
    let pade = 2.85337998765781918463568869 * num / den;
    let c0 = 0.0001139406357036744;
    let c1 = -0.0001365435269792533;
    let c2 = 0.0001067287169183665;
    let c3 = -0.0000693271800931282;
    let c4 = 0.0000407220927867950;
    let eps5 = eps * eps * eps * eps * eps;
    let corr = eps5 * (c0 + eps * (c1 + eps * (c2 + eps * (c3 + c4 * eps))));

    let val = eps * (pade + corr);
    let err = 2.0 * f64::EPSILON * val.abs();

    SpecFunResult { val, err }
}

/// Compute the digamma function
fn digamma_int_e(n: usize) -> SpecFunResult<f64> {
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
        let mut err = std::f64::EPSILON * (nf.ln().abs() + (0.5 / nf).abs() + ser.abs());
        err += std::f64::EPSILON * val.abs();
        SpecFunResult { val, err }
    }
}

/// Compute the trigamma function
fn trigamma_int_e(n: usize) -> SpecFunResult<f64> {
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

pub fn gammastar_ser_e(x: f64) -> SpecFunResult<f64> {
    let y = 1.0 / (x * x);
    let c0 = 1.0 / 12.0;
    let c1 = -1.0 / 360.0;
    let c2 = 1.0 / 1260.0;
    let c3 = -1.0 / 1680.0;
    let c4 = 1.0 / 1188.0;
    let c5 = -691.0 / 360360.0;
    let c6 = 1.0 / 156.0;
    let c7 = -3617.0 / 122400.0;
    let ser = c0 + y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * c7))))));

    let val = (ser / x).exp();
    let err = 2.0 * f64::EPSILON * val * (ser / x).max(1.0);

    SpecFunResult { val, err }
}

/// Compute Gamma(x) for x >= 1/2
pub fn gamma_x_gt_half_e(x: f64) -> SpecFunResult<f64> {
    if (x - 0.5).abs() < f64::EPSILON {
        // Error term
        let val = 1.77245385090551602729817;
        let err = f64::EPSILON * val;
        SpecFunResult { val, err }
    } else if x <= 1.0 + (fact_table.len() as f64) - 1.0 && (x - x.floor()).abs() < f64::EPSILON {
        let n = x.floor() as usize;
        let val = fact_table[n - 1];
        let err = val * f64::EPSILON;
        SpecFunResult { val, err }
    } else if (x - 1.0).abs() < 0.01 {
        let eps = x - 1.0;
        let c1 = 0.4227843350984671394;
        let c2 = -0.01094400467202744461;
        let c3 = 0.09252092391911371098;
        let c4 = -0.018271913165599812664;
        let c5 = 0.018004931096854797895;
        let c6 = -0.006850885378723806846;
        let c7 = 0.003998239557568466030;
        let err = f64::EPSILON;
        let val = 1.0 * x.recip()
            + eps
                * (c1 + eps * (c2 + eps * (c3 + eps * (c4 + eps * (c5 + eps * (c6 + eps * c7))))));
        SpecFunResult { val, err }
    } else if (x - 2.0).abs() < 0.01 {
        let eps = x - 2.0;
        let c1 = 0.4227843350984671394;
        let c2 = 0.4118403304264396948;
        let c3 = 0.08157691924708626638;
        let c4 = 0.07424901075351389832;
        let c5 = -0.00026698206874501476832;
        let c6 = 0.011154045718130991049;
        let c7 = -0.002852645821155340816;
        let c8 = 0.0021039333406973880085;

        let err = f64::EPSILON;
        let val = 1.0
            + eps
                * (c1
                    + eps
                        * (c2
                            + eps
                                * (c3
                                    + eps
                                        * (c4 + eps * (c5 + eps * (c6 + eps * (c7 + eps * c8)))))));
        SpecFunResult { val, err }
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
        let c = cheb_eval_e(t, &gamma_5_10_data, -1.0, 1.0);
        let val = c.val.exp() * gamma_8;
        let err = val * c.err + 2.0 * f64::EPSILON * val;
        SpecFunResult { val, err }
    } else if x < GSL_SF_GAMMA_XMAX {
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
        warn!(
            "SpecFunOverflowError: Overflow occured at gamma_x_gt_half_e({})",
            x
        );
        SpecFunResult { val: 0.0, err: 0.0 }
    }
}

/// Compute the digamma function for either positive or negative `x`.
pub fn psi_x_e(x: f64) -> SpecFunResult<f64> {
    let y = x.abs();

    if y < f64::EPSILON || (x + 1.0).abs() < f64::EPSILON || (x + 2.0).abs() < f64::EPSILON {
        warn!("SpecFunDomainError: Tried to evaluate Digamma at {}", x);
        SpecFunResult {
            val: std::f64::NAN,
            err: std::f64::NAN,
        }
    } else if y >= 2.0 {
        let t = 8.0 / (y * y) - 1.0;
        let result_c = cheb_eval_e(t, &apsics_data, -1.0, 1.0);
        if x < 0.0 {
            let s = (x * std::f64::consts::PI).sin();
            let c = (x * std::f64::consts::PI).cos();
            if s.abs() < 2.0 * f64::MIN_POSITIVE.sqrt() {
                warn!("SpecFunDomainError: Tried to evaluate Digamma at {}", x);
                SpecFunResult {
                    val: std::f64::NAN,
                    err: std::f64::NAN,
                }
            } else {
                let val = y.ln() - 0.5 / x + result_c.val - std::f64::consts::PI * c / s;
                let mut err = std::f64::consts::PI * x.abs() * f64::EPSILON / (s * s);
                err += result_c.err;
                err += f64::EPSILON * val.abs();
                SpecFunResult { val, err }
            }
        } else {
            let val = y.ln() - 0.5 / x + result_c.val;
            let mut err = result_c.err;
            err += f64::EPSILON * val;
            SpecFunResult { val, err }
        }
    } else {
        if x < -1.0 {
            let v = x + 2.0;
            let t1 = 1.0 / x;
            let t2 = 1.0 / (x + 1.0);
            let t3 = 1.0 / v;
            let result_c = cheb_eval_e(2.0 * v - 1.0, &psics_data, -1.0, 1.0);

            let val = -(t1 + t2 + t3) + result_c.val;
            let mut err = f64::EPSILON * (t1.abs() + (x / (t2 * t2)).abs() + (x / (t3 * t3)).abs());
            err += result_c.err;
            err += f64::EPSILON * val.abs();
            SpecFunResult { val, err }
        } else if x < 0.0 {
            let v = x + 1.0;
            let t1 = 1.0 / x;
            let t2 = 1.0 / v;
            let result_c = cheb_eval_e(2.0 * v - 1.0, &psics_data, -1.0, 1.0);

            let val = -(t1 + t2) + result_c.val;
            let mut err = f64::EPSILON * (t1.abs() + (x / (t2 * t2)).abs());
            err += result_c.err;
            err += f64::EPSILON * val.abs();
            SpecFunResult { val, err }
        } else if x < 1.0 {
            let t1 = 1.0 / x;
            let result_c = cheb_eval_e(2.0 * x - 1.0, &psics_data, -1.0, 1.0);

            let val = -t1 + result_c.val;
            let mut err = f64::EPSILON * t1;
            err += result_c.err;
            err += f64::EPSILON * val;
            SpecFunResult { val, err }
        } else {
            let v = x - 1.0;
            cheb_eval_e(2.0 * v - 1.0, &psics_data, -1.0, 1.0)
        }
    }
}

/// Compute psi(z) for large |z| in the right-half plane
/// ref: [Abramowitz + Stegun, 6.3.18]
pub fn psi_complex_asymp(z: Complex<f64>) -> Complex<f64> {
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
    sum = sum + 1.0;
    sum = sum * (c2 / c1);
    sum = sum * w;
    sum = sum + 1.0;
    sum = sum * c1;
    sum = sum * w;
    sum = sum + 1.0;

    // Correction added to log(z)
    let mut cs = sum * w;
    cs = cs * (-1.0 / 12.0);
    cs = cs + zinv * (-0.5);

    cs + z.ln()
}

/// Compute Psi(z) for complex z in the right-half plane
pub fn psi_complex_rhp(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut n_recurse: usize = 0;
    let mut res = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
    };

    if z.re.abs() < std::f64::EPSILON && z.im.abs() < std::f64::EPSILON {
        warn!(
            "SpecFunDomainError: Tried to evaluate Diagamma at z = {}",
            z
        );
        return SpecFunResult {
            val: Complex::new(std::f64::NAN, std::f64::NAN),
            err: Complex::new(std::f64::NAN, std::f64::NAN),
        };
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

    return res;
}

/// Compute Psi^{(n)}(x) for x > 0
pub fn psi_n_xg0(n: usize, x: f64) -> SpecFunResult<f64> {
    if n == 0 {
        psi_x_e(x)
    } else {
        // Abramowitz + Stegun 6.4.10
        let hz = (n as f64 + 1.0).hzeta_e(x);
        let lnnf = lnfact_int_e(n);
        let mut result = exp_mult_err_e(lnnf.val, lnnf.err, hz.val, hz.err);

        if n % 2 == 0 {
            result.val *= -1.0;
        }

        result
    }
}

// ------------------
// External functions
// ------------------

/// ln(gamma(x)) where x is not a negative integer
///
// Uses real Lanczos method.
// Returns the real part of ln(gamma(x)) when x < 0,
// i.e. ln(|gamma(x)|).
pub fn lngamma_e(x: f64) -> SpecFunResult<f64> {
    if (x - 1.0).abs() < 0.01 {
        // Note that we must amplify the errors
        // from the Pade evaluations because of
        // the way we must pass the argument, i.e.
        // writing (1-x) is a loss of precision
        // when x is near 1.
        let mut result = lngamma_1_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 1.0).abs());
        return result;
    } else if (x - 2.0).abs() < 0.01 {
        let mut result = lngamma_2_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 2.0).abs());
        return result;
    } else if x >= 0.5 {
        return lngamma_lanczos_e(x);
    } else if x.abs() < f64::EPSILON {
        warn!("SpecFunDomainError: Tried to evaluate lngamma at x = {}", x);
        return SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
        };
    } else if x.abs() < 0.02 {
        return lngamma_sgn_0_e(x).0;
    } else if x > -0.5 / (f64::EPSILON * std::f64::consts::PI) {
        // Try tp extract a fractional part from x
        let z = 1.0 - x;
        let s = (std::f64::consts::PI * z).sin();
        let abss = s.abs();
        if abss < f64::EPSILON {
            warn!("SpecFunDomainError: Tried to evaluate lngamma at x = {}", x);
            return SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
            };
        } else if abss < std::f64::consts::PI * 0.015 {
            // x is near a negative integer
            if x < std::i32::MIN as f64 + 2.0 {
                warn!(
                    "SpecFunRoundoffError: Tried to evaluate lngamma at x = {}",
                    x
                );
                return SpecFunResult { val: 0.0, err: 0.0 };
            } else {
                let n = -(x - 0.5) as i32;
                let eps = x + n as f64;
                return lngamma_sgn_sing_e(n, eps).0;
            }
        } else {
            let mut result = lngamma_lanczos_e(z);
            result.val = LNPI - (abss.ln() + result.val);
            result.err = 2.0 * f64::EPSILON * result.val.abs() + result.err;
            return result;
        }
    } else {
        // |x| was too large to extract any fraction part
        warn!("SpecFunRoundoffError: |x| was too large to extract any fraction part in lngamma at x = {}",x);
        return SpecFunResult { val: 0.0, err: 0.0 };
    }
}

/// Compute ln(gamma(x)) where x is not a negative integer.
///
// Uses real Lanczos method. Determines
// the sign of gamma(x) as well as ln(|gamma(x)|) for x < 0.
// So gamma(x) = sgn * exp(lngamma_sgn_e(x).res).
pub fn lngamma_sgn_e(x: f64) -> (SpecFunResult<f64>, f64) {
    if (x - 1.0).abs() < 0.01 {
        // Note that we must amplify the errors
        // from the Pade evaluations because of
        // the way we must pass the argument, i.e.
        // writing (1-x) is a loss of precision
        // when x is near 1.
        let mut result = lngamma_1_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 1.0).abs());
        return (result, 1.0);
    } else if (x - 2.0).abs() < 0.01 {
        let mut result = lngamma_2_pade_e(x - 1.0);
        result.err *= 1.0 / (f64::EPSILON + (x - 2.0).abs());
        return (result, 1.0);
    } else if x >= 0.5 {
        return (lngamma_lanczos_e(x), 1.0);
    } else if x.abs() < f64::EPSILON {
        warn!("SpecFunDomainError: Tried to evaluate lngamma at x = {}", x);
        return (
            SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
            },
            0.0,
        );
    } else if x.abs() < 0.02 {
        return lngamma_sgn_0_e(x);
    } else if x > -0.5 / (f64::EPSILON * std::f64::consts::PI) {
        // Try to extract a fractional part from x
        let z = 1.0 - x;
        let s = (std::f64::consts::PI * z).sin();
        let abss = s.abs();
        if abss < f64::EPSILON {
            warn!(
                "SpecFunDomainError: Tried to evaluate lngamma_sgn at x = {}",
                x
            );
            return (
                SpecFunResult {
                    val: f64::NAN,
                    err: f64::NAN,
                },
                0.0,
            );
        } else if abss < std::f64::consts::PI * 0.015 {
            // x is near a negative integer
            if x < std::i32::MIN as f64 + 2.0 {
                warn!(
                    "SpecFunRoundoffError: Tried to evaluate lngamma_sgn at x = {}",
                    x
                );
                return (SpecFunResult { val: 0.0, err: 0.0 }, 0.0);
            } else {
                let n = -(x - 0.5) as i32;
                let eps = x + n as f64;
                return lngamma_sgn_sing_e(n, eps);
            }
        } else {
            let mut result = lngamma_lanczos_e(z);
            let sgn = if s > 0.0 { 1.0 } else { -1.0 };
            result.val = LNPI - (abss.ln() + result.val);
            result.err = 2.0 * f64::EPSILON * result.val.abs() + result.err;
            return (result, sgn);
        }
    } else {
        // |x| was too large to extract any fraction part
        warn!("SpecFunRoundoffError: |x| was too large to extract any fraction part in lngamma_sgn at x = {}",x);
        return (SpecFunResult { val: 0.0, err: 0.0 }, 0.0);
    }
}
/// Gamma(x), where x is not a negative integer
///
/// Uses real Lanczos method.
pub fn gamma_e(x: f64) -> SpecFunResult<f64> {
    if x < 0.5 {
        let rint_x = (x + 0.5).floor() as i32;
        let f_x = x - rint_x as f64;
        let sgn_gamma = if rint_x % 2 == 0 { 1.0 } else { -1.0 };
        let sin_term = sgn_gamma * (std::f64::consts::PI * f_x).sin() / std::f64::consts::PI;

        if sin_term.abs() < f64::EPSILON {
            warn!(
                "SpecFunDomainError: Invalid argument passed to gamma_e: x = {}",
                x
            );
            return SpecFunResult {
                val: f64::NAN,
                err: f64::NAN,
            };
        } else if x > -169.0 {
            let g = gamma_x_gt_half_e(1.0 - x);
            if sin_term.abs() * g.val * f64::MIN_POSITIVE < 1.0 {
                let val = 1.0 / (sin_term * g.val);
                let mut err = (g.err / g.val).abs() * val.abs();
                err += 2.0 * f64::EPSILON * val.abs();
                return SpecFunResult { val, err };
            } else {
                warn!(
                    "SpecFunUnderflowError: Underflow occured in gamma_e at x = {}",
                    x
                );
                return SpecFunResult { val: 0.0, err: 0.0 };
            }
        } else {
            // It is hard to control it here.
            // We can only exponentiate the
            // logarithm and eat the loss of
            // precision.
            let (lng, sgn) = lngamma_sgn_e(x);
            return exp_mult_err_e(lng.val, lng.err, sgn, 0.0);
        }
    } else {
        return gamma_x_gt_half_e(x);
    }
}

/// Regulated Gamma Function, x > 0
///
/// Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
///            = (1 + 1/(12x) + ...), x->Inf
pub fn gammastar_e(x: f64) -> SpecFunResult<f64> {
    if x <= 0.0 {
        warn!(
            "SpecFunDomainError: Domain error occured in gammastar_e as x = {}",
            x
        );
        SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
        }
    } else if x < 0.5 {
        let lg = lngamma_e(x);
        let lx = x.ln();
        let c = 0.5 * (std::f64::consts::LN_2 + LNPI);
        let lnr_val = lg.val - (x - 0.5) * lx + x - c;
        let lnr_err = lg.err + 2.0 * f64::EPSILON * ((x + 0.5) * lx.abs() + c);
        exp_err_e(lnr_val, lnr_err)
    } else if x < 2.0 {
        let t = 4.0 / 3.0 * (x - 0.5) - 1.0;
        cheb_eval_e(t, &gstar_a_data, -1.0, 1.0)
    } else if x < 10.0 {
        let t = 0.25 * (x - 2.0) - 1.0;
        let c = cheb_eval_e(t, &gstar_b_data, -1.0, 1.0);
        let val = c.val / (x * x) + 1.0 + 1.0 / (12.0 * x);
        let mut err = c.err / (x * x);
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult { val, err }
    } else if x < 1.0 / ROOT4_DBL_EPS {
        gammastar_ser_e(x)
    } else if x < 1.0 / f64::EPSILON {
        let xi = 1.0 / x;
        let val = 1.0
            + xi / 12.0 * (1.0 + xi / 24.0 * (1.0 - xi * (139.0 / 180.0 + 571.0 / 8640.0 * xi)));
        let err = 2.0 * f64::EPSILON * val.abs();
        SpecFunResult { val, err }
    } else {
        SpecFunResult {
            val: 1.0,
            err: x.recip(),
        }
    }
}

/// Compute 1 / gamma(x)
///
/// Uses real Lanczos method.
pub fn gammainv_e(x: f64) -> SpecFunResult<f64> {
    if x <= 0.0 && (x - x.floor()).abs() < f64::EPSILON {
        SpecFunResult { val: 0.0, err: 0.0 }
    } else if x < 0.5 {
        let (lng, sgn) = lngamma_sgn_e(x);
        if lng.val == f64::NAN {
            SpecFunResult { val: 0.0, err: 0.0 }
        } else {
            exp_mult_err_e(-lng.val, lng.err, sgn, 0.0)
        }
    } else {
        let g = gamma_x_gt_half_e(x);
        let val = 1.0 / g.val;
        let mut err = (g.err / g.val).abs() * val.abs();
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult { val, err }
    }
}
