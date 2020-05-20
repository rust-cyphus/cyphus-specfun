use crate::consts::{ROOT5_DBL_EPS, SQRT_DLB_EPS};
use crate::gamma::Gamma;
use crate::result::{SpecFunCode, SpecFunResult};

/// These are of use in calculating the oscillating Bessel functions.
// cos(y - pi/4 + eps)
pub(super) fn bessel_cos_pi4_e(y: f64, eps: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let sy = y.sin();
    let cy = y.cos();
    let s = sy + cy;
    let d = sy - cy;
    let abs_sum = cy.abs() + sy.abs();

    let (seps, ceps) = if eps.abs() < ROOT5_DBL_EPS {
        let e2 = eps * eps;
        (
            eps * (1.0 - e2 / 6.0 * (1.0 - e2 / 20.0)),
            1.0 - e2 / 2.0 * (1.0 - e2 / 12.0),
        )
    } else {
        (eps.sin(), eps.cos())
    };

    result.val = (ceps * s - seps * d) / std::f64::consts::SQRT_2;
    result.err =
        2.0 * f64::EPSILON * (ceps.abs() + seps.abs()) * abs_sum / std::f64::consts::SQRT_2;

    // Try to account for error in evaluation of sin(y), cos(y).
    // This is a little sticky because we don't really know
    // how the library routines are doing their argument reduction.
    // However, we will make a reasonable guess.
    // FIXME ?
    result.err *= if y > 1.0 / f64::EPSILON {
        0.5 * y
    } else if y > 1.0 / SQRT_DLB_EPS {
        256.0 * y * SQRT_DLB_EPS
    } else {
        1.0
    };

    result
}

/// These are of use in calculating the oscillating Bessel functions.
// sin(y - pi/4 + eps)
pub(super) fn bessel_sin_pi4_e(y: f64, eps: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    let sy = y.sin();
    let cy = y.cos();
    let s = sy + cy;
    let d = sy - cy;
    let abs_sum = cy.abs() + sy.abs();

    let (seps, ceps) = if eps.abs() < ROOT5_DBL_EPS {
        let e2 = eps * eps;
        (
            eps * (1.0 - e2 / 6.0 * (1.0 - e2 / 20.0)),
            1.0 - e2 / 2.0 * (1.0 - e2 / 12.0),
        )
    } else {
        (eps.sin(), eps.cos())
    };

    result.val = (ceps * d + seps * s) / std::f64::consts::SQRT_2;
    result.err =
        2.0 * f64::EPSILON * (ceps.abs() + seps.abs()) * abs_sum / std::f64::consts::SQRT_2;

    // Try to account for error in evaluation of sin(y), cos(y).
    // This is a little sticky because we don't really know
    // how the library routines are doing their argument reduction.
    // However, we will make a reasonable guess.
    // FIXME ?
    result.err *= if y > 1.0 / f64::EPSILON {
        0.5 * y
    } else if y > 1.0 / SQRT_DLB_EPS {
        256.0 * y * SQRT_DLB_EPS
    } else {
        1.0
    };

    result
}

pub(super) fn bessel_ij_taylor_e(
    nu: f64,
    x: f64,
    sign: i32,
    kmax: i32,
    threshold: f64,
) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::default();

    if nu < 0.0 || x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.issue_warning(
            "bessel_ij_taylor_e",
            &[nu, x, sign as f64, kmax as f64, threshold],
        );
        result
    } else if x == 0.0 {
        if nu == 0.0 {
            result.val = 1.0;
        } else {
            result.val = 0.0;
        }
        result
    } else {
        let mut prefactor = SpecFunResult::default();
        let mut sum = SpecFunResult::<f64>::default();
        if nu == 0.0 {
            prefactor.val = 1.0;
            prefactor.err = 0.0;
            prefactor.code = SpecFunCode::Success;
        } else if nu < (i32::MAX - 1) as f64 {
            // Separate the integer part and use
            // y^nu/gamma(nu+1) = y^N / N! y^f/(N+1)_f to control the error.
            let n = (nu + 0.5).floor() as i32;
            let f = nu - n as f64;
            let poch_factor = (n as f64 + 1.0).poch_e(f);
            let tc_factor = crate::gamma::mono::taylorcoeff_e(n, 0.5 * x);
            let p = (0.5 * x).powf(f);
            prefactor.val = tc_factor.val * p / poch_factor.val;
            prefactor.err = tc_factor.err * p / poch_factor.val;
            prefactor.err += prefactor.val.abs() / poch_factor.val * poch_factor.err;
            prefactor.err += 2.0 * f64::EPSILON * prefactor.val.abs();
        } else {
            let lg = crate::gamma::mono::lngamma_e(nu + 1.0);
            let term1 = nu * (0.5 * x).ln();
            let term2 = lg.val;
            let ln_pre = term1 - term2;
            let ln_pre_err = f64::EPSILON * (term1.abs() + term2.abs()) + lg.err;
            prefactor = crate::exp::core::exp_err_e(ln_pre, ln_pre_err);
        }

        // Evaluate the sum
        {
            let y = sign as f64 * 0.25 * x * x;
            let mut sumk: f64 = 1.0;
            let mut term: f64 = 1.0;
            let mut k: usize = 0;
            for _k in 1..(kmax as usize + 1) {
                k = _k;
                term *= y / ((nu + k as f64) * k as f64);
                sumk += term;
                if (term / sumk).abs() < threshold {
                    break;
                }
            }

            sum.val = sumk;
            sum.err = threshold * sumk.abs();

            sum.code = if k >= kmax as usize {
                SpecFunCode::MaxIterErr
            } else {
                SpecFunCode::Success
            };
        }
        crate::elementary::multiply_err_e(prefactor.val, prefactor.err, sum.val, sum.err)
    }
}

/// Compute the asymptotic approximation to J_{nu}(x) using Hankel's
/// asymptotic expansion.
pub(super) fn bessel_jnu_asympx_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let mu = 4.0 * nu * nu;
    let chi = x - (0.5 * nu + 0.25) * std::f64::consts::PI;
    let mut p: f64 = 0.0;
    let mut q: f64 = 0.0;
    let mut k: usize = 0;
    let mut t: f64 = 1.0;

    while k < 1000 {
        let kd = k as f64;

        let term = (mu - (2.0 * kd - 1.0) * (2.0 * kd - 1.0)) / (kd * (8.0 * x));

        t *= if k == 0 { 1.0 } else { -term };
        let convp: bool = t.abs() < f64::EPSILON * p.abs();
        p += t;

        k += 1;

        t *= term;
        let convq: bool = t.abs() < f64::EPSILON * q.abs();
        q += t;

        // To preserve the consistency of the series we need to exit when
        // p and q have the same number of terms.
        if convp && convq && kd > nu / 2.0 {
            break;
        }

        k += 1
    }

    {
        let pre = (2.0 / (x * std::f64::consts::PI)).sqrt();
        let c = chi.cos();
        let s = chi.sin();
        result.val = pre * (c * p - s * q);
        result.err =
            pre * f64::EPSILON * ((c * p).abs() + (s * q).abs() + t.abs()) * (1.0 + x.abs());
    }

    result
}
