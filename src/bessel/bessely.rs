use super::data::*;
use crate::gamma::Gamma;
use crate::result::{SpecFunCode, SpecFunResult};

#[allow(dead_code)]
pub(super) fn bessely0_e(x: f64) -> SpecFunResult<f64> {
    let two_over_pi = 2.0 / std::f64::consts::PI;
    let xmax = 1.0 / f64::EPSILON;

    let mut result = SpecFunResult::<f64>::default();

    if x <= 0.0 {
        result.code = SpecFunCode::DomainErr;
    } else if x < 4.0 {
        let j0 = super::besselj::besselj0_e(x);
        let c = (*BY0_CHEB).eval(0.125 * x * x - 1.0);
        result.val = two_over_pi * (-std::f64::consts::LN_2 + x.ln()) * j0.val + 0.375 + c.val;
        result.err = 2.0 * f64::EPSILON * result.val.abs() + c.err;
        return result;
    } else if x < xmax {
        // Leading behaviour of phase is x, which is exact,
        // so the error is bounded.
        let z = 32.0 / (x * x) - 1.0;
        let c1 = (*BESSEL_AMP_PHASE_BM0_CHEB).eval(z);
        let c2 = (*BESSEL_AMP_PHASE_BTH0_CHEB).eval(z);
        let sp = super::helpers::bessel_sin_pi4_e(x, c2.val / x);
        let sqrtx = x.sqrt();
        let ampl = (0.75 + c1.val) / sqrtx;
        result.val = ampl * sp.val;
        result.err = sp.val.abs() * c1.err / sqrtx + ampl.abs() * sp.err;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.code = SpecFunCode::UnderflowErr;
    }
    result
}

#[allow(dead_code)]
pub(super) fn bessely1_e(x: f64) -> SpecFunResult<f64> {
    let two_over_pi = 2.0 / std::f64::consts::PI;
    let xmin = 1.571 * f64::MIN_POSITIVE;
    let x_small = 2.0 * crate::consts::SQRT_DLB_EPS;
    let xmax = 1.0 / f64::EPSILON;

    let mut result = SpecFunResult::<f64>::default();

    if x <= 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < xmin {
        result.code = SpecFunCode::OverflowErr;
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
    } else if x < x_small {
        let lnterm = (0.5 * x).ln();
        let j1 = super::besselj::besselj1_e(x);
        let c = (*BY1_CHEB).eval(-1.0);
        result.val = two_over_pi * lnterm * j1.val + (0.5 + c.val) / x;
        result.err = lnterm.abs() * ((f64::EPSILON * j1.val).abs() + j1.err) + c.err / x;
    } else if x < 4.0 {
        let lnterm = (0.5 * x).ln();
        let c = (*BY1_CHEB).eval(0.125 * x * x - 1.0);
        let j1 = super::besselj::besselj1_e(x);
        result.val = two_over_pi * lnterm * j1.val + (0.5 + c.val) / x;
        result.err = lnterm.abs() * ((f64::EPSILON * j1.val).abs() + j1.err) + c.err / x;
    } else if x < xmax {
        let z = 32.0 / (x * x) - 1.0;
        let ca = (*BESSEL_AMP_PHASE_BM1_CHEB).eval(z);
        let ct = (*BESSEL_AMP_PHASE_BTH1_CHEB).eval(z);
        let cp = super::helpers::bessel_cos_pi4_e(x, ct.val / x);
        let sqrtx = x.sqrt();
        let ampl = (0.75 + ca.val) / sqrtx;
        result.val = -ampl * cp.val;
        result.err = cp.val.abs() * ca.err / sqrtx + ampl.abs() * cp.err;
        result.err += f64::EPSILON * result.val.abs();
    } else {
        result.code = SpecFunCode::UnderflowErr;
        result.err = f64::MIN_POSITIVE;
    }
    result
}

/* assumes n >= 1 */
#[allow(dead_code)]
fn besselyn_small_x(n: usize, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let y = 0.25 * x * x;
    let ln_x_2 = (0.5 * x).ln();

    let ln_nm1_fact = crate::gamma::mono::lnfact_e(n - 1);

    let ln_pre1 = -(n as f64) * ln_x_2 + ln_nm1_fact.val;
    if ln_pre1 > crate::consts::LN_DBL_MAX - 3.0 {
        result.code = SpecFunCode::OverflowErr;
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
    }

    let mut sum1 = 1.0;
    let mut k_term = 1.0;
    for k in 1..n {
        k_term *= y / (k * (n - k)) as f64;
        sum1 += k_term;
    }
    let term1 = -(ln_pre1.exp()) * sum1 / std::f64::consts::PI;

    let pre2 = -((n as f64 * ln_x_2).exp()) / std::f64::consts::PI;
    let term2 = if pre2.abs() > 0.0 {
        let kmax = 20;
        let mut yk = 1.0;
        let mut k_fact = 1.0;
        let mut psi_kp1 = -crate::consts::EUL_GAMMA;
        let psi_n = n.digamma_e();
        let mut npk_fact = crate::gamma::mono::fact_e(n);
        let mut psi_npkp1 = psi_n.val + 1.0 / n as f64;
        let mut sum2 = (psi_kp1 + psi_npkp1 - 2.0 * ln_x_2) / npk_fact.val;

        for k in 1..kmax {
            psi_kp1 += 1.0 / k as f64;
            psi_npkp1 += 1.0 / (n + k) as f64;
            k_fact *= k as f64;
            npk_fact.val *= (n + k) as f64;
            yk *= -y;
            k_term = yk * (psi_kp1 + psi_npkp1 - 2.0 * ln_x_2) / (k_fact * npk_fact.val);
            sum2 += k_term;
        }
        pre2 * sum2
    } else {
        0.0
    };

    result.val = term1 + term2;
    result.err = f64::EPSILON * (ln_pre1.abs() * term1.abs() + term2.abs());
    result.err += 2.0 * f64::EPSILON * result.val.abs();
    result
}

#[allow(dead_code)]
pub(super) fn besselyn_e(n: i32, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let mut sign = 1.0;

    let n = if n < 0 {
        /* reduce to case n >= 0 */
        if n % 2 == 1 {
            sign = -1.0;
        }
        (-n) as usize
    } else {
        n as usize
    };

    if n == 0 {
        result = bessely0_e(x);
        result.val *= sign;
    } else if n == 1 {
        result = bessely1_e(x);
        result.val *= sign;
    } else {
        if x <= 0.0 {
            result.code = SpecFunCode::DomainErr;
            result.val = f64::NAN;
            result.err = f64::NAN;
        }
        if x < 5.0 {
            result = besselyn_small_x(n, x);
            result.val *= sign;
        } else if crate::consts::ROOT3_DBL_EPS * x > ((n * n) as f64 + 1.0) {
            result = super::helpers::besselyv_asympx_e(n as f64, x);
            result.val *= sign;
        } else if n > 50 {
            result = super::olver::besselyv_asymp_olver_e(n as f64, x);
            result.val *= sign;
        } else {
            let two_over_x = 2.0 / x;
            let r_by = bessely1_e(x);
            let r_bym = bessely0_e(x);
            let mut bym = r_bym.val;
            let mut by = r_by.val;

            for j in 1..n {
                let byp = j as f64 * two_over_x * by - bym;
                bym = by;
                by = byp;
            }
            result.val = sign * by;
            result.err =
                result.val.abs() * ((r_by.err / r_by.val).abs() + (r_bym.err / r_bym.val).abs());
            result.err += 2.0 * f64::EPSILON * result.val.abs();
        }
    }
    result
}

#[allow(dead_code)]
pub(super) fn besselyv_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x <= 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
        result
    } else if nu < 0.0 {
        let jres = super::helpers::besseljv_pos_e(-nu, x);
        let yres = super::helpers::besselyv_pos_e(-nu, x);
        let spi = crate::trig::sincos::sin_pi_e(nu);
        let cpi = crate::trig::sincos::cos_pi_e(nu);

        result.val = cpi.val * yres.val - spi.val * jres.val;
        result.err = (cpi.val * yres.err).abs()
            + (spi.val * jres.err).abs()
            + (cpi.err * yres.val).abs()
            + (spi.err * jres.val).abs();
        result
    } else {
        super::helpers::besselyv_pos_e(nu, x)
    }
}
