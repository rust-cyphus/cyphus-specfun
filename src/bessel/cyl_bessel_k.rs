use super::bessel_data::*;
use crate::gamma::Gamma;
use crate::result::{result_smash_e, SpecFunCode, SpecFunResult, SpecFunResultE10};

pub(crate) fn cyl_bessel_k0_scaled_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result
    } else if x < 1.0 {
        let lx = x.ln();
        let ex = x.exp();
        let x2 = x * x;
        result.val = ex * ((*K0_POLY).eval(x2) - lx * (1.0 + 0.25 * x2 * (*K0_POLY).eval(x2)));
        result.err = ex * (1.6 + lx.abs() * 0.6) * f64::EPSILON;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    } else if x <= 8.0 {
        let sx = x.sqrt();
        let c = (*AK0_CHEB).eval((16.0 / x - 9.0) / 7.0);
        result.val = (1.203125 + c.val) / sx;
        result.err = c.err / sx;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        let sx = x.sqrt();
        let c = (*AK02_CHEB).eval(16.0 / x - 1.0);
        result.val = (1.25 + c.val) / sx;
        result.err = (c.err + f64::EPSILON) / sx;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

pub(crate) fn cyl_bessel_k0_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result
    } else if x < 1.0 {
        let lx = x.ln();
        let x2 = x * x;
        result.val = (*K0_POLY).eval(x2) - lx * (1.0 + 0.25 * x2 * (*I0_POLY).eval(x2));
        result.err = (1.6 + lx.abs() * 0.6) * f64::EPSILON;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        let k0_scaled = cyl_bessel_k0_scaled_e(x);
        crate::exp::core::exp_mult_err_e(-x, f64::EPSILON * x.abs(), k0_scaled.val, k0_scaled.err)
    }
}

pub(crate) fn cyl_bessel_k1_scaled_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
    } else if x < 2.0 * f64::MIN_POSITIVE {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
    } else if x < 1.0 {
        let lx = x.ln();
        let ex = x.exp();
        let x2 = x * x;
        let t = 0.25 * x2;
        let i1 = 0.5 * x * (1.0 + t * (0.5 + t * (*I1_POLY).eval(t)));
        result.val = ex * (x2 * (*K1_POLY).eval(x2) + x * lx * i1 + 1.0) / x;
        result.err = ex * (1.6 + lx.abs() * 0.6) * f64::EPSILON;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else if x <= 8.0 {
        let sx = x.sqrt();
        let c = (*AK1_CHEB).eval((16.0 / x - 9.0) / 7.0);
        result.val = (1.375 + c.val) / sx; /* 1.375 = 11/8 */
        result.err = c.err / sx;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else {
        let sx = x.sqrt();
        let c = (*AK12_CHEB).eval(16.0 / x - 1.0);
        result.val = (1.25 + c.val) / sx;
        result.err = c.err / sx;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    }

    result
}

pub(crate) fn cyl_bessel_k1_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
    } else if x < 2.0 * f64::MIN_POSITIVE {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
    } else if x < 1.0 {
        let lx = x.ln();
        let x2 = x * x;
        let t = 0.25 * x2;
        let i1 = 0.5 * x * (1.0 + t * (0.5 + t * (*I1_POLY).eval(t)));
        result.val = (x2 * (*K1_POLY).eval(x2) + x * lx * i1 + 1.0) / x;
        result.err = (1.6 + lx.abs() * 0.6) * f64::EPSILON;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else {
        let k1_scaled = cyl_bessel_k1_scaled_e(x);
        result = crate::exp::core::exp_mult_err_e(-x, 0.0, k1_scaled.val, k1_scaled.err);
        result.err = result.val.abs() * (f64::EPSILON * x.abs() + k1_scaled.err / k1_scaled.val);
    }
    result
}

// [Abramowitz+Stegun, 9.6.11]
// assumes n >= 1
#[allow(dead_code)]
fn cyl_bessel_kn_scaled_small_x(n: i32, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let y = 0.25 * x * x;
    let ln_x_2 = (0.5 * x).ln();
    let ex = x.exp();
    let ln_nm1_fact = crate::gamma::mono::lnfact_e((n - 1) as usize);
    let ln_pre1 = -n as f64 * ln_x_2 + ln_nm1_fact.val;

    if ln_pre1 > crate::consts::LN_DBL_MAX - 3.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
    }

    let mut sum1 = 1.0;
    let mut k_term = 1.0;
    for k in 1..n {
        k_term *= -y / (k * (n - k)) as f64;
        sum1 += k_term;
    }
    let term1 = 0.5 * (ln_pre1).exp() * sum1;

    let pre2 = 0.5 * (n as f64 * ln_x_2).exp();
    let term2 = if pre2 > 0.0 {
        let kmax = 20;
        let mut yk = 1.0;
        let mut k_fact = 1.0;
        let mut psi_kp1 = -crate::consts::EUL_GAMMA;
        let psi_n = n.digamma_e();
        let mut npk_fact = crate::gamma::mono::fact_e(n as usize);
        let mut psi_npkp1 = psi_n.val + 1.0 / n as f64;
        let mut sum2 = (psi_kp1 + psi_npkp1 - 2.0 * ln_x_2) / npk_fact.val;
        for k in 1..kmax {
            psi_kp1 += 1.0 / k as f64;
            psi_npkp1 += 1.0 / (n + k) as f64;
            k_fact *= k as f64;
            npk_fact.val *= (n + k) as f64;
            yk *= y;
            k_term = yk * (psi_kp1 + psi_npkp1 - 2.0 * ln_x_2) / (k_fact * npk_fact.val);
            sum2 += k_term;
        }
        pre2 * sum2 * (if n % 2 == 1 { -1.0 } else { 1.0 })
    } else {
        0.0
    };

    result.val = ex * (term1 + term2);
    result.err = ex * f64::EPSILON * ((ln_pre1).abs() * (term1).abs() + (term2).abs());
    result.err += 2.0 * f64::EPSILON * result.val.abs();
    result
}
#[allow(dead_code)]
pub(crate) fn cyl_bessel_kn_scaled_e(n: i32, x: f64) -> SpecFunResult<f64> {
    let n = n.abs(); /* K(-n, z) = K(n, z) */
    let mut result = SpecFunResult::<f64>::default();

    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result
    } else if n == 0 {
        return cyl_bessel_k0_scaled_e(x);
    } else if n == 1 {
        return cyl_bessel_k1_scaled_e(x);
    } else if x <= 5.0 {
        return cyl_bessel_kn_scaled_small_x(n, x);
    } else if crate::consts::ROOT3_DBL_EPS * x > 0.25 * (n * n + 1) as f64 {
        return super::bessel_helpers::besselkv_scaled_asympx_e(n as f64, x);
    } else if (0.29 / (n * n) as f64).min(0.5 / ((n * n) as f64 + x * x))
        < crate::consts::ROOT3_DBL_EPS
    {
        return super::bessel_helpers::besselkv_scaled_asymp_unif_e(n as f64, x);
    } else {
        /* Upward recurrence. [Gradshteyn + Ryzhik, 8.471.1] */
        let two_over_x = 2.0 / x;
        let r_b_jm1 = cyl_bessel_k0_scaled_e(x);
        let r_b_j = cyl_bessel_k1_scaled_e(x);
        let mut b_jm1 = r_b_jm1.val;
        let mut b_j = r_b_j.val;

        for j in 1..n {
            let b_jp1 = b_jm1 + j as f64 * two_over_x * b_j;
            b_jm1 = b_j;
            b_j = b_jp1;
        }

        result.val = b_j;
        result.err = n as f64
            * (b_j.abs() * ((r_b_jm1.err / r_b_jm1.val).abs() + (r_b_j.err / r_b_j.val).abs()));
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

#[allow(dead_code)]
pub(crate) fn cyl_bessel_kn_e(n: i32, x: f64) {
    let mut result = cyl_bessel_kn_scaled_e(n, x);
    let ex = (-x).exp();
    result.val *= ex;
    result.err *= ex;
    result.err += x * f64::EPSILON * result.val.abs();
}

fn cyl_bessel_kv_scaled_e10_e(nu: f64, x: f64) -> SpecFunResultE10<f64> {
    let mut result = SpecFunResultE10::<f64>::default();
    if x <= 0.0 || nu < 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
    } else {
        let N = (nu + 0.5) as i32;
        let mu = nu - N as f64; // -1/2 <= mu <= 1/2
        let e10 = 0;

        let (mut K_mu, mut K_mup1, mut Kp_mu) = if x < 2.0 {
            super::bessel_helpers::besselk_scaled_temme(mu, x)
        } else {
            super::bessel_helpers::besselk_scaled_steed_temme_cf2(mu, x)
        };

        /* recurse forward to obtain K_num1, K_nu */
        let mut K_nu = K_mu.val;
        let mut K_nup1 = K_mup1.val;

        for n in 0..N {
            let K_num1 = K_nu;
            K_nu = K_nup1;
            /* rescale the recurrence to avoid overflow */
            if K_nu.abs() > crate::consts::SQRT_DBL_MAX {
                let p = (((K_nu).abs()).ln() / std::f64::consts::LN_10).floor() as i32;
                let factor = 10f64.powi(p);
                K_num1 /= factor;
                K_nu /= factor;
                e10 += p as i32;
            }
            K_nup1 = 2.0 * (mu + n as f64 + 1.0) / x * K_nu + K_num1;
        }

        result.val = K_nu;
        result.err = 2.0 * f64::EPSILON * (N + 4) as f64 * result.val.abs();
        result.e10 = e10;
    }
    result
}

#[allow(dead_code)]
pub(crate) fn cyl_bessel_kv_scaled_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x <= 0.0 || nu < 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
    } else {
        let mut result_e10 = cyl_bessel_kv_scaled_e10_e(nu, x);
        result_smash_e(&mut result_e10, &mut result);
    }
    result
}

#[allow(dead_code)]
pub(crate) fn cyl_bessel_kv_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let b = cyl_bessel_kv_scaled_e(nu, x);
    crate::exp::core::exp_mult_err_e(-x, 0.0, b.val, b.err)
}
