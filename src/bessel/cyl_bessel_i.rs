use super::bessel_data::*;
use crate::result::{SpecFunCode, SpecFunResult};

/// Compute the modified bessel function of the first kind of order 0 scaled
/// by e^-x along with an error estimate.
pub(crate) fn cyl_bessel_i0_scaled_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let y = x.abs();

    if y < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - y;
        result.err = 0.5 * y * y;
    } else if y <= 3.0 {
        let ey = (-y).exp();
        let c = (*BI0_CHEB).eval(y * y / 4.5 - 1.0);
        result.val = ey * (2.75 + c.val);
        result.err = f64::EPSILON * result.val.abs() + ey * c.err;
    } else if y <= 8.0 {
        let sy = y.sqrt();
        let c = (*AI0_CHEB).eval((48.0 / y - 11.0) / 5.0);
        result.val = (0.375 + c.val) / sy;
        result.err = 2.0 * f64::EPSILON * (0.375 + c.val.abs()) / sy;
        result.err += c.err / sy;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else {
        let sy = y.sqrt();
        let c = (*AI02_CHEB).eval(16.0 / y - 1.0);
        result.val = (0.375 + c.val) / sy;
        result.err = 2.0 * f64::EPSILON * (0.375 + c.val.abs()) / sy;
        result.err += c.err / sy;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    }

    result
}

/// Compute the modified bessel function of the first kind of order 0 along
/// with an error estimate.
pub(crate) fn cyl_bessel_i0_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let y = x.abs();

    if y < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0;
        result.err = 0.5 * y * y;
    } else if y <= 3.0 {
        let c = (*BI0_CHEB).eval(y * y / 4.5 - 1.0);
        result.val = 2.75 + c.val;
        result.err = f64::EPSILON * (2.75 + c.val.abs());
        result.err += c.err;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else if y < crate::consts::LN_DBL_MAX - 1.0 {
        let ey = y.exp();
        let b_scaled = cyl_bessel_i0_scaled_e(x);
        result.val = ey * b_scaled.val;
        result.err = ey * b_scaled.err + y * f64::EPSILON * result.val.abs();
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
    }
    result
}

/// Compute the modified bessel function of the first kind of order 1 scaled
/// by e^-x along with an error estimate.
pub(crate) fn cyl_bessel_i1_scaled_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let xmin = 2.0 * f64::MIN_POSITIVE;
    let x_small = 2.828427124746190 * crate::consts::SQRT_DLB_EPS;
    let y = x.abs();

    if y == 0.0 {
        result.val = 0.0;
        result.err = 0.0;
    } else if y < xmin {
        result.val = 0.0;
        result.err = f64::MIN_POSITIVE;
        result.code = SpecFunCode::UnderflowErr;
    } else if y < x_small {
        result.val = 0.5 * x;
        result.err = 0.0;
    } else if y <= 3.0 {
        let ey = (-y).exp();
        let c = (*BI1_CHEB).eval(y * y / 4.5 - 1.0);
        result.val = x * ey * (0.875 + c.val);
        result.err = ey * c.err + y * f64::EPSILON * (result.val).abs();
        result.err += 2.0 * f64::EPSILON * (result.val).abs();
    } else if y <= 8.0 {
        let sy = (y).sqrt();
        let c = (*AI1_CHEB).eval((48.0 / y - 11.0) / 5.0);
        let b = (0.375 + c.val) / sy;
        let s = if x > 0.0 { 1.0 } else { -1.0 };
        result.val = s * b;
        result.err = c.err / sy;
        result.err += 2.0 * f64::EPSILON * (result.val).abs();
    } else {
        let sy = (y).sqrt();
        let c = (*AI12_CHEB).eval(16.0 / y - 1.0);
        let b = (0.375 + c.val) / sy;
        let s = if x > 0.0 { 1.0 } else { -1.0 };
        result.val = s * b;
        result.err = c.err / sy;
        result.err += 2.0 * f64::EPSILON * (result.val).abs();
    }

    result
}

/// Compute the modified bessel function of the first kind of order 1 along
/// with an error estimate.
pub(crate) fn cyl_bessel_i1_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let xmin = 2.0 * f64::MIN_POSITIVE;
    let x_small = 2.828427124746190 * crate::consts::SQRT_DLB_EPS;
    let y = (x).abs();

    if y == 0.0 {
        result.val = 0.0;
        result.err = 0.0;
    } else if y < xmin {
        result.val = 0.0;
        result.err = f64::MIN_POSITIVE;
        result.code = SpecFunCode::UnderflowErr;
    } else if y < x_small {
        result.val = 0.5 * x;
        result.err = 0.0;
    } else if y <= 3.0 {
        let c = (*BI1_CHEB).eval(y * y / 4.5 - 1.0);
        result.val = x * (0.875 + c.val);
        result.err = y * c.err;
        result.err += 2.0 * f64::EPSILON * (result.val).abs();
    } else if y < crate::consts::LN_DBL_MAX {
        let ey = (y).exp();
        let i1_scaled = cyl_bessel_i1_scaled_e(x);
        result.val = ey * i1_scaled.val;
        result.err = ey * i1_scaled.err + y * f64::EPSILON * result.val.abs();
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
    }

    result
}

/// Compute the modified bessel function of the first kind of order n scaled
/// by e^-x along with an error estimate.
pub(crate) fn cyl_bessel_in_scaled_e(n: i32, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let ax = (x).abs();

    let n = n.abs(); /* I(-n, z) = I(n, z) */

    if n == 0 {
        return cyl_bessel_i0_scaled_e(x);
    } else if n == 1 {
        return cyl_bessel_i1_scaled_e(x);
    } else if x == 0.0 {
        result.val = 0.0;
        result.err = 0.0;
    } else if x * x < 10.0 * (n as f64 + 1.0) / std::f64::consts::E {
        let ex = (-ax).exp();
        let t = super::bessel_helpers::bessel_ij_taylor_e(n as f64, ax, 1, 50, f64::EPSILON);
        result.val = t.val * ex;
        result.err = t.err * ex;
        result.err += 2.0 * f64::EPSILON * (result.val).abs();
        if x < 0.0 && n % 2 == 1 {
            result.val = -result.val
        };
    } else if n < 150 && ax < 1e7 {
        let i0_scaled = cyl_bessel_i0_scaled_e(ax);
        let rat = super::bessel_helpers::besseli_cf1_ser(n as f64, ax);
        let mut ikp1 = rat.val * crate::consts::SQRT_DBL_MIN;
        let mut ik = crate::consts::SQRT_DBL_MIN;
        for k in (1..(n + 1)).rev() {
            let ikm1 = ikp1 + 2.0 * k as f64 / ax * ik;
            ikp1 = ik;
            ik = ikm1;
        }
        result.val = i0_scaled.val * (crate::consts::SQRT_DBL_MIN / ik);
        result.err = i0_scaled.err * (crate::consts::SQRT_DBL_MIN / ik);
        result.err += 2.0 * f64::EPSILON * (result.val).abs();
        if x < 0.0 && n % 2 == 1 {
            result.val = -result.val;
        }
    } else if (0.29 / (n * n) as f64).min(0.5 / ((n * n) as f64 + x * x))
        < 0.5 * crate::consts::ROOT3_DBL_EPS
    {
        result = super::bessel_helpers::besseliv_scaled_asymp_unif_e(n as f64, ax);
        if x < 0.0 && n % 2 == 1 {
            result.val = -result.val
        };
    } else {
        let nhi = 2 + (1.2 / crate::consts::ROOT6_DBL_EPS) as i32;
        let r_ikp1 = super::bessel_helpers::besseliv_scaled_asymp_unif_e(nhi as f64 + 1.0, ax);
        let r_ik = super::bessel_helpers::besseliv_scaled_asymp_unif_e(nhi as f64, ax);
        let mut ikp1 = r_ikp1.val;
        let mut ik = r_ik.val;
        for k in ((n + 1)..(nhi + 1)).rev() {
            let ikm1 = ikp1 + 2.0 * k as f64 / ax * ik;
            ikp1 = ik;
            ik = ikm1;
        }
        result.val = ik;
        result.err = ik * (r_ikp1.err / r_ikp1.val + r_ik.err / r_ik.val);
        if x < 0.0 && n % 2 == 1 {
            result.val = -result.val;
        }
    }

    result
}

/// Compute the modified bessel function of the first kind of order n along
/// with an error estimate.
pub(crate) fn cyl_bessel_in_e(n: i32, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let ax = (x).abs();
    let n = n.abs(); /* I(-n, z) = I(n, z) */
    let in_scaled = cyl_bessel_in_scaled_e(n, ax);

    /* In_scaled is always less than 1,
     * so this overflow check is conservative.
     */
    if ax > crate::consts::LN_DBL_MAX - 1.0 {
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
    } else {
        let ex = (ax).exp();
        result.val = ex * in_scaled.val;
        result.err = ex * in_scaled.err;
        result.err += ax * f64::EPSILON * (result.val).abs();
        if x < 0.0 && n % 2 == 1 {
            result.val = -result.val;
        }
    }

    result
}

/// Compute the modified bessel function of the first kind of fractional order
/// nu scaled by e^-x along with an error estimate.
pub(crate) fn cyl_bessel_iv_scaled_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    /* CHECK_POINTER(result) */

    if x < 0.0 || nu < 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
    } else if x * x < 10.0 * (nu + 1.0) {
        let ex = (-x).exp();
        let b = super::bessel_helpers::bessel_ij_taylor_e(nu, x, 1, 100, f64::EPSILON);
        result.val = b.val * ex;
        result.err = b.err * ex;
        result.err += 2.0 * f64::EPSILON * result.val.abs();
    } else if 0.5 / (nu * nu + x * x) < crate::consts::ROOT3_DBL_EPS {
        return super::bessel_helpers::besseliv_scaled_asymp_unif_e(nu, x);
    } else {
        let nn = (nu + 0.5) as i32;
        let mu = nu - nn as f64; /* -1/2 <= mu <= 1/2 */

        /* obtain K_mu, K_mup1 */
        let (kmu, kmup1, _) = if x < 2.0 {
            super::bessel_helpers::besselk_scaled_temme(mu, x)
        } else {
            super::bessel_helpers::besselk_scaled_steed_temme_cf2(mu, x)
        };

        /* recurse forward to obtain K_num1, K_nu */
        let mut knu = kmu.val;
        let mut knup1 = kmup1.val;

        for n in 0..nn {
            let knum1 = knu;
            knu = knup1;
            knup1 = 2.0 * (mu + n as f64 + 1.0) / x * knu + knum1;
        }

        /* calculate I_{nu+1}/I_nu */
        let inu_ratio = super::bessel_helpers::besseli_cf1_ser(nu, x);

        /* solve for I_nu */
        result.val = 1.0 / (x * (knup1 + inu_ratio.val * knu));
        result.err = f64::EPSILON * (0.5 * nn as f64 + 2.0) * result.val.abs();
    }

    result
}

/// Compute the modified bessel function of the first kind of fractional order
/// nu along with an error estimate.
pub(crate) fn cyl_bessel_iv_e(nu: f64, x: f64) -> SpecFunResult<f64> {
    let b = cyl_bessel_iv_scaled_e(nu, x);
    crate::exp::core::exp_mult_err_e(x, (x * f64::EPSILON).abs(), b.val, b.err)
}
