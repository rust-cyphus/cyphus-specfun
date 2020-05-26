use crate::cheb::cheb_eval_e;
use crate::consts::{ROOT5_DBL_EPS, SQRT_DLB_EPS};
use crate::result::{SpecFunCode, SpecFunResult};

use super::data::*;
use super::helpers::{bessel_cos_pi4_e, bessel_sin_pi4_e};

/// Compute the Bessel function of the first kind of order 0 with error
/// estimate.
///
/// # Example
/// ```
/// assert!((besselj0_e(1.0).val - 0.76519768655796655).abs() < 1e-10);
/// ```
#[allow(dead_code)]
pub fn besselj0_e(x: f64) -> SpecFunResult<f64> {
    let y = x.abs();
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if y <= 2.0 * SQRT_DLB_EPS {
        result.val = 1.0;
        result.err = y * y;
        result
    } else if y <= 4.0 {
        cheb_eval_e(0.125 * y * y - 1.0, &BJ0_DATA, -1.0, 1.0)
    } else {
        let z = 32.0 / (y * y) - 1.0;
        let ca = cheb_eval_e(z, &BM0_DATA, -1.0, 1.0);
        let ct = cheb_eval_e(z, &BTH0_DATA, -1.0, 1.0);
        let cp = bessel_cos_pi4_e(y, ct.val / y);
        let sqrty = y.sqrt();
        let ampl = (0.75 + ca.val) / sqrty;
        result.val = ampl * cp.val;
        result.err = cp.val.abs() * ca.err / sqrty + ampl.abs() * cp.err;
        result.err += f64::EPSILON * result.val.abs();

        result
    }
}

/// Compute the Bessel function of the first kind of order 1 with error
/// estimate.
///
/// # Example
/// ```
/// assert!((besselj1_e(1.0).val - 0.44005058574493352).abs() < 1e-10);
/// ```
#[allow(dead_code)]
pub fn besselj1_e(x: f64) -> SpecFunResult<f64> {
    let y = x.abs();
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if y <= 2.0 * SQRT_DLB_EPS {
        result.val = 1.0;
        result.err = y * y;
        result
    } else if y < 2.0 * std::f64::consts::SQRT_2 * SQRT_DLB_EPS {
        result.val = 0.5 * x;
        result.err = 0.0;
        result
    } else if y <= 4.0 {
        let c = cheb_eval_e(0.125 * y * y - 1.0, &BJ1_DATA, -1.0, 1.0);
        result.val = x * (0.25 + c.val);
        result.err = (x * c.err).abs();
        result
    } else {
        let z = 32.0 / (y * y) - 1.0;
        let ca = cheb_eval_e(z, &BM1_DATA, -1.0, 1.0);
        let ct = cheb_eval_e(z, &BTH1_DATA, -1.0, 1.0);
        let sp = bessel_sin_pi4_e(y, ct.val / y);
        let sqrty = y.sqrt();
        let ampl = (0.75 + ca.val) / sqrty;
        result.val = x.signum() * ampl * sp.val;
        result.err = sp.val.abs() * ca.err / sqrty + ampl.abs() * sp.err;
        result.err += f64::EPSILON * result.val.abs();

        result
    }
}

/// Compute the Bessel function of the first kind of order n with error
/// estimate.
///
/// # Example
/// ```
/// assert!((besseljn_e(4, 1.0).val - 0.0024766389641099550).abs() < 1e-10);
/// ```
#[allow(dead_code)]
pub fn besseljn_e(nn: i32, xx: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::default();

    let mut sign = 1.0f64;
    let mut n = nn;
    let mut x = xx;

    if n < 0 {
        // reduce to case n >= 0
        n = -n;
        if n % 2 == 1 {
            sign *= -1.0;
        }
    }

    if xx < 0.0 {
        // reduce to case x >= 0
        x = -xx;
        if n % 2 == 1 {
            sign *= -1.0;
        }
    }

    if n == 0 {
        let j0 = besselj0_e(x);
        result.val = sign as f64 * j0.val;
        result.err = j0.err;
    } else if n == 1 {
        let j1 = besselj1_e(x);
        result.val = sign as f64 * j1.val;
        result.err = j1.err;
    } else if x == 0.0 {
        return result;
    } else if x * x < 10.0 * (n + 1) as f64 * ROOT5_DBL_EPS {
        let b = super::helpers::bessel_ij_taylor_e(n as f64, x, -1, 50, f64::EPSILON);
        result.val = b.val * sign as f64;
        result.err = b.err;
        result.err += f64::EPSILON * result.val.abs();
        return result;
    } else if x * crate::consts::ROOT4_DBL_EPS > ((n * n) as f64 + 1.0) {
        result = super::helpers::besseljv_asympx_e(n as f64, x);
        result.val *= sign as f64;
        return result;
    } else if n > 50 {
        result = super::olver::besseljv_asymp_olver_e(n as f64, x);
        result.val *= sign as f64;
        return result;
    } else if x > 1000.0 {
        // We need this to avoid feeding large x to CF1; note that due to
        // the above check, we know that n <= 50.
        result = super::helpers::besseljv_asympx_e(n as f64, x);
        result.val *= sign as f64;
        return result;
    } else {
        let (ratio, _) = super::helpers::besselj_cf1(n as f64, x);

        /* backward recurrence */
        let mut jkp1 = crate::consts::SQRT_DBL_MIN * ratio.val;
        let mut jk = crate::consts::SQRT_DBL_MIN;

        for k in (1..(n + 1)).rev() {
            let jkm1 = 2.0 * k as f64 / x * jk - jkp1;
            jkp1 = jk;
            jk = jkm1;
        }

        let (ans, err) = if jkp1.abs() > jk.abs() {
            let b1 = besselj1_e(x);
            (
                b1.val / jkp1 * crate::consts::SQRT_DBL_MIN,
                b1.err / jkp1 * crate::consts::SQRT_DBL_MIN,
            )
        } else {
            let b0 = besselj0_e(x);
            (
                b0.val / jk * crate::consts::SQRT_DBL_MIN,
                b0.err / jk * crate::consts::SQRT_DBL_MIN,
            )
        };

        result.val = sign * ans;
        result.err = err.abs();
    }

    result
}
