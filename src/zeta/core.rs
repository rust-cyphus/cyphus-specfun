use crate::cheb::cheb_eval_e;
use crate::result::{SpecFunCode, SpecFunResult};

use crate::zeta::data::*;

/// Compute the Riemann-Zeta function for s >= 0 (for internal use only)
pub(crate) fn riemann_zeta_sgt0(s: f64) -> SpecFunResult<f64> {
    if s < 1.0 {
        let mut result = cheb_eval_e(2.0 * s - 1.0, &ZETA_XLT1_DATA, -1.0, 1.0);

        result.val *= (s - 1.0).recip();
        result.err = result.err / (s - 1.0).abs() + std::f64::EPSILON * result.val.abs();
        result
    } else if s <= 20.0 {
        let x = (2.0 * s - 21.0) / 19.0;
        let mut result = cheb_eval_e(x, &ZETA_XGT1_DATA, -1.0, 1.0);
        result.val = result.val / (s - 1.0);
        result.err = result.err / (s - 1.0) + std::f64::EPSILON * result.val.abs();
        result
    } else {
        let f2 = 1.0 - 2.0_f64.powf(-s);
        let f3 = 1.0 - 3.0_f64.powf(-s);
        let f5 = 1.0 - 5.0_f64.powf(-s);
        let f7 = 1.0 - 7.0_f64.powf(-s);
        let val = 1.0 / (f2 * f3 * f5 * f7);
        let err = 3.0 * std::f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

/// Compute the Riemann-Zeta function for s < 0 (for internal use only)
pub(crate) fn riemann_zeta1ms_slt0(s: f64) -> SpecFunResult<f64> {
    if s > -19.0 {
        let x = (-19.0 - 2.0 * s) / 19.0;
        let c = cheb_eval_e(x, &ZETA_XGT1_DATA, -1.0, 1.0);
        let val = c.val / (-s);
        let err = c.err / (-s) + std::f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        let f2 = 1.0 - 2.0_f64.powf(-(1.0 - s));
        let f3 = 1.0 - 3.0_f64.powf(-(1.0 - s));
        let f5 = 1.0 - 5.0_f64.powf(-(1.0 - s));
        let f7 = 1.0 - 7.0_f64.powf(-(1.0 - s));
        let val = 1.0 / (f2 * f3 * f5 * f7);
        let err = 3.0 * std::f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

/// Compute the Riemann-Zeta function for 5 < s < 15 (for internal use only)
pub(crate) fn riemann_zeta_minus_1_intermediate_s(s: f64) -> SpecFunResult<f64> {
    let t = (s - 10.0) / 5.0;
    let c = cheb_eval_e(t, &ZETAM1_INTER_DATA, -1.0, 1.0);
    let val = c.val.exp() + 2.0_f64.powf(-s);
    let err = (c.err + 2.0 * std::f64::EPSILON) * val;
    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}

/// Compute the Riemann-Zeta function for large, positive s
pub(crate) fn riemann_zeta_minus1_large_s(s: f64) -> SpecFunResult<f64> {
    let a = 2.0_f64.powf(-s);
    let b = 3.0_f64.powf(-s);
    let c = 5.0_f64.powf(-s);
    let d = 7.0_f64.powf(-s);
    let e = 11.0_f64.powf(-s);
    let f = 13.0_f64.powf(-s);
    let t1 = a + b + c + d + e + f;
    let t2 = a * (b + c + d + e + f) + b * (c + d + e + f) + c * (d + e + f) + d * (e + f) + e * f;

    let numt = t1 - t2 /* + t3 - t4 + t5 - t6 */;
    let zeta = 1.0 / ((1.0 - a) * (1.0 - b) * (1.0 - c) * (1.0 - d) * (1.0 - e) * (1.0 - f));
    let val = numt * zeta;
    let err = (15.0 / s + 1.0) * 6.0 * std::f64::EPSILON * val;

    SpecFunResult {
        val,
        err,
        code: SpecFunCode::Success,
    }
}
