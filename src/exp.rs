use crate::consts::{LN_DBL_MAX, LN_DBL_MIN};
use crate::result::{SpecFunCode, SpecFunResult};

/// Exponentiate and multiply by a given factor: y * Exp(x)
pub fn exp_mult_err_e(x: f64, dx: f64, y: f64, dy: f64) -> SpecFunResult<f64> {
    let ay = y.abs();

    if ay < f64::EPSILON {
        SpecFunResult {
            val: 0.0,
            err: (dy * x.exp()).abs(),
            code: SpecFunCode::Success,
        }
    } else if (x < 0.5 * f64::MAX.ln() && x > 0.5 * f64::MIN_POSITIVE.ln())
        && (ay < 0.8 * f64::MAX.sqrt() && ay > 1.0 * f64::MIN_POSITIVE.sqrt())
    {
        let ex = x.exp();
        let val = y * ex;
        let mut err = ex * (dy.abs() + (y * dx).abs());
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        let ly = ay.ln();
        let lnr = x + ly;

        if lnr > f64::MAX.ln() - 0.01 {
            let result = SpecFunResult {
                val: f64::INFINITY,
                err: f64::INFINITY,
                code: SpecFunCode::OverflowErr,
            };
            result.issue_warning("exp_mult_err_e", &[x, dx, y, dy]);
            result
        } else if lnr < f64::MIN_POSITIVE.ln() + 0.01 {
            let result = SpecFunResult {
                val: 0.0,
                err: 0.0,
                code: SpecFunCode::UnderflowErr,
            };
            result.issue_warning("exp_mult_err_e", &[x, dx, y, dy]);
            result
        } else {
            let sy = y.signum();
            let M = x.floor();
            let N = ly.floor();
            let a = x - M;
            let b = ly - N;
            let eMN = (M + N).exp();
            let eab = (a + b).exp();
            let val = sy * eMN * eab;
            let mut err = eMN * eab * 2.0 * f64::EPSILON;
            err += eMN * eab * (dy / y).abs();
            err += eMN * eab * dx.abs();
            SpecFunResult {
                val,
                err,
                code: SpecFunCode::Success,
            }
        }
    }
}

/// Exponentiate a quantity with an associated error.
pub fn exp_err_e(x: f64, dx: f64) -> SpecFunResult<f64> {
    let adx = dx.abs();

    if x + adx > LN_DBL_MAX {
        let result = SpecFunResult {
            val: f64::INFINITY,
            err: f64::INFINITY,
            code: SpecFunCode::OverflowErr,
        };
        result.issue_warning("exp_err_e", &[x, dx]);
        result
    } else if x - adx < LN_DBL_MIN {
        let result = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::UnderflowErr,
        };
        result.issue_warning("exp_err_e", &[x, dx]);
        result
    } else {
        let ex = x.exp();
        let edx = adx.exp();
        let val = ex;
        let mut err = ex * f64::EPSILON.min(edx - 1.0 / edx);
        err += 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}
