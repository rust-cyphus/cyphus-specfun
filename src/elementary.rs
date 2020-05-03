use crate::consts::SQRT_DLB_EPS;
use crate::result::{SpecFunCode, SpecFunResult};

/// Multiply x and y and return result along with error estimate
pub(crate) fn multiply_e(x: f64, y: f64) -> SpecFunResult<f64> {
    let ax = x.abs();
    let ay = y.abs();

    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if x == 0.0 || y == 0.0 {
        // It is necessary to eliminate this immediately.
        result
    } else if (ax <= 1.0 && ay >= 1.0) || (ay <= 1.0 && ax >= 1.0) {
        // Straddling 1.0 is always safe.
        result.val = x * y;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        let f = 1.0 - 2.0 * f64::EPSILON;
        let min = ax.min(ay);
        let max = ax.max(ay);

        if max < 0.9 * SQRT_DLB_EPS || min < (f * f64::MAX) / max {
            result.val = x * y;
            result.err = 2.0 * f64::EPSILON * result.val.abs();
            if result.val.abs() < f64::MIN {
                result.code = SpecFunCode::UnderflowErr;
                result.issue_warning("multiply_e", &[x, y]);
            }
        } else {
            result.val = f64::INFINITY;
            result.err = f64::INFINITY;
            result.code = SpecFunCode::OverflowErr;
            result.issue_warning("multiply_e", &[x, y]);
        }
        result
    }
}

/// Multiply x and y and return result along with error estimate given errors
/// dx and dy in x and y respectively.
pub(crate) fn multiply_err_e(x: f64, dx: f64, y: f64, dy: f64) -> SpecFunResult<f64> {
    let mut result = multiply_e(x, y);
    result.err += (dx * y).abs() + (dy * x).abs();
    result
}
