use crate::consts::{ROOT5_DBL_EPS, SQRT_DLB_EPS};
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
