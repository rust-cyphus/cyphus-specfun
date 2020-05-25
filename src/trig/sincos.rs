use crate::result::SpecFunResult;

/// routine computing sin(pi*x) using a Taylor expansion around the origin
/// and otherwise the library function.
fn sin_pi_taylor(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if 16.0 * x.abs() < 1.0 {
        let y = std::f64::consts::PI * x;
        let a = y * y;
        result.val = y
            * (1f64)
            .mul_add(-a / 110.0, 1.0)
            .mul_add(-a / 72.0, 1.0)
            .mul_add(-a / 42.0, 1.0)
            .mul_add(-a / 20.0, 1.0)
            .mul_add(-a / 6.0, 1.0);
    } else {
        result.val = (std::f64::consts::PI * x).sin();
    }

    result.err = f64::EPSILON * result.val.abs();
    result
}

fn cos_pi_taylor(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if 20.0 * x.abs() < 1.0 {
        let y = std::f64::consts::PI * x;
        let a = y * y;
        result.val = (1f64)
            .mul_add(-a / 90.0, 1.0)
            .mul_add(-a / 56.0, 1.0)
            .mul_add(-a / 30.0, 1.0)
            .mul_add(-a / 12.0, 1.0)
            .mul_add(-0.5 * a, 1.0);
    } else {
        result.val = (std::f64::consts::PI * x).cos();
    }

    result.err = f64::EPSILON * result.val.abs();
    result
}

#[allow(dead_code)]
pub(crate) fn sin_pi_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let twobig = 2.0 / f64::EPSILON;

    let intx = x.trunc() as i64;
    let mut fracx = x.fract();

    if fracx.abs() == 0.0 || (intx as f64).abs() >= twobig {
        return result;
    }

    let mut sign: i32 = if intx % 2 != 0 { -1 } else { 1 };

    if fracx.abs() == 0.5 {
        // probably unnecessary
        if fracx < 0.0 {
            sign *= -1;
        };
        result.val = if sign != 1 { -1.0 } else { 1.0 };
        return result;
    }
    if fracx.abs() > 0.5 {
        sign *= -1;
        fracx += if fracx > 0.0 { -1.0 } else { 1.0 };
    }

    result = if fracx > 0.25 {
        cos_pi_taylor(fracx - 0.5)
    } else if fracx < -0.25 {
        sign *= -1;
        cos_pi_taylor(fracx + 0.5)
    } else {
        sin_pi_taylor(fracx)
    };

    if sign != 1 {
        result.val *= -1.0;
    };

    result
}

#[allow(dead_code)]
pub(crate) fn cos_pi_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let twobig = 2.0 / f64::EPSILON;

    let intx = x.trunc() as i64;
    let mut fracx = x.fract();

    if fracx.abs() == 0.5 {
        return result;
    }
    if (intx as f64).abs() >= twobig {
        result.val = 1.0;
        return result;
    }

    let mut sign: i32 = if intx % 2 != 0 { -1 } else { 1 };

    if fracx.abs() == 0.0 {
        // probably unnecessary
        if fracx < 0.0 {
            sign *= -1;
        };
        result.val = if sign != 1 { -1.0 } else { 1.0 };
        return result;
    }
    if fracx.abs() > 0.5 {
        sign *= -1;
        fracx += if fracx > 0.0 { -1.0 } else { 1.0 };
    }

    result = if fracx > 0.25 {
        sign *= -1;
        sin_pi_taylor(fracx - 0.5)
    } else if fracx < -0.25 {
        sin_pi_taylor(fracx + 0.5)
    } else {
        cos_pi_taylor(fracx)
    };

    if sign != 1 {
        result.val *= -1.0;
    };

    result
}
