use crate::result::{SpecFunCode, SpecFunResult};

pub(crate) fn powi_e(x: f64, n: i32) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let mut value = 1.0;
    let mut count = 0;
    let mut nn = n;
    let mut xx = x;

    if nn < 0 {
        nn = -nn;

        if xx == 0.0 {
            let u = 1.0 / xx;
            result.val = if nn % 2 == 1 { u } else { u * u };
            result.err = f64::INFINITY;
            result.code = SpecFunCode::OverflowErr;
        }
        xx = 1.0 / xx;
    }

    // repeated squaring method. Returns 0.0^0 = 1.0, so continuous in x.
    loop {
        if nn % 2 == 1 {
            value *= xx;
        }
        nn >>= 1;
        xx *= xx;
        count += 1;
        if nn == 0 {
            break;
        }
    }

    result.val = value;
    result.err = 2.0 * f64::EPSILON * (1 + count) as f64 * value.abs();

    result
}
