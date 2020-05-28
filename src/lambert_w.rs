use crate::result::{SpecFunCode, SpecFunResult};
use num::Float;;


pub trait LambertW{
    /// Compute the Lambert W-function on the principal branch along with an
    /// error estimate.
    fn lambert_w0_e(&self) -> SpecFunResult<Self> where Self:Float;
    /// Compute the Lambert W-function on the principal branch.
    fn lambert_w0(&self) -> Self;
    /// Compute the Lambert W-function on the second real-valued branch along
    /// with an error estimate.
    fn lambert_wm1_e(&self) -> SpecFunResult<Self>where Self:Float;
    /// Compute the Lambert W-function on the second real-valued branch.
    fn lambert_wm1(&self) -> Self;
}

impl LambertW for f64{
    fn lambert_w0_e(&self) -> SpecFunResult<Self>{
        lambert_w0_e(*self)
    }
    fn lambert_w0(&self) -> Self{
        lambert_w0_e(*self).val
    }
    fn lambert_wm1_e(&self) -> SpecFunResult<Self>{
        lambert_wm1_e(*self)
    }
    fn lambert_wm1(&self) -> Self{
        lambert_wm1_e(*self).val
    }
}


/// Perform a Halley-iteration to determine root of x = w * exp(w)
#[allow(dead_code)]
fn halley_iteration(x: f64, w_initial: f64, max_iters: usize) -> SpecFunResult<f64> {
    let mut w = w_initial;
    let mut result = SpecFunResult::<f64>::default();

    for _ in 0..max_iters {
        let mut tol = 0.0;
        let mut e = w.exp();
        let mut p = w + 1.0;
        let mut t = w * e - x;

        t = if w > 0.0 {
            (t / p) / e // Newton iteration
        } else {
            t / (e * p - 0.5 * (p + 1.0) * t / p) // Halley iteration
        };

        w -= t;

        tol = 10.0 * f64::EPSILON * w.abs().max(1.0 / (p.abs() * e));

        if t.abs() < tol {
            result.val = w;
            result.err = 2.0 * tol;
            return result;
        }
    }

    // should never get here
    result.val = w;
    result.err = w.abs();
    result.code = SpecFunCode::MaxIterErr;
    result
}

/* series which appears for q near zero;
 * only the argument is different for the different branches
 */
fn series_eval(r: f64) -> f64 {
    27.029044799010561650f64
        .mul_add(r, -18.100697012472442755)
        .mul_add(r, 12.250753501314460424)
        .mul_add(r, -8.401032217523977370984161688514)
        .mul_add(r, 5.858023729874774148815053846119)
        .mul_add(r, -4.175335600258177138854984177460)
        .mul_add(r, 3.066858901050631912893148922704)
        .mul_add(r, -2.353551201881614516821543561516)
        .mul_add(r, 1.936631114492359755363277457668)
        .mul_add(r, -1.812187885639363490240191647568)
        .mul_add(r, 2.331643981597124203363536062168)
        .mul_add(r, -1.0)
}

/// Compute the Lambert W-function on the principal branch along with an error
/// estimate.
pub(crate) fn lambert_w0_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    let q = x + 1.0 / std::f64::consts::E;

    if x == 0.0 {
        result.val = 0.0;
        result.err = 0.0;
        result
    } else if q < 0.0 {
        // Strictly speaking this is an error. But because of the
        // arithmetic operation connecting x and q, I am a little
        // lenient in case of some epsilon overshoot. The following
        // answer is quite accurate in that case. Anyway, we have
        // to return GSL_EDOM.
        result.val = -1.0;
        result.err = (-q).sqrt();
        result.code = SpecFunCode::DomainErr;
        result
    } else if q == 0.0 {
        result.val = -1.0;
        result.err = f64::EPSILON; // cannot error is zero, maybe q == 0 by "accident"
        result
    } else if q < 1.0e-03 {
        /* series near -1/E in sqrt(q) */
        let r = q.sqrt();
        result.val = series_eval(r);
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        let w = if x < 1.0 {
            // obtain initial approximation from series near x=0;
            // no need for extra care, since the Halley iteration
            // converges nicely on this branch
            let p = (2.0 * std::f64::consts::E * q).sqrt();
            (11.0 / 72.0 as f64)
                .mul_add(p, -1.0 / 3.0)
                .mul_add(p, 1.0)
                .mul_add(p, -1.0)
        } else {
            // obtain initial approximation from rough asymptotic
            if x > 3.0 {
                x.ln() - x.ln().ln()
            } else {
                x.ln()
            }
        };
        halley_iteration(x, w, 10)
    }
}

/// Compute the Lambert W-function on the second real-valued branch along with
/// an error estimate.
pub(crate) fn lambert_wm1_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x > 0.0 {
        lambert_w0_e(x)
    } else if x == 0.0 {
        result.val = 0.0;
        result.err = 0.0;
        result
    } else {
        let q = x + 1.0 / std::f64::consts::E;
        let w;

        if q < 0.0 {
            // As in the W0 branch above, return some reasonable answer anyway.
            result.val = -1.0;
            result.err = (-q).sqrt();
            result.code = SpecFunCode::DomainErr;
            return result;
        }

        if x < -1.0e-6 {
            // Obtain initial approximation from series about q = 0,
            // as long as we're not very close to x = 0.
            // Use full series and try to bail out if q is too small,
            // since the Halley iteration has bad convergence properties
            // in finite arithmetic for q very small, because the
            // increment alternates and p is near zero.
            let r = -(q.sqrt());
            w = series_eval(r);
            if q < 3.0e-3 {
                /* this approximation is good enough */
                result.val = w;
                result.err = 5.0 * f64::EPSILON * w.abs();
                return result;
            }
        } else {
            // Obtain initial approximation from asymptotic near zero.
            let l1 = (-x).ln();
            let l2 = (-l1).ln();
            w = l1 - l2 + l2 / l1;
        }
        halley_iteration(x, w, 32)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_check_result_and_code;
    use crate::test_utils::*;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const TOL1: f64 = 16.0 * f64::EPSILON;

    #[test]
    fn test_lambert_w0_e() {
        test_check_result_and_code!(lambert_w0_e, (0.0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            lambert_w0_e,
            (1.0),
            0.567143290409783872999969,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (2.0),
            0.852605502013725491346472,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (20.0),
            2.205003278024059970493066,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (1000.0),
            5.24960285240159622712606,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (1.0e+6),
            11.38335808614005262200016,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (1.0e+12),
            24.43500440493491313826305,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (1.0e+308),
            702.641362034106812081125,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            lambert_w0_e,
            (1.6849341956993852953416990),
            0.775706963944252869680440,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (-1.0 / std::f64::consts::E - f64::EPSILON),
            -1.0,
            TOL0,
            SpecFunCode::DomainErr
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (-1.0 / std::f64::consts::E + 1.0 / (1024.0 * 1024.0 * 1024.0)),
            -0.999928845560308370714970,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (-1.0 / std::f64::consts::E + 1.0 / (1024.0 * 1024.0)),
            -0.997724730359774141620354,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (-1.0 / std::f64::consts::E + 1.0 / 512.0),
            -0.900335676696088773044678,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_w0_e,
            (-1.0 / std::f64::consts::E + 0.25),
            -0.1349044682661213545487599,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_lambert_wm1_e() {
        test_check_result_and_code!(lambert_wm1_e, (0.0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            lambert_wm1_e,
            (1.0),
            0.567143290409783872999969,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_wm1_e,
            (2.0),
            0.852605502013725491346472,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_wm1_e,
            (20.0),
            2.205003278024059970493066,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            lambert_wm1_e,
            (-1.0 / std::f64::consts::E - f64::EPSILON),
            -1.0,
            TOL0,
            SpecFunCode::DomainErr
        );
        test_check_result_and_code!(
            lambert_wm1_e,
            (-1.0 / std::f64::consts::E + 1.0 / (1024.0 * 1024.0 * 1024.0)),
            -1.000071157815154608049055,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_wm1_e,
            (-1.0 / std::f64::consts::E + 1.0 / (1024.0 * 1024.0)),
            -1.002278726118593023934693,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_wm1_e,
            (-1.0 / std::f64::consts::E + 1.0 / 512.0),
            -1.106761200865743124599130,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_wm1_e,
            (-1.0 / std::f64::consts::E + 1.0 / 64.0),
            -1.324240940341812125489772,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lambert_wm1_e,
            (-1.0 / std::f64::consts::E + 0.25),
            -3.345798131120112,
            TOL1,
            SpecFunCode::Success
        );
    }
}
