use crate::airy::airy_data::*;
use crate::result::{SpecFunCode, SpecFunResult};
use crate::trig::{cos_err_e, sin_err_e};

// Airy modulus and phase for x < -1
fn airy_mod_phase(x: f64) -> (SpecFunResult<f64>, SpecFunResult<f64>) {
    let result_m: SpecFunResult<f64>;
    let result_p: SpecFunResult<f64>;
    let mut result_mod = SpecFunResult::<f64>::default();
    let mut result_phase = SpecFunResult::<f64>::default();

    // cheb_eval_e\((.+?),\s(.+?), -1.0, 1.0\);
    // %2.eval($1);

    if x < -2.0 {
        let z = 16.0 / (x * x * x) + 1.0;
        result_m = (*AM21_CHEB).eval(z);
        result_p = (*ATH1_CHEB).eval(z);
    } else if x <= -1.0 {
        let z = (16.0 / (x * x * x) + 9.0) / 7.0;
        result_m = (*AM22_CHEB).eval(z);
        result_p = (*ATH2_CHEB).eval(z);
    } else {
        result_mod.val = 0.0;
        result_mod.err = 0.0;
        result_mod.code = SpecFunCode::DomainErr;
        result_phase.val = 0.0;
        result_phase.err = 0.0;
        result_phase.code = SpecFunCode::DomainErr;
        result_mod.issue_warning("airy_mod_phase", &[x]);
        return (result_mod, result_phase);
    }

    let m = 0.3125 + result_m.val;
    let p = -0.625 + result_p.val;
    let sqx = (-x).sqrt();

    result_mod.val = (m / sqx).sqrt();
    result_mod.err = result_mod.val.abs() * (f64::EPSILON + (result_m.err / result_m.val).abs());
    result_phase.val = std::f64::consts::FRAC_PI_4 - x * sqx * p;
    result_phase.err =
        result_phase.val.abs() * (f64::EPSILON + (result_p.err / result_p.val).abs());
    (result_mod, result_phase)
}

// assumes x >= 1.0
#[inline]
fn airy_aie(x: f64) -> SpecFunResult<f64> {
    let sqx = x.sqrt();
    let z = 2.0 / (x * sqx) - 1.0;
    let y = sqx.sqrt();
    let result_c = (*AIP_CHEB).eval(z);
    let mut result = SpecFunResult::<f64>::default();
    result.val = (0.28125 + result_c.val) / y;
    result.err = result_c.err / y + f64::EPSILON * result.val.abs();
    result
}

// assumes x >= 2.0
fn airy_bie(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    let atr: f64 = 8.750_690_570_848_434;
    let btr: f64 = -2.093_836_321_356_054_2;

    if x < 4.0 {
        let sqx = x.sqrt();
        let z = atr / (x * sqx) + btr;
        let y = sqx.sqrt();
        let result_c = &(*BIP_CHEB).eval(z);
        result.val = (0.625 + result_c.val) / y;
        result.err = result_c.err / y + f64::EPSILON * result.val.abs();
    } else {
        let sqx = x.sqrt();
        let z = 16.0 / (x * sqx) - 1.0;
        let y = sqx.sqrt();
        let result_c = (*BIP2_CHEB).eval(z);
        result.val = (0.625 + result_c.val) / y;
        result.err = result_c.err / y + f64::EPSILON * result.val.abs();
    }
    result
}

/// Compute the Airy function Ai(x).
pub(crate) fn airy_ai_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x < -1.0 {
        let (rmod, theta) = airy_mod_phase(x);
        let cos_result = cos_err_e(theta.val, theta.err);
        result.val = rmod.val * cos_result.val;
        result.err = (rmod.val * cos_result.err).abs() + (cos_result.val * rmod.err).abs();
        result.err += f64::EPSILON * result.val.abs();
        result
    } else if x <= 1.0 {
        let z = x * x * x;
        let result_c0 = (*AIF_CHEB).eval(z);
        let result_c1 = (*AIG_CHEB).eval(z);
        result.val = 0.375 + (result_c0.val - x * (0.25 + result_c1.val));
        result.err = result_c0.err + (x * result_c1.err).abs();
        result.err += f64::EPSILON * result.val.abs();
        result
    } else {
        let x32 = x * x.sqrt();
        let s = (-2.0 * x32 / 3.0).exp();
        let result_aie = airy_aie(x);
        result.val = result_aie.val * s;
        result.err = result_aie.err * s + result.val * x32 * f64::EPSILON;
        result.err += f64::EPSILON * result.val.abs();
        if result.val.abs() < f64::MIN_POSITIVE {
            result.code = SpecFunCode::UnderflowErr;
        }
        result
    }
}

/// Compute the Airy function Ai(x) scaled by S(x), where S(x) is
/// exp(2x^{3/2}/3) for x > 0 and 1 for x < 0;
pub(crate) fn airy_ai_scaled_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x < -1.0 {
        let (rmod, theta) = airy_mod_phase(x);
        let cos_result = cos_err_e(theta.val, theta.err);
        result.val = rmod.val * cos_result.val;
        result.err = (rmod.val * cos_result.err).abs() + (cos_result.val * rmod.err).abs();
        result.err += f64::EPSILON * result.val.abs();
        result
    } else if x <= 1.0 {
        let z = x * x * x;
        let result_c0 = (*AIF_CHEB).eval(z);
        let result_c1 = (*AIG_CHEB).eval(z);
        result.val = 0.375 + (result_c0.val - x * (0.25 + result_c1.val));
        result.err = result_c0.err + (x * result_c1.err).abs();
        result.err += f64::EPSILON * result.val.abs();

        if x > 0.0 {
            let scale = (2.0 / 3.0 * z.sqrt()).exp();
            result.val *= scale;
            result.err *= scale;
        }

        result
    } else {
        airy_aie(x)
    }
}

/// Compute the Airy function Bi(x).
pub(crate) fn airy_bi_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x < -1.0 {
        let (rmod, theta) = airy_mod_phase(x);
        let sin_result = sin_err_e(theta.val, theta.err);
        result.val = rmod.val * sin_result.val;
        result.err = (rmod.val * sin_result.err).abs() + (sin_result.val * rmod.err).abs();
        result.err += f64::EPSILON * result.val.abs();
    } else if x <= 1.0 {
        let z = x * x * x;
        let result_c0 = (*BIF_CHEB).eval(z);
        let result_c1 = (*BIG_CHEB).eval(z);
        result.val = 0.625 + (result_c0.val + x * (0.4375 + result_c1.val));
        result.err = result_c0.err + (x * result_c1.err).abs();
        result.err += f64::EPSILON * result.val.abs();
    } else if x <= 2.0 {
        let z = (2.0 * x * x * x - 9.0) / 7.0;
        let result_c0 = (*BIF2_CHEB).eval(z);
        let result_c1 = (*BIG2_CHEB).eval(z);
        result.val = 1.125 + result_c0.val + x * (0.625 + result_c1.val);
        result.err = result_c0.err + (x * result_c1.err).abs();
        result.err += f64::EPSILON * result.val.abs();
    } else {
        let y = 2.0 * x * x.sqrt() / 3.0;
        let s = y.exp();

        if y > crate::consts::LN_DBL_MAX - 1.0 {
            result.val = f64::INFINITY;
            result.err = f64::INFINITY;
            result.code = SpecFunCode::UnderflowErr;
        } else {
            let result_bie = airy_bie(x);
            result.val = result_bie.val * s;
            result.err = result_bie.err * s + (1.5 * y * (f64::EPSILON * result.val).abs());
            result.err += f64::EPSILON * result.val.abs();
        }
    }
    result
}

/// Compute the Airy function Bi(x) scaled by S(x), where S(x) is
/// exp(-2x^{3/2}/3) for x > 0 and 1 for x < 0;
pub(crate) fn airy_bi_scaled_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x < -1.0 {
        let (rmod, theta) = airy_mod_phase(x);
        let sin_result = sin_err_e(theta.val, theta.err);
        result.val = rmod.val * sin_result.val;
        result.err = (rmod.val * sin_result.err).abs() + (sin_result.val * rmod.err).abs();
        result.err += f64::EPSILON * result.val.abs();
    } else if x <= 1.0 {
        let z = x * x * x;
        let result_c0 = (*BIF_CHEB).eval(z);
        let result_c1 = (*BIG_CHEB).eval(z);
        result.val = 0.625 + (result_c0.val + x * (0.4375 + result_c1.val));
        result.err = result_c0.err + (x * result_c1.err).abs();
        result.err += f64::EPSILON * result.val.abs();
        if x > 0.0 {
            let scale = (-2.0 / 3.0 * z.sqrt()).exp();
            result.val *= scale;
            result.err *= scale;
        }
    } else if x <= 2.0 {
        let x3 = x * x * x;
        let z = (2.0 * x3 - 9.0) / 7.0;
        let s = (-2.0 / 3.0 * x3.sqrt()).exp();
        let result_c0 = (*BIF2_CHEB).eval(z);
        let result_c1 = (*BIG2_CHEB).eval(z);
        result.val = s * (1.125 + result_c0.val + x * (0.625 + result_c1.val));
        result.err = s * (result_c0.err + (x * result_c1.err).abs());
        result.err += f64::EPSILON * result.val.abs();
    } else {
        return airy_bie(x);
    }
    result
}
