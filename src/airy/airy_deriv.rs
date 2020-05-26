use super::airy_data::*;
use crate::exp::core::exp_mult_err_e;
use crate::result::{SpecFunCode, SpecFunResult};

pub(crate) fn airy_deriv_mod_phase(x: f64) -> (SpecFunResult<f64>, SpecFunResult<f64>) {
    let pi34 = 2.356_194_490_192_345;
    let mut res_mod = SpecFunResult::<f64>::default();
    let mut res_phase = SpecFunResult::<f64>::default();

    if x > -1.0 {
        res_mod.val = 0.0;
        res_mod.err = 0.0;
        res_phase.val = 0.0;
        res_phase.err = 0.0;
        res_mod.code = SpecFunCode::DomainErr;
        res_phase.code = SpecFunCode::DomainErr;
    } else {
        let (result_a, result_p) = if x <= -4.0 {
            let z = 128.0 / (x * x * x) + 1.0;
            ((*D_AN20_CHEB).eval(z), (*D_APH0_CHEB).eval(z))
        } else if x <= -2.0 {
            let z = (128.0 / (x * x * x) + 9.0) / 7.0;
            ((*D_AN21_CHEB).eval(z), (*D_APH1_CHEB).eval(z))
        } else {
            let z = (16.0 / (x * x * x) + 9.0) / 7.0;
            ((*D_AN22_CHEB).eval(z), (*D_APH2_CHEB).eval(z))
        };
        let a = 0.3125 + result_a.val;
        let p = -0.625 + result_p.val;

        let sqx = (-x).sqrt();

        res_mod.val = (a * sqx).sqrt();
        res_mod.err = res_mod.val.abs() * (f64::EPSILON + (result_a.err / result_a.val).abs());

        res_phase.val = pi34 - x * sqx * p;
        res_phase.err = res_phase.val.abs() * (f64::EPSILON + (result_p.err / result_p.val).abs());
    }
    (res_mod, res_phase)
}

pub(crate) fn airy_ai_deriv_scaled_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x < -1.0 {
        let (a, p) = airy_deriv_mod_phase(x);
        let c = p.val.cos();
        result.val = a.val * c;
        result.err = (result.val * p.err).abs() + (c * a.err).abs();
        result.err += f64::EPSILON * result.val.abs();
    } else if x <= 1.0 {
        let x2 = x * x;
        let x3 = x * x2;
        let result_c0 = (*D_AIF_CHEB).eval(x3);
        let result_c1 = (*D_AIG_CHEB).eval(x3);

        result.val = (x2 * (0.125 + result_c0.val) - result_c1.val) - 0.25;
        result.err = (x2 * result_c0.val).abs() + result_c1.err;
        result.err += f64::EPSILON * result.val.abs();
        if x > crate::consts::ROOT3_DBL_EPS.powi(2) {
            /* scale only if x is positive */
            let s = (2.0 * x * x.sqrt() / 3.0).exp();
            result.val *= s;
            result.err *= s;
        }
    } else if x <= 4.0 {
        let sqrtx = x.sqrt();
        let z = (16.0 / (x * sqrtx) - 9.0) / 7.0;
        let s = sqrtx.sqrt();
        let result_c0 = (*D_AIP1_CHEB).eval(z);
        result.val = -(0.28125 + result_c0.val) * s;
        result.err = result_c0.err * s;
        result.err += f64::EPSILON * result.val.abs();
    } else {
        let sqrtx = x.sqrt();
        let z = 16.0 / (x * sqrtx) - 1.0;
        let s = sqrtx.sqrt();
        let result_c0 = (*D_AIP2_CHEB).eval(z);
        result.val = -(0.28125 + result_c0.val) * s;
        result.err = result_c0.err * s;
        result.err += f64::EPSILON * result.val.abs();
    }

    result
}

pub(crate) fn airy_ai_deriv_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();
    if x < -1.0 {
        let (a, p) = airy_deriv_mod_phase(x);
        let c = p.val.cos();
        result.val = a.val * c;
        result.err = (result.val * p.err).abs() + (c * a.err).abs();
        result.err += f64::EPSILON * result.val.abs();
        result
    } else if x < 1.0 {
        let x3 = x * x * x;
        let result_c1 = (*D_AIF_CHEB).eval(x3);
        let result_c2 = (*D_AIG_CHEB).eval(x3);
        result.val = (x * x * (0.125 + result_c1.val) - result_c2.val) - 0.25;
        result.err = (x * x * result_c1.err).abs() + result_c2.err;
        result.err += f64::EPSILON * result.val.abs();
        result
    } else if x * x * x < 9.0 / 4.0 * crate::consts::LN_DBL_MIN.powi(2) {
        let arg = -2.0 * x * x.sqrt() / 3.0;
        let result_aps = airy_ai_deriv_scaled_e(x);
        exp_mult_err_e(
            arg,
            1.5 * (arg * f64::EPSILON).abs(),
            result_aps.val,
            result_aps.err,
        )
    } else {
        result.code = SpecFunCode::UnderflowErr;
        result
    }
}

pub(crate) fn airy_bi_deriv_scaled_e(x: f64) -> SpecFunResult<f64> {
    let atr: f64 = 8.750_690_570_848_434; // 16./(sqrt(8)-1)
    let btr: f64 = -2.093_836_321_356_054_2; // -(sqrt(8)+1)/(sqrt(8)-1)

    let mut result = SpecFunResult::<f64>::default();

    if x < -1.0 {
        let (a, p) = airy_deriv_mod_phase(x);
        let s = p.val.sin();
        result.val = a.val * s;
        result.err = (result.val * p.err).abs() + (s * a.err).abs();
        result.err += f64::EPSILON * (result.val).abs();
    } else if x < 1.0 {
        let x2 = x * x;
        let x3 = x * x2;
        let result_c1 = (*D_BIF_CHEB).eval(x3);
        let result_c2 = (*D_BIG_CHEB).eval(x3);
        result.val = x2 * (result_c1.val + 0.25) + result_c2.val + 0.5;
        result.err = x2 * result_c1.err + result_c2.err;
        result.err += f64::EPSILON * result.val.abs();

        if x > crate::consts::ROOT3_DBL_EPS * crate::consts::ROOT3_DBL_EPS {
            // scale only if x is positive
            let s = (-2.0 * x * x.sqrt() / 3.0).exp();
            result.val *= s;
            result.err *= s;
        }
    } else if x < 2.0 {
        let z = (2.0 * x * x * x - 9.0) / 7.0;
        let s = (-2.0 * x * x.sqrt() / 3.0).exp();
        let result_c0 = (*D_BIF2_CHEB).eval(z);
        let result_c1 = (*D_BIG2_CHEB).eval(z);
        result.val = s * (x * x * (0.25 + result_c0.val) + 0.5 + result_c1.val);
        result.err = s * (x * x * result_c0.err + result_c1.err);
        result.err += f64::EPSILON * result.val.abs();
    } else if x < 4.0 {
        let sqrtx = x.sqrt();
        let z = atr / (x * sqrtx) + btr;
        let s = sqrtx.sqrt();
        let result_c0 = (*D_BIP1_CHEB).eval(z);
        result.val = s * (0.625 + result_c0.val);
        result.err = s * result_c0.err;
        result.err += f64::EPSILON * (result.val).abs();
    } else {
        let sqrtx = x.sqrt();
        let z = 16.0 / (x * sqrtx) - 1.0;
        let s = sqrtx.sqrt();
        let result_c0 = (*D_BIP2_CHEB).eval(z);
        result.val = s * (0.625 + result_c0.val);
        result.err = s * result_c0.err;
        result.err += f64::EPSILON * result.val.abs();
    }
    result
}

pub(crate) fn airy_bi_deriv_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if x < -1.0 {
        let (a, p) = airy_deriv_mod_phase(x);
        let s = p.val.sin();
        result.val = a.val * s;
        result.err = (result.val * p.err).abs() + (s * a.err).abs();
        result.err += f64::EPSILON * result.val.abs();
        result
    } else if x < 1.0 {
        let x2 = x * x;
        let x3 = x * x2;
        let result_c1 = (*D_BIF_CHEB).eval(x3);
        let result_c2 = (*D_BIG_CHEB).eval(x3);
        result.val = x2 * (result_c1.val + 0.25) + result_c2.val + 0.5;
        result.err = x2 * result_c1.err + result_c2.err;
        result.err += f64::EPSILON * result.val.abs();
        result
    } else if x < 2.0 {
        let z = (2.0 * x * x * x - 9.0) / 7.0;
        let result_c1 = (*D_BIF2_CHEB).eval(z);
        let result_c2 = (*D_BIG2_CHEB).eval(z);
        result.val = x * x * (result_c1.val + 0.25) + result_c2.val + 0.5;
        result.err = x * x * result_c1.err + result_c2.err;
        result.err += f64::EPSILON * result.val.abs();
        result
    } else if x < crate::consts::ROOT3_DBL_MAX * crate::consts::ROOT3_DBL_MAX {
        let arg = 2.0 * (x * x.sqrt() / 3.0);
        let result_bps = airy_bi_deriv_scaled_e(x);
        exp_mult_err_e(
            arg,
            1.5 * (arg * f64::EPSILON).abs(),
            result_bps.val,
            result_bps.err,
        )
    } else {
        result.code = SpecFunCode::OverflowErr;
        result.val = f64::INFINITY;
        result.err = f64::INFINITY;
        result
    }
}
