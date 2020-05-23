use super::airy_data::*;
use crate::result::{SpecFunCode, SpecFunResult};

fn zero_f(z: f64) -> f64 {
    let pre = z.powf(2.0 / 3.0);
    let zi2 = 1.0 / (z * z);
    let zi4 = zi2 * zi2;
    let t1 = 5.0 / 48.0 * zi2;
    let t2 = -5.0 / 36.0 * zi4;
    let t3 = 77125.0 / 82944.0 * zi4 * zi2;
    let t4 = -108056875.0 / 6967296.0 * zi4 * zi4;
    pre * (1.0 + t1 + t2 + t3 + t4)
}

fn zero_g(z: f64) -> f64 {
    let pre = z.powf(2.0 / 3.0);
    let zi2 = 1.0 / (z * z);
    let zi4 = zi2 * zi2;
    let t1 = -7.0 / 48.0 * zi2;
    let t2 = 35.0 / 288.0 * zi4;
    let t3 = -181223.0 / 207360.0 * zi4 * zi2;
    let t4 = 18683371.0 / 1244160.0 * zi4 * zi4;
    pre * (1.0 + t1 + t2 + t3 + t4)
}

/// Compute the nth zero of Ai(x)
pub(crate) fn airy_zero_ai_e(n: usize) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if n == 0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if n < (*AI_ZEROS).len() {
        result.val = (*AI_ZEROS)[n];
        result.err = f64::EPSILON * result.val.abs();
    } else {
        let z = 3.0 * std::f64::consts::PI / 8.0 * (4 * n - 1) as f64;
        let f = zero_f(z);
        result.val = -f;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    }

    result
}

/// Compute the nth zero of Bi(x)
pub(crate) fn airy_zero_bi_e(n: usize) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if n == 0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if n < (*BI_ZEROS).len() {
        result.val = (*BI_ZEROS)[n];
        result.err = f64::EPSILON * result.val.abs();
    } else {
        let z = 3.0 * std::f64::consts::PI / 8.0 * (4 * n - 3) as f64;
        let f = zero_f(z);
        result.val = -f;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    }

    result
}

/// Compute the nth zero of Ai'(x)
pub(crate) fn airy_zero_ai_deriv_e(n: usize) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if n == 0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if n < (*AIP_ZEROS).len() {
        result.val = (*AIP_ZEROS)[n];
        result.err = f64::EPSILON * result.val.abs();
    } else {
        let z = 3.0 * std::f64::consts::PI / 8.0 * (4 * n - 3) as f64;
        let g = zero_g(z);
        result.val = -g;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    }

    result
}

/// Compute the nth zero of Bi'(x)
pub(crate) fn airy_zero_bi_deriv_e(n: usize) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if n == 0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if n < (*BIP_ZEROS).len() {
        result.val = (*BIP_ZEROS)[n];
        result.err = f64::EPSILON * result.val.abs();
    } else {
        let z = 3.0 * std::f64::consts::PI / 8.0 * (4 * n - 1) as f64;
        let g = zero_g(z);
        result.val = -g;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
    }

    result
}
