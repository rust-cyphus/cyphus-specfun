use crate::consts::{EUL_GAMMA, ROOT5_DBL_EPS};
use crate::result::{SpecFunCode, SpecFunResult};

use crate::elementary::multiply_e;
use crate::exp::core::exp_e;

use super::data::*;
use super::riemann::{zeta_e, zeta_int_e};

use std::f64::consts::LN_2;

pub(crate) fn eta_e(s: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if s > 100.0 {
        result.val = 1.0;
        result.err = f64::EPSILON;
        result
    } else if (s - 1.0).abs() < 10.0 * ROOT5_DBL_EPS {
        let del = s - 1.0;
        let c0 = LN_2;
        let c1 = LN_2 * (EUL_GAMMA - 0.5 * LN_2);
        let c2 = -0.032_686_296_279_449_3_f64;
        let c3 = 0.001_568_991_705_415_515_f64;
        let c4 = 0.000_749_872_421_120_475_4_f64;
        result.val = c0 + del * (c1 + del * (c2 + del * (c3 + del * c4)));
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        let z = zeta_e(s);
        let p = exp_e((1.0 - s) * LN_2);
        let _m = multiply_e(1.0 - p.val, z.val);
        result.err = (p.err * (LN_2 * (1.0 - s)) * z.val).abs() + z.err * p.val.abs();
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

pub(crate) fn eta_int_e(n: i64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if n > ETA_POS_TABLE_NMAX as i64 {
        result.val = 1.0;
        result.err = f64::EPSILON;
        result
    } else if n >= 0 {
        result.val = ETA_POS_INT_TABLE[n as usize];
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else if n % 2 == 0 {
        result
    } else if n > -(ETA_NEG_TABLE_NMAX as i64) {
        result.val = ETA_NEG_INT_TABLE[(-(n + 1) / 2) as usize];
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        let z = zeta_int_e(n);
        let p = exp_e((1.0 - n as f64) * LN_2);
        result = multiply_e(-p.val, z.val);
        result.err = (p.err * (LN_2 * (1.0 - n as f64)) * z.val).abs() + z.err * p.val.abs();
        result.err += 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}
