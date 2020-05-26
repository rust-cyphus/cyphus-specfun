use crate::result::{SpecFunCode, SpecFunResult};

use crate::gamma::mono::gamma_e;
use crate::zeta::{data::*, internal::*};
use std::num::FpCategory;

pub(crate) fn zeta_e(s: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult::<f64>::default();

    if (s - 1.0).abs() < std::f64::EPSILON {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result.issue_warning("zeta_e", &[s]);
        result
    } else if s >= 0.0 {
        riemann_zeta_sgt0(s)
    } else {
        // reflection formula, [Abramowitz+Stegun, 23.2.5]
        let zeta_one_minus_s = riemann_zeta1ms_slt0(s);
        let fmod = s - (s / 2.0).trunc() * 2.0;
        let sin_term = if fmod.classify() == FpCategory::Zero {
            0.0
        } else {
            let fmod = s - (s / 4.0).trunc() * 4.0;
            (0.5 * std::f64::consts::PI * fmod).sin() / std::f64::consts::PI
        };

        if sin_term.classify() == FpCategory::Zero {
            result
        } else if s > -170.0 {
            // We have to be careful about losing digits
            // in calculating pow(2 Pi, s). The gamma
            // function is fine because we were careful
            // with that implementation.
            // We keep an array of (2 Pi)^(10 n).
            let twopi_pow: [f64; 18] = [
                1.0,
                9.589_560_061_550_902e7_f64,
                9.195_966_217_409_212e15_f64,
                8.818_527_036_583_87e23_f64,
                8.456_579_467_173_15e31_f64,
                8.109_487_671_573_504e39_f64,
                7.776_641_909_496_07e47_f64,
                7.457_457_466_828_644e55_f64,
                7.151_373_628_461_452e63_f64,
                6.857_852_693_272_23e71_f64,
                6.576_379_029_540_266e79_f64,
                6.306_458_169_130_021e87_f64,
                6.047_615_938_853_066e95_f64,
                5.799_397_627_482_402e103_f64,
                5.561_367_186_955_83e111_f64,
                5.333_106_466_365_131e119_f64,
                5.114_214_477_385_392e127_f64,
                4.904_306_689_854_036_5e135_f64,
            ];
            let n = (-s / 10.0).floor();
            let fs = s + 10.0 * n;
            let p = (2.0 * std::f64::consts::PI).powf(fs) / twopi_pow[n as usize];
            let g = gamma_e(1.0 - s);

            result.val = p * g.val * sin_term * zeta_one_minus_s.val;
            result.err = (p * g.val * sin_term).abs() * zeta_one_minus_s.err;
            result.err += (p * sin_term * zeta_one_minus_s.val).abs() * g.err;
            result.err += f64::EPSILON * (s.abs() + 2.0) * result.val.abs();
            result
        } else {
            //  The actual zeta function may or may not
            // overflow here. But we have no easy way
            // to calculate it when the prefactor(s)
            // overflow. Trying to use log's and exp
            // is no good because we loose a couple
            // digits to the exp error amplification.
            // When we gather a little more patience,
            // we can implement something here. Until
            // then just give up.
            result.val = f64::NAN;
            result.err = f64::NAN;
            result.code = SpecFunCode::OverflowErr;
            result.issue_warning("zeta_e", &[s]);
            result
        }
    }
}

pub(crate) fn zeta_m1_e(s: f64) -> SpecFunResult<f64> {
    if s <= 5.0 {
        let mut result = zeta_e(s);
        result.val -= 1.0;
        result
    } else if s < 15.0 {
        riemann_zeta_minus_1_intermediate_s(s)
    } else {
        riemann_zeta_minus1_large_s(s)
    }
}

pub(crate) fn zeta_int_e(n: i64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if n < 0 {
        if n % 2 == 0 {
            // Exactly zero at even negative integers
            result
        } else if n > -(ZETA_NEG_TABLE_NMAX as i64) {
            result.val = ZETA_NEG_INT_TABLE[(-(n + 1) / 2) as usize];
            result.err = 2.0 * f64::EPSILON * result.val.abs();
            result
        } else {
            zeta_e(n as f64)
        }
    } else if n == 1 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result.issue_warning("zeta_e", &[n as f64]);
        result
    } else if n <= ZETA_POS_TABLE_NMAX as i64 {
        result.val = 1.0 + ZETAM1_POS_INT_TABLE[n as usize];
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        result.val = 1.0;
        result.err = f64::EPSILON;
        result
    }
}

pub(crate) fn zeta_m1_int_e(n: i64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if n < 0 {
        if n % 2 == 0 {
            result.val = -1.0;
            result.err = 0.0;
            result
        } else if n > -(ZETA_NEG_TABLE_NMAX as i64) {
            result.val = ZETA_NEG_INT_TABLE[(-(n + 1) / 2) as usize] - 1.0;
            result.err = 2.0 * f64::EPSILON * result.val.abs();
            result
        } else {
            zeta_e(n as f64)
        }
    } else if n == 1 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result.issue_warning("zeta_m1_e", &[n as f64]);
        result
    } else if n <= ZETA_POS_TABLE_NMAX as i64 {
        result.val = ZETAM1_POS_INT_TABLE[n as usize];
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        zeta_m1_e(n as f64)
    }
}
