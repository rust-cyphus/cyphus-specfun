// Author: Logan A. Morrison
// Original Author: Gerard Jungman (GSL)
//
// Any errors in this code are due to my own translation errors and not that
// of the original author G. Jungman.
//
// Copyrights from GSL:
//
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

//! Functionality for the Hurwitz and Riemann-zeta functions. These functions
//! are defined as:
//!        zeta(s) = sum_{n=1}^{\infty} 1/n^s
//!     hzeta(s,q) = sum_{n=0}^{\infty} 1 / (q + n)^s (term with q+n =0 excluded)
//!
//! # Examples
//! Hurwitz-zeta function
//! ```
//! let s = 4.0_f64;
//! let q = 3.0_f64;
//! let analytic = std::f64::consts::PI.powi(4) / 90.0 - 17.0 / 16.0;
//! assert!((s.hzeta(q) - analytic).abs() < 1e-10);
//! ```
//! Riemann-zeta function
//! ```
//! let s = 2.0_f64;
//! let analytic = std::f64::consts::PI.powi(2) / 6.0;
//! assert!((s.zeta() - analytic).abs() < 1e-10);
//! ```

use std::f64::consts::{LN_2};

use crate::consts::{EUL_GAMMA, ROOT5_DBL_EPS};
use crate::result::{SpecFunCode, SpecFunResult};

mod core;
mod data;

use crate::elementary::multiply_e;
use crate::exp::core::exp_e;
use crate::gamma::Gamma;
use crate::zeta::{core::*, data::*};

/// Implementation for the Riemann- and Hurwitz-zeta functions.
pub trait Zeta {
    /// Compute the Hurwitz-zeta function with error estimate.
    fn hzeta_e(&self, q: f64) -> SpecFunResult<f64>;
    /// Compute the Hurwitz-zeta function.
    fn hzeta(&self, q: f64) -> f64;
    /// Compute the Riemann-zeta function with error estimate.
    fn zeta_e(&self) -> SpecFunResult<f64>;
    /// Compute the Riemann-zeta function
    fn zeta(&self) -> f64;
    /// Compute zeta(self) - 1 with error estimate
    fn zeta_m1_e(&self) -> SpecFunResult<f64>;
    /// Compute zeta(self) - 1
    fn zeta_m1(&self) -> f64;
    /// Compute the Dirichlet eta function with error estimate
    fn eta_e(&self) -> SpecFunResult<f64>;
    /// Compute the Dirichlet eta function
    fn eta(&self) -> f64;
}

impl Zeta for f64 {
    fn hzeta_e(&self, q: f64) -> SpecFunResult<f64> {
        let s = *self;
        if s <= 1.0 || q <= 0.0 {
            let result = SpecFunResult {
                val: std::f64::NAN,
                err: std::f64::NAN,
                code: SpecFunCode::DomainErr,
            };
            result.issue_warning("hzeta_e", &[*self, q]);
            result
        } else {
            let max_bits = 54.0;
            let ln_term0 = -s * q.ln();

            if ln_term0 < std::f64::MIN_POSITIVE.ln() + 1.0 {
                let result = SpecFunResult {
                    val: 0.0,
                    err: 0.0,
                    code: SpecFunCode::UnderflowErr,
                };
                result.issue_warning("hzeta_e", &[*self, q]);
                result
            } else if ln_term0 > std::f64::MAX.ln() - 1.0 {
                let result = SpecFunResult {
                    val: 0.0,
                    err: 0.0,
                    code: SpecFunCode::OverflowErr,
                };
                result.issue_warning("hzeta_e", &[*self, q]);
                result
            } else if (s > max_bits && q < 1.0) || (s > 0.5 * max_bits && q < 0.25) {
                let val = q.powf(-s);
                let err = 2.0 * std::f64::EPSILON * val.abs();
                SpecFunResult {
                    val,
                    err,
                    code: SpecFunCode::Success,
                }
            } else if s > 0.5 * max_bits && q < 1.0 {
                let p1 = q.powf(-s);
                let p2 = (q / (1.0 + q)).powf(s);
                let p3 = (q / (2.0 + q)).powf(s);
                let val = p1 * (1.0 + p2 + p3);
                let err = std::f64::EPSILON * (0.5 * s + 2.0) * val.abs();
                SpecFunResult {
                    val,
                    err,
                    code: SpecFunCode::Success,
                }
            } else {
                // Euler-Maclaurin summation formula
                // [Moshier, p. 400, with several typo corrections]
                let jmax = 12;
                let kmax = 10;
                let pmax = (kmax as f64 + q).powf(-s);
                let mut scp = s;
                let mut pcp = pmax / (kmax as f64 + q);
                let mut ans = pmax * ((kmax as f64 + q) / (s - 1.0) + 0.5);

                ans += (0..kmax).fold(0.0, |acc, k| acc + (k as f64 + q).powf(-s));

                for j in 0..(jmax + 1) {
                    let delta = HZETA_C[j + 1] * scp * pcp;
                    ans += delta;
                    if (delta / ans).abs() < 0.5 * std::f64::EPSILON {
                        break;
                    }
                    scp *= (s + (2 * j) as f64 + 1.0) * (s + (2 * j) as f64 + 2.0);
                    pcp /= (kmax as f64 + q) * (kmax as f64 + q);
                }
                let err = 2.0 * (jmax as f64 + 1.0) * std::f64::EPSILON * ans.abs();
                SpecFunResult {
                    val: ans,
                    err,
                    code: SpecFunCode::Success,
                }
            }
        }
    }
    fn hzeta(&self, q: f64) -> f64 {
        self.hzeta_e(q).val
    }
    fn zeta_e(&self) -> SpecFunResult<Self> {
        let s = *self;
        let mut result = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };

        if (s - 1.0).abs() < std::f64::EPSILON {
            result.val = f64::NAN;
            result.err = f64::NAN;
            result.code = SpecFunCode::DomainErr;
            result.issue_warning("zeta_e", &[*self]);
            result
        } else if s >= 0.0 {
            riemann_zeta_sgt0(s)
        } else {
            // reflection formula, [Abramowitz+Stegun, 23.2.5]
            let zeta_one_minus_s = riemann_zeta1ms_slt0(s);
            let fmod = s - (s / 2.0).floor() * 2.0;
            let sin_term = if fmod.abs() < std::f64::EPSILON {
                0.0
            } else {
                let fmod = s - (s / 4.0).floor() * 4.0;
                (0.5 * std::f64::consts::PI * fmod).sin() / std::f64::consts::PI
            };

            if sin_term.abs() < std::f64::EPSILON {
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
                let n = (-s / 10.0).floor() as usize;
                let fs = s + 10.0 * n as f64;
                let p = (2.0 * std::f64::consts::PI).powf(fs) / twopi_pow[n];
                let g = (1.0 - s).gamma_e();

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
    fn zeta(&self) -> f64 {
        self.zeta_e().val
    }
    fn zeta_m1_e(&self) -> SpecFunResult<f64> {
        let s = *self;

        if s <= 5.0 {
            let mut result = s.zeta_e();
            result.val -= 1.0;
            result
        } else if s < 15.0 {
            riemann_zeta_minus_1_intermediate_s(s)
        } else {
            riemann_zeta_minus1_large_s(s)
        }
    }
    fn zeta_m1(&self) -> f64 {
        self.zeta_m1_e().val
    }
    fn eta_e(&self) -> SpecFunResult<f64> {
        let s = *self;
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
            let z = s.zeta_e();
            let p = exp_e((1.0 - s) * LN_2);
            let _m = multiply_e(1.0 - p.val, z.val);
            result.err = (p.err * (LN_2 * (1.0 - s)) * z.val).abs() + z.err * p.val.abs();
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            result
        }
    }
    fn eta(&self) -> f64 {
        (*self).eta_e().val
    }
}

impl Zeta for i32 {
    fn zeta_e(&self) -> SpecFunResult<f64> {
        let n = *self;
        let mut result = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };

        if n < 0 {
            if n % 2 == 0 {
                // Exactly zero at even negative integers
                result
            } else if n > -(ZETA_NEG_TABLE_NMAX as i32) {
                result.val = ZETA_NEG_INT_TABLE[(-(n + 1) / 2) as usize];
                result.err = 2.0 * f64::EPSILON * result.val.abs();
                result
            } else {
                (n as f64).zeta_e()
            }
        } else if n == 1 {
            result.val = f64::NAN;
            result.err = f64::NAN;
            result.code = SpecFunCode::DomainErr;
            result.issue_warning("zeta_e", &[n as f64]);
            result
        } else if n <= ZETA_POS_TABLE_NMAX as i32 {
            result.val = 1.0 + ZETAM1_POS_INT_TABLE[n as usize];
            result.err = 2.0 * f64::EPSILON * result.val.abs();
            result
        } else {
            result.val = 1.0;
            result.err = f64::EPSILON;
            result
        }
    }
    fn zeta(&self) -> f64 {
        self.zeta_e().val
    }
    fn hzeta_e(&self, q: f64) -> SpecFunResult<f64> {
        (*self as f64).hzeta_e(q)
    }
    fn hzeta(&self, q: f64) -> f64 {
        (*self as f64).hzeta(q)
    }
    fn zeta_m1_e(&self) -> SpecFunResult<f64> {
        let n = *self;
        let mut result = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };

        if n < 0 {
            if n % 2 == 0 {
                result
            } else if n > -(ZETA_NEG_TABLE_NMAX as i32) {
                result.val = ZETA_NEG_INT_TABLE[(-(n + 1) / 2) as usize] - 1.0;
                result.err = 2.0 * f64::EPSILON * result.val.abs();
                result
            } else {
                (n as f64).zeta_e()
            }
        } else if n == 1 {
            result.val = f64::NAN;
            result.err = f64::NAN;
            result.code = SpecFunCode::DomainErr;
            result.issue_warning("zeta_m1_e", &[n as f64]);
            result
        } else if n <= ZETA_POS_TABLE_NMAX as i32 {
            result.val = ZETAM1_POS_INT_TABLE[n as usize];
            result.err = 2.0 * f64::EPSILON * result.val.abs();
            result
        } else {
            (n as f64).zeta_m1_e()
        }
    }
    fn zeta_m1(&self) -> f64 {
        self.zeta_m1_e().val
    }
    fn eta_e(&self) -> SpecFunResult<f64> {
        let mut result = SpecFunResult {
            val: 0.0,
            err: 0.0,
            code: SpecFunCode::Success,
        };

        let n = *self;

        if n > ETA_POS_TABLE_NMAX as i32 {
            result.val = 1.0;
            result.err = f64::EPSILON;
            result
        } else if n >= 0 {
            result.val = ETA_POS_INT_TABLE[n as usize];
            result.err = 2.0 * f64::EPSILON * result.val.abs();
            result
        } else if n % 2 == 0 {
            result
        } else if n > -(ETA_NEG_TABLE_NMAX as i32) {
            result.val = ETA_NEG_INT_TABLE[(-(n + 1) / 2) as usize];
            result.err = 2.0 * f64::EPSILON * result.val.abs();
            result
        } else {
            let z = n.zeta_e();
            let p = exp_e((1.0 - n as f64) * LN_2);
            result = multiply_e(-p.val, z.val);
            result.err =
                (p.err * (LN_2 * (1.0 - n as f64)) * z.val).abs() + z.err * p.val.abs();
            result.err += 2.0 * f64::EPSILON * result.val.abs();
            result
        }
    }
    fn eta(&self) -> f64 {
        self.eta_e().val
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::{SQRT_DLB_EPS};
    use crate::test_utils::*;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131072.0 * f64::EPSILON;

    #[test]
    fn test_zeta_int_e() {
        test_sf_check_result_and_code(
            (-61).zeta_e(),
            -3.306_608_987_657_758e34,
            TOL0,
            SpecFunCode::Success,
        );

        test_sf_check_result_and_code((-8).zeta_e(), 0.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-6).zeta_e(), 0.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code(
            (-5).zeta_e(),
            -0.003_968_253_968_253_968,
            TOL0,
            SpecFunCode::Success,
        );

        test_sf_check_result_and_code((-4).zeta_e(), 0.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-3).zeta_e(), 1.0 / 120.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-2).zeta_e(), 0.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-1).zeta_e(), -1.0 / 12.0, TOL0, SpecFunCode::Success);

        test_sf_check_result_and_code(
            5.zeta_e(),
            1.036_927_755_143_37,
            TOL0,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            31.zeta_e(),
            1.000_000_000_465_662_8,
            TOL0,
            SpecFunCode::Success,
        );

        test_sf_check_result_and_code(
            (-151).zeta_e(),
            8.195_215_221_831_378e143,
            TOL2,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            (-51).zeta_e(),
            9.689_957_887_463_594e24,
            TOL1,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            (-5).zeta_e(),
            -0.003_968_253_968_253_968,
            TOL1,
            SpecFunCode::Success,
        );
    }
    #[test]
    fn test_hzeta() {
        let s = 4.0_f64;
        let q = 3.0_f64;
        let analytic = std::f64::consts::PI.powi(4) / 90.0 - 17.0 / 16.0;
        assert!((s.hzeta(q) - analytic).abs() < 1e-10);
    }
    #[test]
    fn test_zeta() {
        let s = 2.0_f64;
        let analytic = std::f64::consts::PI.powi(2) / 6.0;
        println!("{:?}", s.zeta());
        assert!((s.zeta() - analytic).abs() < 1e-10);
    }

    #[test]
    fn test_zeta_i32() {
        let s: i32 = 2;
        let analytic = std::f64::consts::PI.powi(2) / 6.0;
        println!("{:?}", s.zeta());
        assert!((s.zeta() - analytic).abs() < 1e-10);
    }
}
