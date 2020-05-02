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

use crate::result::SpecFunResult;

mod data;
pub(crate) mod eta;
pub(crate) mod hurwitz;
mod internal;
pub(crate) mod riemann;

use crate::zeta::{eta::*, hurwitz::*, riemann::*};

use num::Float;

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

macro_rules! impl_zeta_integer {
    ($T:ty) => {
        impl Zeta for $T {
            fn zeta_e(&self) -> SpecFunResult<f64> {
                zeta_int_e(*self as i64)
            }
            fn zeta(&self) -> f64 {
                zeta_int_e(*self as i64).val
            }
            fn hzeta_e(&self, q: f64) -> SpecFunResult<f64> {
                hzeta_e(*self as f64, q)
            }
            fn hzeta(&self, q: f64) -> f64 {
                hzeta_e(*self as f64, q).val
            }
            fn zeta_m1_e(&self) -> SpecFunResult<f64> {
                zeta_m1_int_e(*self as i64)
            }
            fn zeta_m1(&self) -> f64 {
                zeta_m1_int_e(*self as i64).val
            }
            fn eta_e(&self) -> SpecFunResult<f64> {
                eta_int_e(*self as i64)
            }
            fn eta(&self) -> f64 {
                eta_int_e(*self as i64).val
            }
        }
    };
}

macro_rules! impl_zeta_float {
    ($T:ty) => {
        impl Zeta for $T {
            fn zeta_e(&self) -> SpecFunResult<f64> {
                zeta_e(*self as f64)
            }
            fn zeta(&self) -> f64 {
                zeta_e(*self as f64).val
            }
            fn hzeta_e(&self, q: f64) -> SpecFunResult<f64> {
                hzeta_e(*self as f64, q)
            }
            fn hzeta(&self, q: f64) -> f64 {
                hzeta_e(*self as f64, q).val
            }
            fn zeta_m1_e(&self) -> SpecFunResult<f64> {
                zeta_m1_e(*self as f64)
            }
            fn zeta_m1(&self) -> f64 {
                zeta_m1_e(*self as f64).val
            }
            fn eta_e(&self) -> SpecFunResult<f64> {
                eta_e(*self as f64)
            }
            fn eta(&self) -> f64 {
                eta_e(*self as f64).val
            }
        }
    };
}

impl_zeta_integer!(u8);
impl_zeta_integer!(i8);
impl_zeta_integer!(i16);
impl_zeta_integer!(u16);
impl_zeta_integer!(i32);
impl_zeta_integer!(u32);
impl_zeta_integer!(i64);

impl_zeta_float!(f32);
impl_zeta_float!(f64);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::SQRT_DLB_EPS;
    use crate::result::SpecFunCode;
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

        test_sf_check_result_and_code(5.zeta_e(), 1.036_927_755_143_37, TOL0, SpecFunCode::Success);
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
    fn test_zetam1_int_e() {
        test_sf_check_result_and_code((-8).zeta_m1_e(), -1.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-6).zeta_m1_e(), -1.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-4).zeta_m1_e(), -1.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-3).zeta_m1_e(), -119.0 / 120.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-2).zeta_m1_e(), -1.0, TOL0, SpecFunCode::Success);
        test_sf_check_result_and_code((-1).zeta_m1_e(), -13.0 / 12.0, TOL0, SpecFunCode::Success);

        test_sf_check_result_and_code(
            5.zeta_m1_e(),
            0.0369277551433699263313655,
            TOL0,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            31.zeta_m1_e(),
            0.0000000004656629065033784,
            TOL0,
            SpecFunCode::Success,
        );
    }

    #[test]
    fn test_zeta_e() {
        test_sf_check_result_and_code(
            (-151.0).zeta_e(),
            8.195215221831378294e+143,
            TOL2,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            (-51.0).zeta_e(),
            9.68995788746359406565e+24,
            TOL1,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            (-5.0).zeta_e(),
            -0.003968253968253968253968,
            TOL1,
            SpecFunCode::Success,
        );

        test_sf_check_result_and_code((-8.0).zeta_e(), 0.0, TOL1, SpecFunCode::Success);
        test_sf_check_result_and_code((-6.0).zeta_e(), 0.0, TOL1, SpecFunCode::Success);
        test_sf_check_result_and_code((-4.0).zeta_e(), 0.0, TOL1, SpecFunCode::Success);
        test_sf_check_result_and_code((-3.0).zeta_e(), 1.0 / 120.0, TOL1, SpecFunCode::Success);
        test_sf_check_result_and_code((-2.0).zeta_e(), 0.0, TOL1, SpecFunCode::Success);
        test_sf_check_result_and_code((-1.0).zeta_e(), -1.0 / 12.0, TOL1, SpecFunCode::Success);

        test_sf_check_result_and_code(
            (-0.5).zeta_e(),
            -0.207886224977354566017307,
            TOL1,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            (-1e-10).zeta_e(),
            -0.49999999990810614668948,
            TOL1,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(0.0.zeta_e(), -0.5, TOL1, SpecFunCode::Success);
        test_sf_check_result_and_code(
            (1e-10).zeta_e(),
            -0.50000000009189385333058,
            TOL1,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            0.5.zeta_e(),
            -1.460354508809586812889499,
            TOL1,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            (1.0 - 1.0 / 1024.0).zeta_e(),
            -1023.4228554489429787,
            TOL1,
            SpecFunCode::Success,
        );
        test_sf_check_result_and_code(
            (1.0 + 1.0 / 1048576.0).zeta_e(),
            1.0485765772157343441e+06,
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
