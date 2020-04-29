// Copyright from GSL
/* specfunc/trig.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

//! Implementations of various trig functions and their assoiciated errors.

use crate::cheb::cheb_eval_e;
use crate::consts::{LN_DBL_EPS, LN_DBL_MAX, ROOT4_DBL_EPS, SQRT_DLB_EPS};
use crate::logarithm::{complex_ln_e, ln_p1_e};
use crate::result::{SpecFunCode, SpecFunResult};
use num::Complex;
mod data;
use data::*;

pub trait Trig {
    /// Compute sine of a number along with the error estimate.
    ///
    /// ## Examples
    /// ```
    /// let x = 1.0;
    ///
    /// assert!((x.sin_e().val - 0.841470984807897).abs() < 1e-10)
    /// ```
    fn sin_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute cosine of a number along with the error estimate.
    ///
    /// ## Examples
    /// ```
    /// let x = 1.0;
    ///
    /// assert!((x.cos_e().val - 0.540302305868140).abs() < 1e-10)
    /// ```
    fn cos_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute sqrt(x^2+y^2) along with the error estimate.
    ///
    /// ## Examples
    /// ```
    /// let x = 1.0;
    /// let y = 2.0;
    ///
    /// assert!((x.hypot_e(y).val - 2.2360679774997897).abs() < 1e-10)
    /// ```
    fn hypot_e(&self, other: Self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute sinc(x) = sin(pi*x) / (pi*x) along with error estimate.
    ///
    /// # Examples
    /// ```
    /// let x = 1.5_f64
    /// assert!((x.sinc_e().val - -0.21220659078919378).abs() < 1e-10)
    /// ```
    fn sinc_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute sinc(x) = sin(pi*x) / (pi*x).
    ///
    /// # Examples
    /// ```
    /// let x = 1.5_f64
    /// assert!((x.sinc() - -0.21220659078919378).abs() < 1e-10)
    /// ```
    fn sinc(&self) -> Self
    where
        Self: Sized;
    /// Compute log(sinh(x)) for a positive real number x and the associated error
    ///
    /// # Examples
    /// ```
    /// let x = 10.0_f64;
    ///
    /// assert!((x.lnsinh_e().val - 9.3068528173789011).abs() < 1e-10);
    /// ```
    fn lnsinh_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute log(sinh(x)) for a positive real number x.
    ///
    /// # Examples
    /// ```
    /// let x = 10.0_f64;
    ///
    /// assert!((x.lnsinh() - 9.3068528173789011).abs() < 1e-10);
    /// ```
    fn lnsinh(&self) -> Self
    where
        Self: Sized;
    /// Compute ln(cosh(x)) for a real number x and the associated error
    ///
    /// # Examples
    /// ```
    /// let x = 10.0_f64;
    ///
    /// assert!((x.lncosh_e().val -9.3068528215012083).abs() < 1e-10);
    /// ```
    fn lncosh_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute ln(cosh(x)) for a real number x.
    ///
    /// # Examples
    /// ```
    /// let x = 10.0_f64;
    ///
    /// assert!((x.lncosh() -9.3068528215012083).abs() < 1e-10);
    /// ```
    fn lncosh(&self) -> Self
    where
        Self: Sized;
}

impl Trig for f64 {
    fn sin_e(&self) -> SpecFunResult<f64> {
        sin_e(*self)
    }
    fn cos_e(&self) -> SpecFunResult<f64> {
        cos_e(*self)
    }
    fn hypot_e(&self, other: f64) -> SpecFunResult<f64> {
        hypot_e(*self, other)
    }
    fn sinc_e(&self) -> SpecFunResult<f64> {
        sin_e(*self)
    }
    fn sinc(&self) -> f64 {
        sin_e(*self).val
    }
    fn lnsinh_e(&self) -> SpecFunResult<f64> {
        lnsinh_e(*self)
    }
    fn lnsinh(&self) -> f64 {
        lnsinh_e(*self).val
    }
    fn lncosh_e(&self) -> SpecFunResult<f64> {
        lncosh_e(*self)
    }
    fn lncosh(&self) -> f64 {
        lncosh_e(*self).val
    }
}

/// Compute sine of a number along with the error estimate.
///
/// ## Examples
/// ```
/// let x = 1.0;
///
/// assert!((sin_e(x).val - 0.841470984807897).abs() < 1e-10)
/// ```
pub fn sin_e(x: f64) -> SpecFunResult<f64> {
    let p1 = 7.853_981_256_484_985e-1;
    let p2 = 3.774_894_707_930_798e-8;
    let p3 = 2.695_151_429_079_059_5e-15;

    let sgn_x = x.signum();
    let abs_x = x.abs();

    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if abs_x < ROOT4_DBL_EPS {
        let x2 = x * x;
        result.val = x * (1.0 - x2 / 6.0);
        result.err = (x * x2 * x2 / 100.0).abs();
        result
    } else {
        let mut sgn_result = sgn_x;
        let mut y = (abs_x / (0.25 * std::f64::consts::PI)).floor();
        let mut octant = (y - (y * (-3_f64).exp2()).floor() * 3_f64.exp2()) as i32;

        if octant % 2 != 0 {
            octant += 1;
            octant &= 0o7;
            y += 1.0;
        }
        if octant > 3 {
            octant -= 4;
            sgn_result *= -1.0;
        }

        let z = ((abs_x - y * p1) - y * p2) - y * p3;

        if octant == 0 {
            let t = 8.0 * z.abs() / std::f64::consts::PI - 1.0;
            result = cheb_eval_e(t, &SIN_DATA, -1.0, 1.0);
            result.val = z * (1.0 + z * z * result.val);
        } else {
            let t = 8.0 * z.abs() / std::f64::consts::PI - 1.0;
            result = cheb_eval_e(t, &COS_DATA, -1.0, 1.0);
            result.val = 1.0 - 0.5 * z * z * (1.0 - z * z * result.val);
        }

        result.val *= sgn_result;

        if abs_x > 1.0 / f64::EPSILON {
            result.err = result.val.abs();
        } else if abs_x > 100.0 / SQRT_DLB_EPS {
            result.err = 2.0 * abs_x * f64::EPSILON * result.val.abs();
        } else if abs_x > 0.1 / SQRT_DLB_EPS {
            result.err = 2.0 * SQRT_DLB_EPS * result.val.abs();
        } else {
            result.err = 2.0 * f64::EPSILON * result.val.abs();
        }
        result
    }
}

/// Compute cosine of a number along with the error estimate.
///
/// ## Examples
/// ```
/// let x = 1.0;
///
/// assert!((cos_e(x).val - 0.540302305868140).abs() < 1e-10)
/// ```
pub fn cos_e(x: f64) -> SpecFunResult<f64> {
    let p1 = 7.853_981_256_484_985e-1;
    let p2 = 3.774_894_707_930_798e-8;
    let p3 = 2.695_151_429_079_059_5e-15;

    let abs_x = x.abs();

    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if abs_x < ROOT4_DBL_EPS {
        let x2 = x * x;
        result.val = 1.0 - 0.5 * x2;
        result.err = (x2 * x2 / 12.0).abs();
        result
    } else {
        let mut sgn_result = 1.0;
        let mut y = (abs_x / (0.25 * std::f64::consts::PI)).floor();
        let mut octant = (y - (y * (-3_f64).exp2()).floor() * 3_f64.exp2()) as i32;

        if octant % 2 != 0 {
            octant += 1;
            octant &= 0o7;
            y += 1.0;
        }
        if octant > 3 {
            octant -= 4;
            sgn_result *= -1.0;
        }
        if octant > 1 {
            sgn_result *= -1.0;
        }

        let z = ((abs_x - y * p1) - y * p2) - y * p3;

        if octant == 0 {
            let t = 8.0 * z.abs() / std::f64::consts::PI - 1.0;
            result = cheb_eval_e(t, &COS_DATA, -1.0, 1.0);
            result.val = 1.0 - 0.5 * z * z * (1.0 - z * z * result.val);
        } else {
            let t = 8.0 * z.abs() / std::f64::consts::PI - 1.0;
            result = cheb_eval_e(t, &SIN_DATA, -1.0, 1.0);
            result.val = z * (1.0 + z * z * result.val);
        }

        result.val *= sgn_result;

        if abs_x > 1.0 / f64::EPSILON {
            result.err = result.val.abs();
        } else if abs_x > 100.0 / SQRT_DLB_EPS {
            result.err = 2.0 * abs_x * f64::EPSILON * result.val.abs();
        } else if abs_x > 0.1 / SQRT_DLB_EPS {
            result.err = 2.0 * SQRT_DLB_EPS * result.val.abs();
        } else {
            result.err = 2.0 * f64::EPSILON * result.val.abs();
        }
        result
    }
}

/// Compute sqrt(x^2+y^2) along with the error estimate.
///
/// ## Examples
/// ```
/// let x = 1.0;
/// let y = 2.0;
///
/// assert!((hypot_e(x, y).val - 2.2360679774997897).abs() < 1e-10)
/// ```
pub fn hypot_e(x: f64, y: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };
    if x.abs() < f64::EPSILON && y.abs() < f64::EPSILON {
        result
    } else {
        let a = x.abs();
        let b = y.abs();
        let min = a.min(b);
        let max = a.max(b);
        let rat = min / max;
        let root_term = (1.0 + rat * rat).sqrt();

        if max < f64::MAX / root_term {
            result.val = max * root_term;
            result.err = 2.0 * f64::EPSILON * result.val.abs();
            result
        } else {
            result.val = f64::INFINITY;
            result.err = f64::INFINITY;
            result.code = SpecFunCode::OverflowErr;
            result.issue_warning("hypot_e", &[x, y]);
            result
        }
    }
}

/// Compute sine of a complex number along with the error estimate.
///
/// ## Examples
/// ```
/// use num::Complex;
///
/// let z = Complex::new(1.0, 1.0);
/// let w = complex_sin_e(z).val;
///
/// assert!((w.re - 1.2984575814159773).abs() < 1e-10);
/// assert!((w.im - 0.6349639147847361).abs() < 1e-10);
/// ```
pub fn complex_sin_e(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut result = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
        code: SpecFunCode::Success,
    };
    if z.im.abs() < 1.0 {
        let sh = sinh_series(z.im);
        let chm1 = cosh_m1_series(z.im);
        result.val.re = z.re.sin() * (chm1 + 1.0);
        result.val.im = z.re.cos() * sh;
        result.err.re = 2.0 * f64::EPSILON * result.val.re.abs();
        result.err.im = 2.0 * f64::EPSILON * result.val.im.abs();
    } else if z.im.abs() < LN_DBL_MAX {
        let ex = z.im.exp();
        let ch = 0.5 * (ex + 1.0 / ex);
        let sh = 0.5 * (ex - 1.0 / ex);
        result.val.re = z.re.sin() * ch;
        result.val.im = z.re.cos() * sh;
        result.err.re = 2.0 * f64::EPSILON * result.val.re.abs();
        result.err.im = 2.0 * f64::EPSILON * result.val.im.abs();
    } else {
        result.val.re = f64::INFINITY;
        result.val.im = f64::INFINITY;
        result.err.re = f64::INFINITY;
        result.err.im = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("complex_sin_e", &[z]);
    }

    result
}

/// Compute cosine of a complex number along with the error estimate.
///
/// ## Examples
/// ```
/// use num::Complex;
///
/// let z = Complex::new(1.0, 1.0);
/// let w = complex_cos_e(z).val;
///
/// assert!((w.re - 0.83373002513114905).abs() < 1e-10);
/// assert!((w.im + 0.98889770576286510).abs() < 1e-10);
/// ```
pub fn complex_cos_e(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut result = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
        code: SpecFunCode::Success,
    };
    if z.im.abs() < 1.0 {
        let sh = sinh_series(z.im);
        let chm1 = cosh_m1_series(z.im);
        result.val.re = z.re.cos() * (chm1 + 1.0);
        result.val.im = -z.re.sin() * sh;
        result.err.re = 2.0 * f64::EPSILON * result.val.re.abs();
        result.err.im = 2.0 * f64::EPSILON * result.val.im.abs();
    } else if z.im.abs() < LN_DBL_MAX {
        let ex = z.im.exp();
        let ch = 0.5 * (ex + 1.0 / ex);
        let sh = 0.5 * (ex - 1.0 / ex);
        result.val.re = z.re.cos() * ch;
        result.val.im = -z.re.sin() * sh;
        result.err.re = 2.0 * f64::EPSILON * result.val.re.abs();
        result.err.im = 2.0 * f64::EPSILON * result.val.im.abs();
    } else {
        result.val.re = f64::INFINITY;
        result.val.im = f64::INFINITY;
        result.err.re = f64::INFINITY;
        result.err.im = f64::INFINITY;
        result.code = SpecFunCode::OverflowErr;
        result.issue_warning("complex_cos_e", &[z]);
    }

    result
}

/// Compute log(sin(z)) for a complex number z
///
/// ## Examples
/// ```
/// use num::Complex;
///
/// let z = Complex::new(1.0, 1.0);
/// let w = complex_lnsin_e(z).val;
///
/// assert!((w.re - 0.36838373142492511).abs() < 1e-10);
/// assert!((w.im - 0.45482023330994990).abs() < 1e-10);
/// ```
pub fn complex_lnsin_e(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut result = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
        code: SpecFunCode::Success,
    };

    if z.im > 60.0 {
        result.val.re = -std::f64::consts::LN_2 + z.im;
        result.val.im = 0.5 * std::f64::consts::PI - z.re;
        result.err.re = 2.0 * std::f64::EPSILON * result.val.re.abs();
        result.err.im = 2.0 * std::f64::EPSILON * result.val.im.abs();
    } else if z.im < -60.0 {
        result.val.re = -std::f64::consts::LN_2 - z.im;
        result.val.im = -0.5 * std::f64::consts::PI + z.re;
        result.err.re = 2.0 * std::f64::EPSILON * result.val.re.abs();
        result.err.im = 2.0 * std::f64::EPSILON * result.val.im.abs();
    } else {
        let sin = complex_sin_e(z);
        result = complex_ln_e(sin.val);
        if result.code == SpecFunCode::DomainErr {
            result.issue_warning("complex_lnsin_e", &[z]);
            return result;
        }
    }
    result.val.im = angle_restrict_symm(result.val.im);
    result
}

/// Compute log(sinh(x)) for a positive real number x as the associated error
///
/// # Examples
/// ```
/// let x = 10.0_f64;
///
/// assert!((lnsinh_e(x).val - 9.3068528173789011).abs() < 1e-10);
/// ```
pub fn lnsinh_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if x <= 0.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result.issue_warning("lnsinh_e", &[x]);
        result
    } else if x.abs() < 1.0 {
        result.val = sinh_series(x).ln();
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else if x < -0.5 * LN_DBL_EPS {
        result.val = x + (0.5 * (1.0 - (-2.0 * x).exp())).ln();
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        result.val = -std::f64::consts::LN_2 + x;
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

/// Compute log(cosh(x)) for a real number x as the associated error
///
/// # Examples
/// ```
/// let x = 10.0_f64;
///
/// assert!((lncosh_e(x).val -9.3068528215012083).abs() < 1e-10);
/// ```
pub fn lncosh_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };
    if x.abs() < 1.0 {
        ln_p1_e(cosh_m1_series(x))
    } else if x < -0.5 * LN_DBL_EPS {
        result.val = x.abs() + (0.5 * (1.0 + (-2.0 * x.abs()).exp())).ln();
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    } else {
        result.val = -std::f64::consts::LN_2 + x.abs();
        result.err = 2.0 * f64::EPSILON * result.val.abs();
        result
    }
}

/// Convert a complex number in polar form: z=r*e^(i*theta) into retangular
/// for, z=x+iy.
///
/// # Examples
/// ```
/// let r = 10.0_f64;
/// let theta = std::f64::consts::PI;
///
/// let z = polar_to_rect(r, theta);
///
/// assert!((z.val.re - r * theta.cos()).abs() < 1e-10);
/// assert!((z.val.im - r * theta.sin()).abs() < 1e-10);
/// ```
pub fn polar_to_rect(r: f64, theta: f64) -> SpecFunResult<Complex<f64>> {
    let mut result = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
        code: SpecFunCode::Success,
    };

    let t = angle_restrict_symm(theta);

    let c = t.cos();
    let s = t.sin();

    result.val.re = r * c;
    result.val.im = r * s;

    result.err.re = r * (s * t * f64::EPSILON).abs();
    result.err.re += 2.0 * f64::EPSILON * result.val.re;

    result.err.im = r * (c * t * f64::EPSILON).abs();
    result.err.im += 2.0 * f64::EPSILON * result.val.im;

    result
}

/// Convert a complex number in retangular form: z=x+iy into polar
/// from z=r*e^(i*theta).
///
/// # Examples
/// ```
/// use num::Complex;
///
/// let z = Complex::new(10.0, -4.0);
///
/// let (r,t) = rect_to_polar(z);
///
/// assert!((r.val - z.re.hypot(z.im)).abs() < 1e-10);
/// assert!((t.val - z.arg()).abs() < 1e-10);
/// ```
pub fn rect_to_polar(z: Complex<f64>) -> (SpecFunResult<f64>, SpecFunResult<f64>) {
    let mut theta = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };
    let r = hypot_e(z.re, z.im);
    if r.val > 0.0 {
        theta.val = z.im.atan2(z.re);
        theta.err = 2.0 * f64::EPSILON * theta.val.abs();
    } else {
        theta.val = f64::NAN;
        theta.err = f64::NAN;
        theta.code = SpecFunCode::DomainErr;
    }
    (r, theta)
}

/// Return a new angle in range (-pi,pi] and the associated error of the
/// computation.
///
/// # Examples
/// ```
/// let theta = 2.5 * std::f64::consts::PI;
/// let clipped = angle_restrict_symm_e(theta);
/// assert!((clipped.val - 0.5 * std::f64::consts::PI).abs() < 1e-10);
/// ```
pub fn angle_restrict_symm_e(theta: f64) -> SpecFunResult<f64> {
    // synthetic extended precision constants
    let p1 = 4.0 * 7.853_981_256_484_985e-1;
    let p2 = 4.0 * 3.774_894_707_930_798e-8;
    let p3 = 4.0 * 2.695_151_429_079_059_5e-15;
    let twopi = 2.0 * (p1 + p2 + p3);

    let y = theta.signum() * 2.0 * (theta.abs() / twopi).floor();
    let mut r = ((theta - y * p1) - y * p2) - y * p3;

    if r > std::f64::consts::PI {
        r = ((r - 2.0 * p1) - 2.0 * p2) - 2.0 * p3;
    } else if r < -std::f64::consts::PI {
        r = ((r + 2.0 * p1) + 2.0 * p2) + 2.0 * p3;
    }

    let val = r;

    if theta.abs() > 0.0625 / f64::EPSILON {
        let result = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::AccLossErr,
        };
        result.issue_warning("angle_restrict_symm_err_e", &[theta]);
        result
    } else if theta.abs() > 0.0625 / SQRT_DLB_EPS {
        let err = 2.0 * f64::EPSILON * (val - theta).abs();
        let code = SpecFunCode::Success;
        SpecFunResult { val, err, code }
    } else {
        let delta = (val - theta).abs();
        let err = 2.0
            * f64::EPSILON
            * if delta < std::f64::consts::PI {
                delta
            } else {
                std::f64::consts::PI
            };
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

/// Return a new angle in range (-pi,pi].
///
/// # Examples
/// ```
/// let theta = 2.5 * std::f64::consts::PI;
/// let clipped = angle_restrict_symm(theta);
/// assert!((clipped - 0.5 * std::f64::consts::PI).abs() < 1e-10);
/// ```
pub fn angle_restrict_symm(theta: f64) -> f64 {
    angle_restrict_symm_e(theta).val
}

/// Return a new angle in range [0,2pi) and the associated error of the
/// computation.
///
/// # Examples
/// ```
/// let theta = 2.5 * std::f64::consts::PI;
/// let clipped = angle_restrict_pos_e(theta);
/// assert!((clipped.val - 0.5 * std::f64::consts::PI).abs() < 1e-10);
/// ```
pub fn angle_restrict_pos_e(theta: f64) -> SpecFunResult<f64> {
    // synthetic extended precision constants
    let p1 = 4.0 * 7.853_981_256_484_985e-1;
    let p2 = 4.0 * 3.774_894_707_930_798e-8;
    let p3 = 4.0 * 2.695_151_429_079_059_5e-15;
    let twopi = 2.0 * (p1 + p2 + p3);

    let y = 2.0 * (theta / twopi).floor();

    let mut r = ((theta - y * p1) - y * p2) - y * p3;

    if r > twopi {
        r = ((r - 2.0 * p1) - 2.0 * p2) - 2.0 * p3;
    } else if r < 0.0 {
        // may happend due to fp rounding
        r = ((r + 2.0 * p1) + 2.0 * p2) + 2.0 * p3;
    }

    let mut val = r;

    if theta.abs() > 0.0625 / f64::EPSILON {
        val = f64::NAN;
        let err = val.abs();
        let result = SpecFunResult {
            val,
            err,
            code: SpecFunCode::AccLossErr,
        };
        result.issue_warning("angle_restrict_symm_err_e", &[theta]);
        result
    } else if theta.abs() > 0.0625 / SQRT_DLB_EPS {
        let err = f64::EPSILON * (val - theta).abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        let delta = (val - theta).abs();
        let err = 2.0
            * f64::EPSILON
            * if delta < std::f64::consts::PI {
                delta
            } else {
                std::f64::consts::PI
            };
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    }
}

/// Return a new angle in range [0,2pi).
///
/// # Examples
/// ```
/// let theta = 2.5 * std::f64::consts::PI;
/// let clipped = angle_restrict_pos(theta);
/// assert!((clipped - 0.5 * std::f64::consts::PI).abs() < 1e-10);
/// ```
pub fn angle_restrict_pos(theta: f64) -> f64 {
    angle_restrict_pos_e(theta).val
}

/// Compute the sin(x) along with the error given an error in x of dx.
///
/// # Examples
/// ```
/// let x = 1.0_f64;
/// assert!((sin_err_e(x).val - 0.84147098480789651).abs()<1e-10);
/// ```
pub fn sin_err_e(x: f64, dx: f64) -> SpecFunResult<f64> {
    let mut result = sin_e(x);
    result.err += (x.cos() * dx).abs();
    result.err += f64::EPSILON * result.val.abs();
    result
}

/// Compute the cos(x) along with the error given an error in x of dx.
///
/// # Examples
/// ```
/// let x = 1.0_f64;
/// assert!((cos_err_e(x).val - 0.54030230586813972).abs()<1e-10);
/// ```
pub fn cos_err_e(x: f64, dx: f64) -> SpecFunResult<f64> {
    let mut result = cos_e(x);
    result.err += (x.sin() * dx).abs();
    result.err += f64::EPSILON * result.val.abs();
    result
}

/// Compute sinc(x) = sin(pi*x) / (pi*x)
///
/// # Examples
/// ```
/// let x = 1.5_f64
/// assert!((sinc_e(x).val - -0.21220659078919378).abs() < 1e-10)
/// ```
pub fn sinc_e(x: f64) -> SpecFunResult<f64> {
    let ax = x.abs();
    if ax < 0.8 {
        // Do not go to the limit of the fit since
        // there is a zero there and the Chebyshev
        // accuracy will go to zero.
        cheb_eval_e(2.0 * ax - 1.0, &SINC_DATA, -1.0, 1.0)
    } else if ax < 100.0 {
        // Small arguments are no problem.
        // We trust the library sin() to
        // roughly machine precision.
        let val = (std::f64::consts::PI * ax).sin() / (ax * std::f64::consts::PI);
        let err = 2.0 * f64::EPSILON * val.abs();
        SpecFunResult {
            val,
            err,
            code: SpecFunCode::Success,
        }
    } else {
        // Large arguments must be handled seperately
        let r = std::f64::consts::PI * ax;
        let mut result = sin_e(r);
        result.val = result.val / r;
        result.err = result.err / r + 2.0 * f64::EPSILON * result.val;
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sin_e() {
        let x = 1.0;
        assert!((sin_e(x).val - 0.841470984807897).abs() < 1e-10);
        assert!((x.sin_e().val - 0.841470984807897).abs() < 1e-10);
    }
    #[test]
    fn test_cos_e() {
        let x = 1.0;
        assert!((cos_e(x).val - 0.54030230586813972).abs() < 1e-10);
        assert!((x.cos_e().val - 0.54030230586813972).abs() < 1e-10);
    }
    #[test]
    fn test_hypot_e() {
        let x = 1.0;
        let y = 2.0;
        assert!((hypot_e(x, y).val - 2.2360679774997897).abs() < 1e-10);
        assert!((x.hypot_e(y).val - 2.2360679774997897).abs() < 1e-10);
    }
    #[test]
    fn test_complex_sin_e() {
        let z1 = Complex::new(1.0, 1.0);
        let z2 = Complex::new(10.0, -4.0);
        let z3 = Complex::new(-20.0, 0.1);
        let z4 = Complex::new(-0.1, -2.0);

        let w1 = complex_sin_e(z1).val;
        let w2 = complex_sin_e(z2).val;
        let w3 = complex_sin_e(z3).val;
        let w4 = complex_sin_e(z4).val;

        assert!((w1.re - 1.2984575814159773).abs() < 1e-10);
        assert!((w1.im - 0.6349639147847361).abs() < 1e-10);

        assert!((w2.re - -14.856255163875255).abs() < 1e-10);
        assert!((w2.im - 22.898192550963759).abs() < 1e-10);

        assert!((w3.re - -0.91751378218801646).abs() < 1e-10);
        assert!((w3.im - 0.04087625387324457).abs() < 1e-10);

        assert!((w4.re - -0.3755928499348538).abs() < 1e-10);
        assert!((w4.im - -3.6087412126897430).abs() < 1e-10);
    }
    #[test]
    fn test_complex_cos_e() {
        let z = Complex::new(1.0, 1.0);
        let w = complex_cos_e(z).val;
        assert!((w.re - 0.83373002513114905).abs() < 1e-10);
        assert!((w.im + 0.98889770576286510).abs() < 1e-10);
    }
    #[test]
    fn test_complex_lnsin_e() {
        let z = Complex::new(2.0, -1.0);
        let w = complex_lnsin_e(z).val;
        // 0.39602537002984730+0.33538186897032753 I
        assert!((w.re - 0.39602537002984730).abs() < 1e-10);
        assert!((w.im - 0.33538186897032753).abs() < 1e-10);
    }
    #[test]
    fn test_lnsinh_e() {
        let x = 10.0_f64;
        assert!((lnsinh_e(x).val - 9.3068528173789011).abs() < 1e-10);
        assert!((x.lnsinh_e().val - 9.3068528173789011).abs() < 1e-10);
        assert!((x.lnsinh() - 9.3068528173789011).abs() < 1e-10);
    }
    #[test]
    fn test_lncosh_e() {
        let x = 10.0_f64;
        assert!((lncosh_e(x).val - 9.3068528215012083).abs() < 1e-10);
        assert!((x.lncosh_e().val - 9.3068528215012083).abs() < 1e-10);
        assert!((x.lncosh() - 9.3068528215012083).abs() < 1e-10);
    }
    #[test]
    fn test_angle_restrict_symm_e() {
        let theta = 2.5 * std::f64::consts::PI;
        let clipped = angle_restrict_symm_e(theta);
        assert!((clipped.val - 0.5 * std::f64::consts::PI).abs() < 1e-10);
    }
    #[test]
    fn test_angle_restrict_symm() {
        let theta = 2.5 * std::f64::consts::PI;
        let clipped = angle_restrict_symm(theta);
        assert!((clipped - 0.5 * std::f64::consts::PI).abs() < 1e-10);
    }
    #[test]
    fn test_angle_restrict_pos_e() {
        let theta = 2.5 * std::f64::consts::PI;
        let clipped = angle_restrict_pos_e(theta);
        assert!((clipped.val - 0.5 * std::f64::consts::PI).abs() < 1e-10);
    }
    #[test]
    fn test_angle_restrict_pos() {
        let theta = 2.5 * std::f64::consts::PI;
        let clipped = angle_restrict_pos(theta);
        assert!((clipped - 0.5 * std::f64::consts::PI).abs() < 1e-10);
    }
    #[test]
    fn test_polar_to_rect() {
        let r = 10.0_f64;
        let theta = std::f64::consts::PI;
        let z = polar_to_rect(r, theta);
        assert!((z.val.re - r * theta.cos()).abs() < 1e-10);
        assert!((z.val.im - r * theta.sin()).abs() < 1e-10);
    }
    #[test]
    fn test_rect_to_polar() {
        let z = Complex::new(10.0, -4.0);

        let (r, t) = rect_to_polar(z);
        assert!((r.val - z.re.hypot(z.im)).abs() < 1e-10);
        assert!((t.val - z.arg()).abs() < 1e-10);
    }
    #[test]
    fn test_sin_err_e() {
        let x = 1.0_f64;
        let dx = 0.0_f64;
        assert!((sin_err_e(x, dx).val - 0.84147098480789651).abs() < 1e-10);
    }
    #[test]
    fn test_cos_err_e() {
        let x = 1.0_f64;
        let dx = 0.0_f64;
        assert!((cos_err_e(x, dx).val - 0.54030230586813972).abs() < 1e-10);
    }
    #[test]
    fn test_sinc_e() {
        let x = 1.5_f64;
        assert!((sinc_e(x).val - -0.21220659078919378).abs() < 1e-10)
    }
}
