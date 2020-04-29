//! Module for implementing logarithms with real and complex arguments as well
//! as the associated errors.
//!
//! # Example
//! ```
//! ```

use crate::cheb::cheb_eval_e;
use crate::consts::{ROOT5_DBL_EPS, ROOT6_DBL_EPS};
use crate::result::{SpecFunCode, SpecFunResult};
use num::Complex;

pub trait Logarithm {
    /// Compute the natural logarithm of a number and the associated error.
    ///
    /// # Example
    /// ```
    /// let x = 1.0;
    ///
    /// assert!(x.ln_e().val.abs() < 1e-10);
    /// ```
    fn ln_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute the natural log of the absolute value of a number and the
    /// associated error (value equivalent to x.abs().ln()).
    ///
    /// # Example
    /// ```
    /// let x = -1.0;
    ///
    /// assert!(x.ln_abs_e().val.abs() < 1e-10);
    /// ```
    fn ln_abs_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute the natural log of the absolute value of a
    /// number
    ///
    /// # Example
    /// ```
    /// let x = -1.0;
    ///
    /// assert!(x.ln_abs().abs() < 1e-10);
    /// ```
    fn ln_abs(&self) -> Self
    where
        Self: Sized;
    /// Compute the natural log of (1 + x) along with error
    /// estimate.
    ///
    /// # Example
    /// ```
    /// let res = 0.0_f64.ln_p1_e();
    ///
    /// assert!(res.abs() < 1e-10);
    /// ```
    fn ln_p1_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute the natural log(1 + x) - x along with error
    /// estimate.
    ///
    /// # Example
    /// ```
    /// let res = 10.0_f64.ln_p1_mx_e();
    ///
    /// assert!((res.val + 7.6021047272016295).abs() < 1e-10);
    /// ```
    fn ln_p1_mx_e(&self) -> SpecFunResult<Self>
    where
        Self: Sized;
    /// Compute the natural log(1 + x) - x.
    ///
    /// # Example
    /// ```
    /// let res = 10.0_f64.ln_p1_mx();
    ///
    /// assert!((res + 7.6021047272016295).abs() < 1e-10);
    /// ```
    fn ln_p1_mx(&self) -> Self
    where
        Self: Sized;
}

impl Logarithm for f64 {
    fn ln_e(&self) -> SpecFunResult<Self> {
        ln_e(*self)
    }
    fn ln_abs_e(&self) -> SpecFunResult<Self> {
        ln_abs_e(*self)
    }
    fn ln_abs(&self) -> Self {
        ln_abs_e(*self).val
    }
    fn ln_p1_e(&self) -> SpecFunResult<Self> {
        ln_p1_e(*self)
    }
    fn ln_p1_mx_e(&self) -> SpecFunResult<Self> {
        ln_p1_mx_e(*self)
    }
    fn ln_p1_mx(&self) -> Self {
        ln_p1_mx_e(*self).val
    }
}

// --------------
// ---- Data ----
// --------------

// Chebyshev expansion for log(1 + x(t))/x(t)
// x(t) = (4t-1)/(2(4-t))
// t(x) = (8x+1)/(2(x+2))
// -1/2 < x < 1/2
// -1 < t < 1
const LOPX_DATA: [f64; 21] = [
    2.166_479_106_643_952_70_f64,
    -0.285_653_985_510_497_42_f64,
    0.015_177_672_556_905_53_f64,
    -0.002_002_159_049_414_15_f64,
    0.000_192_113_751_640_56_f64,
    -0.000_025_532_588_861_05_f64,
    2.900_451_266_040_062e-06_f64,
    -3.887_381_351_705_734e-07_f64,
    4.774_367_872_940_045e-08_f64,
    -6.450_196_977_609_031e-09_f64,
    8.275_197_662_881_238e-10_f64,
    -1.126_049_937_649_204e-10_f64,
    1.484_457_669_227_093e-11_f64,
    -2.032_851_597_246_211e-12_f64,
    2.729_123_122_054_921e-13_f64,
    -3.758_197_783_038_793e-14_f64,
    5.110_734_587_086_167e-15_f64,
    -7.072_215_001_143_327e-16_f64,
    9.708_975_832_824_846e-17_f64,
    -1.349_263_745_752_193e-17_f64,
    1.865_732_791_067_729e-18_f64,
];

// Chebyshev expansion for (log(1 + x(t)) - x(t))/x(t)^2
//
// x(t) = (4t-1)/(2(4-t))
// t(x) = (8x+1)/(2(x+2))
// -1/2 < x < 1/2
// -1 < t < 1
const LOPXMX_DATA: [f64; 20] = [
    -1.121_002_313_237_441_03_f64,
    0.195_534_627_733_793_86_f64,
    -0.014_674_704_538_080_83_f64,
    0.001_666_782_504_743_65_f64,
    -0.000_185_433_561_477_00_f64,
    0.000_022_801_540_217_71_f64,
    -2.803_125_311_663_352e-06_f64,
    3.593_656_887_252_216e-07_f64,
    -4.624_185_704_106_206e-08_f64,
    6.082_263_745_940_399e-09_f64,
    -8.033_982_442_481_579e-10_f64,
    1.075_171_827_749_937e-10_f64,
    -1.444_531_091_422_461e-11_f64,
    1.957_391_218_061_033e-12_f64,
    -2.661_443_679_679_306e-13_f64,
    3.640_263_431_526_958e-14_f64,
    -4.993_749_592_275_500e-15_f64,
    6.880_289_021_884_680e-16_f64,
    -9.503_412_979_480_427e-17_f64,
    1.317_013_501_305_099e-17_f64,
];

/// Compute the natural logarithm of a number and the associated error.
///
/// # Example
/// ```
/// let x = 1.0;
///
/// assert!(ln_e(x).val.abs() < 1e-10);
/// ```
pub fn ln_e(x: f64) -> SpecFunResult<f64> {
    if x <= 0.0 {
        let res = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        res.issue_warning("ln_e", &[x]);
        res
    } else {
        let val = x.ln();
        SpecFunResult {
            val: val,
            err: 2.0 * f64::EPSILON * val.abs(),
            code: SpecFunCode::Success,
        }
    }
}

/// Compute the natural log of the absolute value of a number and the
/// associated error (value equivalent to x.abs().ln()).
///
/// # Example
/// ```
/// let x = -1.0;
///
/// assert!(ln_abs_e(x).val.abs() < 1e-10);
/// ```
pub fn ln_abs_e(x: f64) -> SpecFunResult<f64> {
    if x.abs() < f64::EPSILON {
        let res = SpecFunResult {
            val: f64::NAN,
            err: f64::NAN,
            code: SpecFunCode::DomainErr,
        };
        res.issue_warning("ln_abs_e", &[x]);
        res
    } else {
        let val = x.abs().ln();
        SpecFunResult {
            val: val,
            err: 2.0 * f64::EPSILON * val.abs(),
            code: SpecFunCode::Success,
        }
    }
}

/// Compute the natural log of a complex number z as (ln(|z|, theta) where:
/// ln(z) = ln(|z|) + i theta
///
/// # Example
/// ```
/// use num::Complex;
///
/// let z = Complex::new(2.0, -1.0);
/// let (lnr, theta) = complex_ln_e(z);
///
/// assert!((lnr.val - 0.92873123336475721).abs() < 1e-10);
/// assert!((theta.val +0.52270629950815609).abs() < 1e-10);
/// ```
pub fn complex_ln_e(z: Complex<f64>) -> SpecFunResult<Complex<f64>> {
    let mut result = SpecFunResult {
        val: Complex::new(0.0, 0.0),
        err: Complex::new(0.0, 0.0),
        code: SpecFunCode::Success,
    };

    if z.re.abs() > f64::EPSILON && z.im.abs() > f64::EPSILON {
        let ax = z.re.abs();
        let ay = z.im.abs();
        let min = ax.min(ay);
        let max = ax.max(ay);
        result.val.re = max.ln() + 0.5 * (1.0 + (min / max) * (min / max)).ln();
        result.err.re = 2.0 * f64::EPSILON * result.val.re.abs();
        result.val.im = z.im.atan2(z.re);
        result.err.im = f64::EPSILON * result.val.re.abs();
        result
    } else {
        result.val = Complex::new(f64::NAN, f64::NAN);
        result.err = Complex::new(f64::NAN, f64::NAN);
        result.code = SpecFunCode::DomainErr;
        result.issue_warning("complex_ln_e", &[z]);
        result
    }
}

/// Compute the natural log of (1 + x).
///
/// # Example
/// ```
/// let res = ln_p1_e(0.0);
///
/// assert!(res.val.abs() < 1e-10);
/// ```
pub fn ln_p1_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if x <= -1.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result.issue_warning("ln_p1_e", &[x]);
        result
    } else if x.abs() < ROOT6_DBL_EPS {
        let c1: f64 = -0.5;
        let c2: f64 = 1.0 / 3.0;
        let c3: f64 = -1.0 / 4.0;
        let c4: f64 = 1.0 / 5.0;
        let c5: f64 = -1.0 / 6.0;
        let c6: f64 = 1.0 / 7.0;
        let c7: f64 = -1.0 / 8.0;
        let c8: f64 = 1.0 / 9.0;
        let c9: f64 = -1.0 / 10.0;

        let t = c9
            .mul_add(x, c8)
            .mul_add(x, c7)
            .mul_add(x, c6)
            .mul_add(x, c5);

        result.val = t
            .mul_add(x, c4)
            .mul_add(x, c3)
            .mul_add(x, c2)
            .mul_add(x, c1)
            .mul_add(x, 1.0)
            * x;
        result.err = f64::EPSILON * result.val.abs();
        result
    } else if x.abs() < 0.5 {
        let t = 0.5 * x.mul_add(8.0, 1.0) / (x + 2.0);
        let c = cheb_eval_e(t, &LOPX_DATA, -1.0, 1.0);
        result.val = x * c.val;
        result.err = (x * c.err).abs();
        result
    } else {
        result.val = (1.0 + x).ln();
        result.err = f64::EPSILON * result.val.abs();
        result
    }
}

/// Compute the natural log(1 + x) - x.
///
/// # Example
/// ```
/// let res = ln_p1_mx_e(10.0);
///
/// assert!((res.val + 7.6021047272016295).abs() < 1e-10);
/// ```
pub fn ln_p1_mx_e(x: f64) -> SpecFunResult<f64> {
    let mut result = SpecFunResult {
        val: 0.0,
        err: 0.0,
        code: SpecFunCode::Success,
    };

    if x <= -1.0 {
        result.val = f64::NAN;
        result.err = f64::NAN;
        result.code = SpecFunCode::DomainErr;
        result.issue_warning("ln_p1_e", &[x]);
        result
    } else if x.abs() < ROOT5_DBL_EPS {
        let c1: f64 = -0.5;
        let c2: f64 = 1.0 / 3.0;
        let c3: f64 = -1.0 / 4.0;
        let c4: f64 = 1.0 / 5.0;
        let c5: f64 = -1.0 / 6.0;
        let c6: f64 = 1.0 / 7.0;
        let c7: f64 = -1.0 / 8.0;
        let c8: f64 = 1.0 / 9.0;
        let c9: f64 = -1.0 / 10.0;
        let t = c9
            .mul_add(x, c8)
            .mul_add(x, c7)
            .mul_add(x, c6)
            .mul_add(x, c5);
        result.val = t
            .mul_add(x, c4)
            .mul_add(x, c3)
            .mul_add(x, c2)
            .mul_add(x, c1)
            * x
            * x;
        result.err = f64::EPSILON * result.val.abs();
        result
    } else if x.abs() < 0.5 {
        let t = 0.5 * x.mul_add(8.0, 1.0) / (x + 2.0);
        let c = cheb_eval_e(t, &LOPXMX_DATA, -1.0, 1.0);
        result.val = x * x * c.val;
        result.err = x * x * c.err;
        result
    } else {
        let lterm = (1.0 + x).ln();
        result.val = lterm - x;
        result.err = f64::EPSILON * (lterm.abs() + x.abs());
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ln_e() {
        let x = 1.0;
        assert!(ln_e(x).val.abs() < 1e-10);
    }
    #[test]
    fn test_ln_abs_e() {
        let x = -1.0;
        assert!(ln_abs_e(x).val.abs() < 1e-10);
    }
    #[test]
    fn test_complex_ln_e() {
        let z1 = Complex::new(2.0, -1.0);
        let z2 = Complex::new(10.0, 1.0);
        let z3 = Complex::new(-4.0, 1.0);
        let ln1 = complex_ln_e(z1);
        let ln2 = complex_ln_e(z2);
        let ln3 = complex_ln_e(z3);
        assert!((ln1.val.re - 0.80471895621705019).abs() < 1e-10);
        assert!((ln1.val.im + 0.4636476090008061).abs() < 1e-10);

        assert!((ln2.val.re - 2.3075602584206297).abs() < 1e-10);
        assert!((ln2.val.im - 0.0996686524911620).abs() < 1e-10);

        assert!((ln3.val.re - 1.4166066720281080).abs() < 1e-10);
        assert!((ln3.val.im - 2.8966139904629291).abs() < 1e-10);
    }
    #[test]
    fn test_ln_p1_e() {
        let res = ln_p1_e(9.0);
        assert!((res.val - 2.3025850929940457).abs() < 1e-10);
    }
    #[test]
    fn test_ln_p1_mx() {
        let res = ln_p1_mx_e(10.0);
        assert!((res.val + 7.6021047272016295).abs() < 1e-10);
    }
}
