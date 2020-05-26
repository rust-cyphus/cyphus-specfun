//! Chebyshev Polynomials
//!
//! Contains implementations for evaluating Chebyshev series.
//!
//! # Examples
//! ```
//! let x = 1.0;
//! let coeffs = [1.0, 2.0, 3.0];
//! assert!((cheb_eval_e(x, &coeffs, -1.0, 1.0).val - 6.0).abs() < 1e-10);
//! ```

use crate::result::{SpecFunCode, SpecFunResult};
use num::Float;

/// Evaluate a Chebyshev series at `x` given coefficients and domain boundaries.
///
/// # Notes
/// For a function f(x), the Nth order Chebyshev approximation to f is given by:
///     f(x) = -1/2 c[0] + Sum[c[k] * T[k,x] ,{k, 0, N-1}]
/// where c[k] is the kth coefficient and T[k,x] is the kth Chebyshev polynomial
/// evaluated at x. The coefficients should be given by:
///     c[k] = 2/N Sum[f[Cos[pi*(k-1/2)/N]] * Cos[pi*j*(k-1/2)/N],{k,1,N}]
///
/// # Examples
/// 6th order Chebyshev approximation to sin(x):
/// ```
/// let x = 0.5;
/// let coeffs = [0.,0.880101,0.,-0.0391267,0.,0.00050252];
/// let sinhalf = 0.5_f64.sin();
/// let cheb = cheb_eval_e(x, &coeffs, -1.0, 1.0).val;
///
/// assert!((sinhalf - cheb).abs() < 1e-5);
/// ```
pub fn cheb_eval_e<T: Float>(x: T, coeffs: &[T], lb: T, ub: T) -> SpecFunResult<T> {
    let mut d = T::zero();
    let mut dd = T::zero();
    let two = T::from(2).unwrap();

    let y = (two * x - lb - ub) / (ub - lb);
    let y2 = two * y;

    let mut e = T::zero();

    for j in (1..(coeffs.len())).rev() {
        let temp = d;
        d = y2 * d - dd + coeffs[j];
        e = e + (y2 * temp).abs() + dd.abs() + coeffs[j].abs();
        dd = temp;
    }

    {
        let temp = d;
        d = y * d - dd + coeffs[0] / two;
        e = e + (y * temp).abs() + dd.abs() + T::from(0.5).unwrap() * coeffs[0].abs();
    }

    let val = d;
    let err = T::epsilon() * e + coeffs.last().unwrap().abs();
    let code = SpecFunCode::Success;
    SpecFunResult { val, err, code }
}

/// Data for a Chebyshev series over a given interval
pub struct ChebSeries<T: Float> {
    /// Coefficients
    pub coeffs: Vec<T>,
    /// Lower interval point
    pub a: T,
    /// Upper interval point
    pub b: T,
}

impl<T: Float> ChebSeries<T> {
    pub fn new(coeffs: Vec<T>, a: T, b: T) -> ChebSeries<T> {
        ChebSeries { coeffs, a, b }
    }
    /// Evaluate the Chebyshev series at x and return result and error.
    pub fn eval(&self, x: T) -> SpecFunResult<T> {
        let mut d = T::zero();
        let mut dd = T::zero();
        let two = T::from(2).unwrap();
        let lb: T = self.a;
        let ub: T = self.b;

        let y = (two * x - lb - ub) / (ub - lb);
        let y2 = two * y;

        let mut e = T::zero();

        for j in (1..(self.coeffs.len())).rev() {
            let temp = d;
            d = y2 * d - dd + self.coeffs[j];
            e = e + (y2 * temp).abs() + dd.abs() + self.coeffs[j].abs();
            dd = temp;
        }

        {
            let temp = d;
            d = y * d - dd + self.coeffs[0] / two;
            e = e + (y * temp).abs() + dd.abs() + T::from(0.5).unwrap() * self.coeffs[0].abs();
        }

        let val = d;
        let err = T::epsilon() * e + self.coeffs.last().unwrap().abs();
        let code = SpecFunCode::Success;
        SpecFunResult { val, err, code }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cheby() {
        let x = 0.5;
        let coeffs = [0., 0.880_101, 0., -0.039_126_7, 0., 0.000_502_52];
        let sinhalf = 0.5_f64.sin();
        let cheb = cheb_eval_e(x, &coeffs, -1.0, 1.0).val;
        assert!((sinhalf - cheb).abs() < 1e-5);

        let cheb = ChebSeries {
            coeffs: vec![0.0, 0.880_101, 0.0, -0.039_126_7, 0.0, 0.000_502_52],
            a: -1.0,
            b: 1.0,
        };
        assert!((sinhalf - cheb.eval(x).val).abs() < 1e-5);
    }
}
