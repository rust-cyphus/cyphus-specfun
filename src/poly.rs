use num::Float;

pub struct Polynomial<T: Float> {
    coef: Vec<T>,
}

impl<T: Float> Polynomial<T> {
    pub fn new(coef: Vec<T>) -> Polynomial<T> {
        Polynomial { coef: coef.clone() }
    }
    /// Evaluate a polynomial at a point.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::poly::Polynomial;
    /// let poly = Polynomial::new([1.0,2.0,3.0].to_vec());
    /// assert_eq!(poly.eval(-1.0), 2.0);
    /// ```
    pub fn eval(&self, x: T) -> T {
        // Make sure that there is at least one element in the coeff list. If not,
        // return zero.
        match self.coef.last() {
            Some(cn) => {
                let mut z = *cn;
                // Construct the result in Horner form.
                for val in self.coef.iter().rev().skip(1) {
                    z = z.mul_add(x, *val);
                }
                z
            }
            None => T::zero(),
        }
    }
    /// Compute all derivatives of a polynomial at a point (including the
    /// zeroth derivative.)
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::poly::Polynomial;
    /// let poly = Polynomial::new([1.0,2.0,3.0].to_vec());
    /// assert_eq!(poly.derivs(-1.0), vec![2.0, -4.0, 6.0]);
    /// ```
    pub fn derivs(&self, x: T) -> Vec<T> {
        let mut n = 0;
        let mut nmax = 0;
        let lenc = self.coef.len();
        let mut res = vec![T::zero(); lenc];

        for i in 0..lenc {
            if n < lenc {
                res[i] = self.coef[lenc - 1];
                nmax = n;
                n += 1;
            } else {
                res[i] = T::zero();
            }
        }

        for i in 0..(lenc - 1) {
            let k = (lenc - 1) - i;
            res[0] = x.mul_add(res[0], self.coef[k - 1]);
            let lmax = if nmax < k { nmax } else { k - 1 };
            for l in 1..(lmax + 1) {
                res[l] = x.mul_add(res[l], res[l - 1]);
            }
        }

        let mut f = T::one();
        for i in 2..(nmax + 1) {
            f = f * T::from(i).unwrap();
            res[i] = res[i] * f;
        }
        res
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_poly_eval_f64() {
        let poly = Polynomial::new([1.0, 2.0, 3.0].to_vec());
        let x = -1.0;
        println!("{:?}", poly.eval(x));
        println!("{:?}", poly.derivs(x));
        assert_eq!(poly.derivs(-1.0), vec![2.0, -4.0, 6.0]);
    }
    #[test]
    fn test_poly_eval_f32() {
        let poly = Polynomial::new([1f32, 2.0, 3.0].to_vec());
        let x = -1.0;
        println!("{:?}", poly.eval(x));
        println!("{:?}", poly.derivs(x));
        assert_eq!(poly.derivs(-1.0), vec![2.0, -4.0, 6.0]);
    }
}
