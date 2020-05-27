pub struct Polynomial {
    coef: Vec<f64>,
}

impl Polynomial {
    pub fn new(coef: Vec<f64>) -> Polynomial {
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
    pub fn eval(&self, x: f64) -> f64 {
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
            None => 0.0,
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
    pub fn derivs(&self, x: f64) -> Vec<f64> {
        let mut n = 0;
        let mut nmax = 0;
        let lenc = self.coef.len();
        let mut res = vec![0.0; lenc];

        for i in 0..lenc {
            if n < lenc {
                res[i] = self.coef[lenc - 1];
                nmax = n;
                n += 1;
            } else {
                res[i] = 0.0;
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

        let mut f = 1.0;
        for i in 2..(nmax + 1) {
            f *= i as f64;
            res[i] *= f;
        }
        res
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_poly_eval() {
        let poly = Polynomial::new([1.0, 2.0, 3.0].to_vec());
        let x = -1.0;
        println!("{:?}", poly.eval(x));
        println!("{:?}", poly.derivs(x));
        assert_eq!(poly.derivs(-1.0), vec![2.0, -4.0, 6.0]);
    }
}
