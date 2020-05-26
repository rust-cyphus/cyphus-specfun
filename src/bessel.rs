pub(crate) mod cyl_bessel_j;
pub(crate) mod cyl_bessel_y;
pub(crate) mod bessel_data;
mod bessel_helpers;
mod olver;

use crate::bessel::cyl_bessel_j::cyl_bessel_j0_e;
use crate::result::SpecFunResult;

use cyl_bessel_j::*;
use cyl_bessel_y::*;


pub trait CylBesselJ {
    /// Compute the cylindrical Bessel-function of the first kind of order 0
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselJ;
    /// assert!((1f64.cyl_bessel_j0_e().val - 0.7651976865579666).abs() < 1e-10);
    /// ```
    fn cyl_bessel_j0_e(&self) -> SpecFunResult<Self> where Self: num::Num;
    /// Compute the cylindrical Bessel-function of the first kind of order 0.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselJ;
    /// assert!((1f64.cyl_bessel_j0() - 0.7651976865579666).abs() < 1e-10);
    /// ```
    fn cyl_bessel_j0(&self) -> Self;
    /// Compute the cylindrical Bessel-function ofof the first kind  order 1
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselJ;
    /// assert!((1f64.cyl_bessel_j1_e().val - 0.4400505857449335).abs() < 1e-10);
    /// ```
    fn cyl_bessel_j1_e(&self) -> SpecFunResult<Self> where Self: num::Num;
    /// Compute the cylindrical Bessel-function of the first kind of order 1.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselJ;
    /// assert!((1f64.cyl_bessel_j1() - 0.4400505857449335).abs() < 1e-10);
    /// ```
    fn cyl_bessel_j1(&self) -> Self;
    /// Compute the cylindrical Bessel-function of the first kind of order n
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselJ;
    /// assert!((1f64.cyl_bessel_jn_e(5).val - 0.002476638964109955).abs() < 1e-10);
    /// ```
    fn cyl_bessel_jn_e(&self, n: i32) -> SpecFunResult<Self> where Self: num::Num;
    /// Compute the cylindrical Bessel-function of the first kind of order n.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselJ;
    /// assert!((1f64.cyl_bessel_jn(5) - 0.002476638964109955).abs() < 1e-10);
    /// ```
    fn cyl_bessel_jn(&self, n: i32) -> Self;
    /// Compute the cylindrical Bessel-function of the first kind of order nu
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselJ;
    /// assert!((1f64.cyl_bessel_jv_e(0.5).val - 0.6713967071418031).abs() < 1e-10);
    /// ```
    fn cyl_bessel_jv_e(&self, nu: Self) -> SpecFunResult<Self> where Self: num::Num;
    /// Compute the cylindrical Bessel-function of the first kind of order nu.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselJ;
    /// assert!((1f64.cyl_bessel_jv(0.5) - 0.6713967071418031).abs() < 1e-10);
    /// ```
    fn cyl_bessel_jv(&self, nu: Self) -> Self;
}

impl CylBesselJ for f64 {
    fn cyl_bessel_j0_e(&self) -> SpecFunResult<Self> {
        cyl_bessel_j0_e(*self)
    }
    fn cyl_bessel_j0(&self) -> Self {
        cyl_bessel_j0_e(*self).val
    }
    fn cyl_bessel_j1_e(&self) -> SpecFunResult<Self> {
        cyl_bessel_j1_e(*self)
    }
    fn cyl_bessel_j1(&self) -> Self {
        cyl_bessel_j0_e(*self).val
    }
    fn cyl_bessel_jn_e(&self, n: i32) -> SpecFunResult<Self> {
        cyl_bessel_jn_e(n, *self)
    }
    fn cyl_bessel_jn(&self, n: i32) -> Self {
        cyl_bessel_jn_e(n, *self).val
    }
    fn cyl_bessel_jv_e(&self, nu: Self) -> SpecFunResult<Self> {
        cyl_bessel_jv_e(nu, *self)
    }
    fn cyl_bessel_jv(&self, nu: Self) -> Self {
        cyl_bessel_jv_e(nu, *self).val
    }
}

pub trait CylBesselY {
    /// Compute the cylindrical Bessel-function of the second kind of order 0
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselY;
    /// assert!((1f64.cyl_bessel_y0_e().val - 0.08825696421567696).abs() < 1e-10);
    /// ```
    fn cyl_bessel_y0_e(&self) -> SpecFunResult<Self> where Self: num::Num;
    /// Compute the cylindrical Bessel-function of the second kind of order 0.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselY;
    /// assert!((1f64.cyl_bessel_y0() - 0.08825696421567696).abs() < 1e-10);
    /// ```
    fn cyl_bessel_y0(&self) -> Self;
    /// Compute the cylindrical Bessel-function of the second kind of order 1
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselY;
    /// assert!((1f64.cyl_bessel_y1_e().val + 0.7812128213002887).abs() < 1e-10);
    /// ```
    fn cyl_bessel_y1_e(&self) -> SpecFunResult<Self> where Self: num::Num;
    /// Compute the cylindrical Bessel-function of the second kind of order 1
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselY;
    /// assert!((1f64.cyl_bessel_y1() + 0.7812128213002887).abs() < 1e-10);
    /// ```
    fn cyl_bessel_y1(&self) -> Self;
    /// Compute the cylindrical Bessel-function of the second kind of order n
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselY;
    /// assert!((1f64.cyl_bessel_yn_e(5).val + 260.4058666258122).abs() < 1e-10);
    /// ```
    fn cyl_bessel_yn_e(&self, n: i32) -> SpecFunResult<Self> where Self: num::Num;
    /// Compute the cylindrical Bessel-function of the second kind of order n.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselY;
    /// assert!((1f64.cyl_bessel_yn(5) + 260.4058666258122).abs() < 1e-10);
    /// ```
    fn cyl_bessel_yn(&self, n: i32) -> Self;
    /// Compute the cylindrical Bessel-function of the second kind of order nu
    /// along with an error estimate.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselY;
    /// assert!((1f64.cyl_bessel_yv_e().val + 0.4310988680183761).abs() < 1e-10);
    /// ```
    fn cyl_bessel_yv_e(&self, nu: Self) -> SpecFunResult<Self> where Self: num::Num;
    /// Compute the cylindrical Bessel-function of the second kind of order nu.
    ///
    /// # Examples
    /// ```rust
    /// # use cyphus_specfun::bessel::CylBesselY;
    /// assert!((1f64.cyl_bessel_yv() + 0.4310988680183761).abs() < 1e-10);
    /// ```
    fn cyl_bessel_yv(&self, nu: Self) -> Self;
}

impl CylBesselY for f64 {
    fn cyl_bessel_y0_e(&self) -> SpecFunResult<Self> {
        cyl_bessel_y0_e(*self)
    }
    fn cyl_bessel_y0(&self) -> Self {
        cyl_bessel_y0_e(*self).val
    }
    fn cyl_bessel_y1_e(&self) -> SpecFunResult<Self> {
        cyl_bessel_y1_e(*self)
    }
    fn cyl_bessel_y1(&self) -> Self {
        cyl_bessel_y1_e(*self).val
    }
    fn cyl_bessel_yn_e(&self, n: i32) -> SpecFunResult<Self> {
        cyl_bessel_yn_e(n, *self)
    }
    fn cyl_bessel_yn(&self, n: i32) -> Self {
        cyl_bessel_yn_e(n, *self).val
    }
    fn cyl_bessel_yv_e(&self, nu: Self) -> SpecFunResult<Self> {
        cyl_bessel_yv_e(nu, *self)
    }
    fn cyl_bessel_yv(&self, nu: Self) -> Self {
        cyl_bessel_yv_e(nu, *self).val
    }
}

#[cfg(test)]
mod test {
    use crate::bessel::cyl_bessel_j::*;
    use crate::consts::SQRT_DLB_EPS;
    use crate::test_utils::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL4: f64 = 16384.0 * f64::EPSILON;

    #[test]
    fn test_cyl_bessel_j0_e() {
        test_check_result_and_code!(
            cyl_bessel_j0_e,
            (0.1),
            0.997_501_562_066_04,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_j0_e,
            (2.0),
            0.223_890_779_141_235_67,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_j0_e,
            (100.0),
            0.019_985_850_304_223_122,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_j0_e,
            (1.0e+10),
            2.175_591_750_246_892e-6,
            SQRT_TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_cyl_bessel_j1_e() {
        test_check_result_and_code!(
            cyl_bessel_j1_e,
            (0.1),
            0.049_937_526_036_242,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_j1_e,
            (2.0),
            0.576_724_807_756_873_4,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_j1_e,
            (100.0),
            -0.077_145_352_014_112_16,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_j1_e,
            (1.0e+10),
            -7.676_508_175_684_158e-6,
            TOL4,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_cyl_bessel_jn_e() {
        test_check_result_and_code!(
            cyl_bessel_jn_e,
            (4, 0.1),
            2.602_864_854_568_403e-7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jn_e,
            (5, 2.0),
            0.007_039_629_755_871_685,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jn_e,
            (10, 20.0),
            0.186_482_558_023_945_1,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jn_e,
            (100, 100.0),
            0.096_366_673_295_861_56,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_cyl_bessel_jn_e_long() {
        test_check_result_and_code!(
            cyl_bessel_jn_e,
            (45, 900.0),
            0.02562434700634278108,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_cyl_bessel_jv_e() {
        test_check_result_and_code!(cyl_bessel_jv_e, (0.0001, 1.0),         0.7652115411876708497,  TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (0.0001, 10.0),       -0.2459270166445205,     TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (0.0009765625, 10.0), -0.2458500798634692,     TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (0.75, 1.0),           0.5586524932048917478,  TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (0.75, 10.0),         -0.04968928974751508135, TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, ( 1.0, 0.001), 0.0004999999375000026,     TOL0, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, ( 1.0,   1.0), 0.4400505857449335160,     TOL0, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, ( 1.75,  1.0), 0.168593922545763170103,     TOL1, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (30.0,   1.0), 3.482869794251482902e-42,  TOL0, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (30.0, 100.0), 0.08146012958117222297,    TOL1, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (10.0,   1.0), 2.6306151236874532070e-10, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (10.0, 100.0), -0.05473217693547201474,   TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (10.2, 100.0), -0.03548919161046526864,   TOL2, SpecFunCode::Success);
    }
}
