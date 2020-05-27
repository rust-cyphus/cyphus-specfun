pub(crate) mod bessel_data;
mod bessel_helpers;
pub(crate) mod cyl_bessel_j;
pub(crate) mod cyl_bessel_k;
pub(crate) mod cyl_bessel_y;
mod olver;

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
    fn cyl_bessel_j0_e(&self) -> SpecFunResult<Self>
    where
        Self: num::Num;
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
    fn cyl_bessel_j1_e(&self) -> SpecFunResult<Self>
    where
        Self: num::Num;
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
    fn cyl_bessel_jn_e(&self, n: i32) -> SpecFunResult<Self>
    where
        Self: num::Num;
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
    fn cyl_bessel_jv_e(&self, nu: Self) -> SpecFunResult<Self>
    where
        Self: num::Num;
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
    fn cyl_bessel_y0_e(&self) -> SpecFunResult<Self>
    where
        Self: num::Num;
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
    fn cyl_bessel_y1_e(&self) -> SpecFunResult<Self>
    where
        Self: num::Num;
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
    fn cyl_bessel_yn_e(&self, n: i32) -> SpecFunResult<Self>
    where
        Self: num::Num;
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
    fn cyl_bessel_yv_e(&self, nu: Self) -> SpecFunResult<Self>
    where
        Self: num::Num;
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
    use crate::bessel::cyl_bessel_y::*;
    use crate::consts::SQRT_DLB_EPS;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_utils::*;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
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
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (0.0001, 1.0),
            0.7652115411876708497,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (0.0001, 10.0),
            -0.2459270166445205,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (0.0009765625, 10.0),
            -0.2458500798634692,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (0.75, 1.0),
            0.5586524932048917478,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (0.75, 10.0),
            -0.04968928974751508135,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (1.0, 0.001),
            0.0004999999375000026,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (1.0, 1.0),
            0.4400505857449335160,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (1.75, 1.0),
            0.168593922545763170103,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (30.0, 1.0),
            3.482869794251482902e-42,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (30.0, 100.0),
            0.08146012958117222297,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (10.0, 1.0),
            2.6306151236874532070e-10,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (10.0, 100.0),
            -0.05473217693547201474,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (10.2, 100.0),
            -0.03548919161046526864,
            TOL2,
            SpecFunCode::Success
        );

        /* related to BUG#3 problem */
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (2.0, 900.0),
            -0.019974345269680646400,
            TOL4,
            SpecFunCode::Success
        );
        // FIXME: Test fails
        //test_check_result_and_code!(cyl_bessel_jv_e, (2.0, 15000.0), -0.0020455820181216382666,TOL4, SpecFunCode::Success);

        /* Jnu for negative integer nu */

        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-13.0, 4.3966088395743909188e-1),
            -4.4808688873945295916e-19,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-43.0, 3.2953184511924683369e-1),
            -3.4985116870357066158e-87,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-8.0, 3.5081119129046332101e-1),
            2.2148387185334257376e-11,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-17.0, 1.6717234250796879247e-1),
            -1.3336696368048418208e-33,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-20.0, 1.0085435516296562067e-1),
            4.6463938513369299663e-45,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 1.0423882905853389231e-1),
            -7.0211488877617616588e-69,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-10.0, 2.014635069337132169e-1),
            2.9614218635575180136e-17,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-23.0, 3.3470953660313309239e-1),
            -5.3786691977970313226e-41,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-13.0, 1.796487688043502395e-2),
            -3.9796275104402232636e-37,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-30.0, 4.3851505455005663023e-1),
            6.3723816221242007043e-53,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-31.0, 2.9134673992671269075e-1),
            -1.4108323502121482121e-60,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-28.0, 7.6696967407852829575e-1),
            7.2211466723046636958e-42,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-6.0, 1.9829576302527197691e-2),
            1.3193561299244410609e-15,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 4.1513019703157674847e-1),
            4.362149422492811209e-45,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-19.0, 1.3033500482689031834e-2),
            -2.4071043287214877014e-59,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.0, 7.2050154387915780891e-1),
            6.8377203227324865568e-4,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-36.0, 9.058436592111083308e-1),
            1.1063535538518134862e-54,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-37.0, 4.4893454860671838425e-2),
            -7.1436404620419702996e-105,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-25.0, 1.4253284590439372148e-1),
            -1.3532982149810467304e-54,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-35.0, 6.8075538970323047806e-1),
            -4.0087398199591044218e-57,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-9.0, 6.9978823056579836732e-2),
            -2.1657760307485265867e-19,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-30.0, 8.3903642499214990225e-1),
            1.8046736761082165751e-44,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-43.0, 6.543891152683782553e-1),
            -2.2618181452671187564e-74,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.0, 2.557673774862201033e-1),
            -1.2684081505646766845e-1,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-21.0, 2.3961944728579204361e-1),
            -8.7037441094405713341e-40,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-13.0, 3.7732868931080794955e-1),
            -6.1458679923678486123e-20,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-33.0, 3.8478064679256785588e-1),
            -2.7471851206170506902e-61,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-50.0, 8.3045728127372512419e-1),
            2.6904315506244885678e-84,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-33.0, 3.4391840270211572065e-1),
            -6.7604474418522862789e-63,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-23.0, 4.4699489157302132678e-1),
            -4.1674520864249550489e-38,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-23.0, 8.6894302936681315656e-1),
            -1.8080036072447165775e-31,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 2.3997647040322801696e-1),
            1.2775332144705955465e-46,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.0, 9.7688590680055575385e-1),
            2.260676935437026869e-3,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-20.0, 8.9090003293395823104e-1),
            3.8525697232471553813e-26,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-37.0, 9.4057990133775869779e-1),
            -5.4483632006108552908e-56,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.0, 2.6641070579761597867e-1),
            -1.3202706692617745567e-1,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-7.0, 9.8878610880180753494e-1),
            -1.3892626129288932172e-6,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-42.0, 3.1042396387236895292e-1),
            7.4291433893935351027e-86,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-18.0, 6.8722270094498481586e-1),
            6.9210950943508574463e-25,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-16.0, 1.6113707818439830316e-2),
            1.5066992425677873707e-47,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-31.0, 3.0720392023079679622e-1),
            -7.2938615192106070364e-60,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-8.0, 3.4957196590626821859e-1),
            2.153068469114375124e-11,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.0, 9.9805150837159144851e-1),
            2.4578749047519916152e-3,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-8.0, 8.4538260733781281419e-1),
            2.4776262290872395403e-8,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-21.0, 7.8935512795174026579e-1),
            -6.4516652247824702208e-29,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-32.0, 3.4859401756229795902e-1),
            1.9985100102959434086e-60,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-28.0, 9.2530929240293814377e-1),
            1.3798987861611642315e-39,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-32.0, 5.7276907866586953775e-1),
            1.5890407703711501803e-53,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-16.0, 8.6364032267280037696e-1),
            6.9097436254460436398e-20,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 8.1065458967451038451e-2),
            6.2316950805227520253e-58,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-17.0, 6.946415511078842723e-1),
            -4.3495405983099048944e-23,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-14.0, 2.7619565270505476723e-1),
            1.0511271746088180311e-23,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-7.0, 8.579791601466317093e-2),
            -5.3039175204559594309e-14,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-50.0, 7.4681490716561041968e-1),
            1.3331741219055661824e-86,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-30.0, 9.6631424224612452556e-1),
            1.2465815577059566852e-42,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-14.0, 2.6745819874317492788e-1),
            6.7024985048974981757e-24,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-22.0, 3.6947309321414497419e-1),
            6.4975710352073715961e-38,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-49.0, 2.3389602803779083655e-2),
            -3.5281224467005794333e-158,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-48.0, 7.7482504299905127354e-1),
            1.3711857024107478291e-81,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-30.0, 6.5337950709230049969e-1),
            9.9749347191786840601e-48,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 5.4117488569210959438e-1),
            -3.8843890347204927665e-48,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-50.0, 9.4038774092645791075e-1),
            1.3455212768163778034e-81,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.0, 5.627399919447310892e-1),
            -2.703780920058947009e-1,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-35.0, 7.9925999507829895011e-2),
            -1.1060656306690664139e-89,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-10.0, 2.8875439067692705135e-3),
            1.0844833648855202142e-35,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-31.0, 2.7645446745852278572e-1),
            -2.7740931384652548129e-61,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-12.0, 8.5374968606378578252e-1),
            7.5366644001760905687e-14,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-11.0, 2.2997159465452075155e-1),
            -1.1626579415654480052e-18,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-47.0, 5.8616043503503357359e-1),
            -3.4270544836018351049e-85,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 8.4884849730214263521e-1),
            1.8679667366331791405e-33,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 5.3943176031162412158e-1),
            3.5300994115064965228e-38,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.0, 9.3387179309343112648e-1),
            -1.6062668879215016378e-2,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-20.0, 5.5061917642540109891e-1),
            2.5629166990385734238e-30,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 2.9315675715515822781e-1),
            5.1505442739915338391e-49,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-6.0, 8.9079715253593128884e-1),
            1.0539758479683914316e-5,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-30.0, 9.1220250697048036241e-2),
            2.2291985671052015068e-73,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-49.0, 3.2934552845357396486e-1),
            -6.7628009346158039786e-102,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-20.0, 6.2008856084878849548e-3),
            2.7691703032688977815e-69,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-23.0, 7.1620177578706903472e-1),
            -2.124755495743639839e-33,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 5.1124413385957329246e-1),
            -7.462204354283390559e-49,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-23.0, 6.7405321505832900803e-1),
            -5.2689769204867120958e-34,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-31.0, 7.6828173704641609615e-2),
            -1.600162678935095954e-78,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.0, 7.8546641070654814244e-1),
            9.610534504452577567e-4,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-39.0, 7.8708523270381710591e-1),
            -7.8237015591341025437e-63,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-8.0, 1.1973778137874968338e-1),
            4.0918170073170305674e-15,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-22.0, 7.9699790929323855652e-1),
            1.4309765990568672215e-30,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-11.0, 4.4800161319259052168e-1),
            -1.7773520138988139169e-15,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-39.0, 9.2637892745877041811e-1),
            -4.4944192865215894733e-60,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 4.6475617572623309956e-1),
            -4.7026155154566357963e-50,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-25.0, 4.6435749048262147537e-1),
            -8.9952935327704632698e-42,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-28.0, 9.6511079209639008185e-1),
            4.4842768640542298547e-39,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-36.0, 3.8411791760130458261e-2),
            4.296403338839098526e-104,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-6.0, 2.891145018145052606e-1),
            1.2636012998902815302e-8,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-43.0, 3.0866754076112185885e-1),
            -2.1010663709383681337e-88,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-44.0, 8.3972445789951961023e-1),
            9.7813493638738341846e-72,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-22.0, 8.2307293342676789357e-1),
            2.9041436661554638719e-30,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-12.0, 9.8857250941935924585e-1),
            4.357505306871049656e-13,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-47.0, 9.9901153675544108048e-1),
            -2.6090278787375958472e-74,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-44.0, 2.8286777063559952163e-1),
            1.5832896249242322418e-92,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-10.0, 9.7638349947419439863e-1),
            2.0735913941162063742e-10,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-36.0, 4.6452631984005212745e-1),
            4.0190510760125038996e-65,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-25.0, 8.2954737598010312336e-1),
            -1.7855017422981225812e-35,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 4.2326977999795604681e-1),
            -3.1249531389372439734e-51,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-9.0, 6.339167576518494852e-1),
            -8.8082994072251150317e-11,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-27.0, 7.4535258299077637954e-1),
            -2.4368032520208918805e-40,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-45.0, 8.0998976704970278418e-1),
            -1.8024726219976677441e-74,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 5.3437414610693284794e-1),
            2.8159486268019058346e-38,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-19.0, 3.6987646459831328184e-1),
            -9.7200835901252686409e-32,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-49.0, 4.7308684164506190021e-1),
            -3.4438316393709775799e-94,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-42.0, 8.4616849424460901479e-1),
            1.4469107537397962381e-67,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-39.0, 8.6819804427642418616e-1),
            -3.5837055310896411954e-61,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-28.0, 4.1490163307722590922e-1),
            2.448154551562710141e-49,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-25.0, 6.312651668273163642e-1),
            -1.9374739456138691106e-38,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-32.0, 5.3525470388026220677e-1),
            1.8191207037881755277e-54,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-13.0, 9.2605599515408535679e-3),
            -7.2212651931099339154e-41,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-37.0, 4.8798589647640992562e-1),
            -1.5614996577959031121e-66,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-37.0, 2.2096551301466541766e-2),
            -2.9050812034192451665e-116,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-42.0, 1.8092873560926168207e-1),
            1.0575388628113044972e-95,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-22.0, 4.2954143428969324083e-1),
            1.7857221060719548205e-36,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-41.0, 1.7704659657689619579e-1),
            -2.017674513721498041e-93,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.0, 8.8755155653793053987e-1),
            -1.3862799262620708632e-2,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.0, 6.7617657704638521874e-1),
            -3.1913053561511127823e-1,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-34.0, 6.4464304038438048178e-1),
            6.4622314654558520448e-56,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-41.0, 7.2625991432244163042e-1),
            -2.7337118024791538344e-68,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.0, 6.357184822423154937e-1),
            -5.2186206749718741014e-3,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-6.0, 8.4999778579632777264e-1),
            7.9757193518631922586e-6,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-47.0, 4.9270801771231759268e-1),
            -9.7743102162756354643e-89,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-44.0, 3.0382431918975667824e-1),
            3.6749344196700669304e-91,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-21.0, 8.3709075395163832807e-1),
            -2.2120291813240017972e-28,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-13.0, 2.2926361048823468174e-1),
            -9.4684470339645238166e-23,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.0, 7.7595731076113971064e-1),
            7.1557648504788709828e-2,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-13.0, 8.3155666144468474085e-1),
            -1.760187305034807398e-15,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 8.9229820511590331545e-1),
            1.8952192929209682492e-36,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-9.0, 4.3091918823985208523e-1),
            -2.7448142505229531657e-12,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-45.0, 7.6232573842954975111e-1),
            -1.177044451508791199e-75,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-47.0, 8.5370192783467777094e-1),
            -1.6168449888151601463e-77,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-17.0, 7.0917926374340919579e-1),
            -6.1835780259077271289e-23,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-48.0, 5.3177634068823620691e-1),
            1.9544175507861216812e-89,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-20.0, 8.6236563774769292261e-1),
            2.0102896278761019978e-26,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-47.0, 1.9398034049650767996e-1),
            -9.1928516850183758066e-108,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-22.0, 1.9059407555598432968e-1),
            3.0818305203000064476e-44,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-28.0, 1.177497340192600765e-1),
            1.1853471604859571152e-64,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-44.0, 5.4411719229619710879e-1),
            5.0099597095462268336e-80,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-22.0, 8.4096736569723091858e-1),
            4.6598891827094030103e-30,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-14.0, 5.7336031747513286455e-1),
            2.8892068362904543886e-19,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 4.3701161637218456507e-1),
            2.2572337854881401663e-40,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-16.0, 5.5950502438849852065e-1),
            6.6952184074205206877e-23,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-15.0, 7.2135709344094709909e-1),
            -1.724936657760223936e-19,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-42.0, 8.9503358993826252397e-1),
            1.5285906153951174029e-66,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-34.0, 8.8842662562391397862e-1),
            3.5115599466639988263e-51,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 7.1820429322328243568e-1),
            3.3906734873299682373e-35,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-42.0, 6.9983283703407621113e-1),
            4.9840979255532732037e-71,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 3.4956023377394930027e-1),
            -1.2169918834421082937e-53,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 3.1780802172157483391e-1),
            1.0816638066322224274e-43,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-42.0, 6.6747313617191922586e-1),
            6.8258422034194326368e-72,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 9.8281798409069473598e-1),
            -1.2641970706335426485e-40,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-48.0, 9.3933236478420901552e-1),
            1.4126090419557425431e-77,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-33.0, 7.9969562813573452357e-1),
            -8.3521436185818779821e-51,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.0, 1.749832037277650746e-1),
            2.4377505307831309992e-6,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-12.0, 7.4928517324325015606e-1),
            1.578984980873870348e-14,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-50.0, 6.8668643284952988717e-1),
            2.0060641365741957198e-88,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 4.628453813124784986e-1),
            7.3802979097358757591e-44,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-17.0, 1.3398796901656815413e-1),
            -3.1005014564142463477e-35,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-21.0, 4.3589737907768555406e-1),
            -2.4897620016130572781e-34,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-49.0, 9.3471214811043877923e-1),
            -1.0635172631183545319e-79,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-12.0, 6.2842914244476056474e-1),
            1.919020727030761691e-15,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-44.0, 3.9902397275715061045e-1),
            5.9381636725852262522e-86,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-38.0, 7.4844742354073981342e-1),
            1.1452741249870059158e-61,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-22.0, 1.7024143837455386901e-1),
            2.5696433432212310518e-45,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-35.0, 2.4663778723047911984e-1),
            -1.4846380265631517846e-72,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-42.0, 3.4675532474808813305e-1),
            7.7576502065450687145e-84,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-42.0, 9.7015007293565236242e-1),
            4.5073367850753509233e-65,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-18.0, 5.8064934044500015204e-1),
            3.3392049842730194442e-26,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.0, 8.9712869996285071984e-1),
            1.6201373653351794808e-3,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-14.0, 4.7509593095794841351e-1),
            2.0819777331649343154e-20,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 1.9971440277743206076e-1),
            1.5567772398263651029e-48,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.0, 3.1326041354180815314e-1),
            1.2166506426643152762e-2,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-21.0, 7.1780328035546876532e-1),
            -8.7820775346034440302e-30,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-40.0, 3.4783624449708255548e-1),
            5.0330128895858604038e-79,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 4.6025606383259078342e-1),
            7.8278566159609528116e-40,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-37.0, 6.8329209514074967672e-1),
            -4.0018049034521909865e-61,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 3.1644447969459523952e-1),
            3.757803139189680586e-48,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-49.0, 3.8814659649103508174e-1),
            -2.1178321069354253493e-98,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-33.0, 2.2586340634075651258e-1),
            -6.3690605699693325702e-69,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-31.0, 1.1833425544176035201e-1),
            -1.0457450400633015896e-72,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-7.0, 4.716233225505345007e-1),
            -7.9892591292002198427e-9,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-28.0, 4.0216249780484400656e-1),
            1.0224868057823447281e-49,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-39.0, 2.1728515555094074309e-1),
            -1.2424793343150735922e-84,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-10.0, 1.5286805658821812372e-1),
            1.8744497113230339685e-18,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.0, 4.2012489177724585853e-1),
            2.1740379601921820516e-2,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.0, 5.4168160735476556955e-1),
            -3.2509626190739798323e-3,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 6.0999547254418081401e-1),
            9.6515399655293906821e-41,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-31.0, 9.3599472437451867441e-1),
            -7.236873645285246215e-45,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.0, 8.9238535456317991508e-2),
            9.9477909077321557346e-4,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-41.0, 3.3286697432119768766e-1),
            -3.5168501713472066379e-82,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 1.3103200887095798302e-2),
            4.1630610278945554747e-84,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-7.0, 6.8545576155223653312e-1),
            -1.0860095433456484207e-7,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-33.0, 7.4597656684747976078e-1),
            -8.4232256181114982289e-52,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-35.0, 9.5265851504353628581e-1),
            -5.1260638475279101845e-52,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-5.0, 1.9993324513780069188e-2),
            -8.319296787329444617e-13,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-15.0, 7.291071285552115835e-2),
            -2.0411952376466778385e-34,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-5.0, 4.307852573603263607e-1),
            -3.8336545021181925733e-6,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-8.0, 3.0719264454074178501e-1),
            7.6627991262305533713e-12,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-7.0, 2.9261260328577001029e-1),
            -2.8395431884068098274e-10,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-36.0, 3.4285828911893011822e-1),
            7.1807133181285014617e-70,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-22.0, 2.1687887653368606307e-1),
            5.2860475778514174667e-43,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 7.2816755037040510323e-1),
            4.7187086299885949165e-35,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-18.0, 2.0826276232560462604e-2),
            3.2368741843295077202e-52,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-9.0, 6.8082174052201672454e-1),
            -1.6719854980400483279e-10,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 1.1998114825675920571e-1),
            4.2119340347581322841e-59,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-9.0, 5.9197600088654039875e-1),
            -4.7631865156190005935e-11,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-30.0, 1.2367288101522705215e-1),
            2.0588316029270585207e-69,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-7.0, 8.3266930505292536647e-1),
            -4.2096524602233328394e-7,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-18.0, 4.360196631312459384e-1),
            1.9281550661128359168e-28,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-13.0, 9.8501660515145040901e-1),
            -1.5833136710018445626e-14,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 9.3281324180154613247e-1),
            -2.7828395455119501545e-41,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-50.0, 8.9831019278310796215e-1),
            1.3646576617083900982e-82,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-34.0, 8.1153956230382506493e-1),
            1.6192058088789772833e-52,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-15.0, 9.5908894233909785027e-1),
            -1.2293883538807523551e-17,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-28.0, 4.7478353398916835208e-1),
            1.0667274838088242221e-47,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-26.0, 3.1425431663890729964e-1),
            3.137014371489532261e-48,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 6.8861856868877100233e-1),
            -4.2000859317520628674e-45,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-16.0, 6.9807655407582355921e-1),
            2.3026948238804970982e-21,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-15.0, 1.9223628937777433793e-1),
            -4.2201734817670464106e-28,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-24.0, 7.91209811831343471e-1),
            3.458241440092889033e-34,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-12.0, 2.8881796002183274727e-1),
            1.7143390913163291276e-19,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-34.0, 3.6891378721647167497e-1),
            3.7139793083014182422e-64,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-20.0, 8.4841223828616526898e-1),
            1.4510812436551651903e-26,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-23.0, 2.2490195812594682131e-1),
            -5.7468688920782767025e-45,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-18.0, 2.2504144134842978217e-1),
            1.3048322081397375779e-33,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-41.0, 8.9085721717774902491e-1),
            -1.1841910084346823163e-64,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-14.0, 3.5776817082613807574e-1),
            3.9325021938284721675e-22,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-17.0, 4.6898364389788035003e-1),
            -5.492570150236103145e-26,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-6.0, 7.4085998179632632531e-1),
            3.5186865149767756957e-6,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-29.0, 8.1247678941673114604e-1),
            -5.0783189409391835819e-43,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-6.0, 1.7590874156867732351e-2),
            6.4299450335557031571e-16,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-6.0, 4.1122931368227349961e-2),
            1.0494595145859932593e-13,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-18.0, 3.4357780610013843947e-2),
            2.6519427947417673311e-48,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-8.0, 7.2901630769663700817e-1),
            7.6159005881302088369e-9,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-40.0, 6.2434787548617066655e-1),
            7.297739135890827566e-69,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-37.0, 2.5717302633809380684e-1),
            -7.9671811532819427071e-77,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-12.0, 4.4871088635019795091e-1),
            3.3823657137507787902e-17,
            TOL2,
            SpecFunCode::Success
        );

        /* Jnu for negative nu */

        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.7408502287167772557e+1, 4.3891036254061936271e-2),
            -1.5118966152679114528e+42,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.0369289750688261206e+1, 8.617077861907621132e-1),
            1.3458137581188368176e+61,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.6398263821779503158, 5.7108775125890870648e-1),
            -1.1073338178875224514e+2,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.0743837616631487936e+1, 2.2372123946196517647e-1),
            1.3987244312547157439e+37,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.7297952956642177128e+1, 2.8440280848533973972e-1),
            -5.5832331287880973953e+27,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.3507072230824139103e+1, 9.7008738913549467403e-1),
            -9.9108981827284991851e+12,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.4663106025697138284e+1, 9.6739483983540323655e-1),
            2.5305841722999151766e+67,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.7450767905208691173e+1, 5.2288309478685515663e-3),
            -3.4541920228396234708e+180,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.8736914274754280581e+1, 4.4318942692773250336e-1),
            1.2783291424623089209e+27,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.4560791710276042462e+1, 5.6595841783863809163e-1),
            1.7364781221343685309e+56,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.8723165418996042381e+1, 2.4003201116391659975e-1),
            8.229650479070536446e+54,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.1780816480084454714e+1, 3.5000033756039229564e-1),
            -8.9158096963672457852e+56,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-8.3572537806135923827, 8.994859317735841446e-1),
            2.4471474432717596765e+6,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.9179499652113791027e+1, 7.3466044408596358731e-1),
            -1.0446080588162613503e+82,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.0204535104461498585e+1, 1.3316368076269799096e-1),
            1.6723180404777533538e+93,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.9597169419874779916e+1, 8.3895848736941180651e-1),
            -1.9885394381212418599e+80,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.2228113163441444643e+1, 7.2360070820350912315e-1),
            3.7033741801434815187e+12,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.3252347726870680973e+1, 4.7042383176138740544e-2),
            -2.2524439080867705956e+121,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.8363275112470146292e+1, 3.5406545166014987925e-1),
            7.0915928368505654149e+95,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.4031313414732363799e+1, 1.7077658851102128503e-1),
            4.2707681524978432343e+46,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.4994220322161396103e+1, 7.6004498361697560048e-2),
            8.3491267575601512811e+85,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.5842451866339926986e+1, 7.1655842470093644641e-1),
            -1.4956016465164596038e+18,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.3303969013414183336e+1, 3.4261981308825774017e-1),
            -1.7313464383524562686e+84,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.9620214546900609668e+1, 8.9559048893071335969e-2),
            -6.5439278427877993956e+69,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.9267775155252921091e+1, 4.9173394917277714574e-1),
            -2.7879726255003962141e+68,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.085805422677201981e+1, 7.1844728658692273953e-2),
            1.7330833141098585591e+21,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.3205222720657600449e+1, 2.0938011160695718605e-1),
            -1.2855953290893319419e+93,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.1761832688338363731e-1, 4.0692821730479031008e-1),
            1.0616114810207300625,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.4176152886105712673e+1, 7.3083748089885469319e-1),
            2.3708170326600879156e+51,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.7495316392965243753e+1, 3.6374757654471397858e-1),
            -2.4050181588419272908e+26,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.4363999308889822481e+1, 2.2964635433363436307e-1),
            1.4858445128296594446e+23,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.5945782535187393546e+1, 5.5004988492253222851e-1),
            -5.3196903529172738965e+19,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.6539439893219625339e+1, 5.4022213494661757887e-1),
            3.4606719889912371036e+40,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.9104819423206076833e+1, 8.7581079115207012305e-1),
            -8.3806822204633705862e+57,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.1344422530419629852e+1, 1.8412292138063948156e-1),
            -1.3032526695489281999e+18,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.2479004807998664153e+1, 7.628052499161305046e-1),
            3.8593605090529013702e+67,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.7234292208462683278e+1, 5.6929354834881763134e-1),
            -1.3560087920173394791e+41,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.9852923491811760298e+1, 2.1764339448558498796e-2),
            -3.1065720979758225192e+88,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.9470737712526065987e+1, 9.1397839227827761003e-1),
            -4.9854244578384505794e+39,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.4879459756439547254e+1, 8.7694253508024822462e-1),
            4.0540656704233640216e+31,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-8.6562240771489400173e-1, 6.4882619355798588927e-1),
            9.5827819637209987022e-2,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.4096305256724256786e+1, 1.1837901256079790395e-1),
            1.5389662008468777023e+26,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.3739240782592000796e+1, 1.1830951764837136725e-1),
            -5.4851415830067607572e+25,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-7.7683272384817803185e-1, 2.0897180603218726575e-1),
            1.3452855819917342033,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.9007825894244022613e+1, 8.0804305816981973945e-1),
            -1.9558153171413640836e+78,
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.224872379992519031, 9.4716012050013901355e-1),
            -8.7507643021583199242e-1,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-5.0807600398060325807e+1, 6.1902119795226148946e-1),
            1.9558680407989173708e+89,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.2130610423140110162e+1, 7.2184881647444607846e-1),
            3.0709609117301881893e+67,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.1186300874143057513e+1, 1.3944550368322800884e-1),
            -1.2405415201132534983e+95,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.1777582295815773824e+1, 7.6178874271077561762e-1),
            -7.1748063501973138583e+27,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.2185912810368133222e+1, 3.959510541558671016e-1),
            1.196451653185591802e+56,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.2976006083946226988e+1, 4.5739390828369379657e-1),
            1.0599129365585919882e+77,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.0412232215064606945e+1, 6.7510758855896918145e-1),
            2.4302601636670462267e+45,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.8388188389747281955e+1, 8.9233979909724013718e-1),
            9.1410432331502484426e+20,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-2.4298840569257253701e+1, 6.8042862360863559591e-1),
            4.0995648850574613979e+33,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.8834735272673504063e+1, 4.2324469410479691518e-1),
            7.0355133597135130631e+69,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.6091840934431606344e+1, 1.7508918790902548661e-1),
            8.7456315137421979067e+103,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.9754264394942756728, 1.558717798933954111e-2),
            -3.551027943747170162e+2,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.7168594342366560374e+1, 6.976373865080374073e-1),
            -1.1036447967023475572e+58,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.5379944292029245754e+1, 5.3150471205968938472e-2),
            -1.0469743921754287032e+126,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-1.153181439556533474e+1, 1.8204366094125528818e-1),
            -4.0986515168430307785e+18,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-4.9939680334734766385e+1, 8.132277926604580844e-1),
            -9.5179038651143567503e+80,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-3.5986707652967270442e+1, 5.5665988728905782182e-1),
            -1.27797927382078249e+58,
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_jv_e,
            (-6.7046620273111013832, 1.059530133767196237e-1),
            3.8106055649273069958e+10,
            TOL2,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_cyl_bessel_y0_e() {
        test_check_result_and_code!(
            cyl_bessel_y0_e,
            (0.1),
            -1.5342386513503668441,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_y0_e,
            (2.0),
            0.5103756726497451196,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_y0_e,
            (256.0),
            -0.03381290171792454909,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_y0_e,
            (4294967296.0),
            3.657903190017678681e-06,
            SQRT_TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_cyl_bessel_y1_e() {
        test_check_result_and_code!(
            cyl_bessel_y1_e,
            (0.1),
            -6.45895109470202698800,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_y1_e,
            (2.0),
            -0.10703243154093754689,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_y1_e,
            (100.0),
            -0.020372312002759793305,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_y1_e,
            (4294967296.0),
            0.000011612249378370766284,
            TOL4,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_cyl_bessel_yn_e() {
        test_check_result_and_code!(
            cyl_bessel_yn_e,
            (4, 0.1),
            -305832.29793353160319,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_yn_e,
            (5, 2.0),
            -9.935989128481974981,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_yn_e,
            (100, 100.0),
            -0.16692141141757650654,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            cyl_bessel_yn_e,
            (100, 4294967296.0),
            3.657889671577715808e-06,
            SQRT_TOL0,
            SpecFunCode::Success
        );
        //FIXME: This test fails
        test_check_result_and_code!(
            cyl_bessel_yn_e,
            (1000, 4294967296.0),
            3.656551321485397501e-06,
            2.0e-05,
            SpecFunCode::Success
        );
    }
}
