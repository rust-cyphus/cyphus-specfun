mod data;
pub(crate) mod inc;
pub(crate) mod mono;
pub(crate) mod poch;
pub(crate) mod psi;

use crate::result::SpecFunResult;

//use inc::*;
use mono::*;
use poch::*;
use psi::*;

/// Implementation for the gamma function
pub trait Gamma {
    /// Gamma-function with error estimate
    fn gamma_e(&self) -> SpecFunResult<f64>;
    /// Gamma-function
    fn gamma(&self) -> f64;
    /// Compute gamma* with error information
    fn gammastar_e(&self) -> SpecFunResult<f64>;
    /// Compute gamma*
    fn gammastar(&self) -> f64;
    /// Compute 1 / gamma(x) with error information
    fn gammainv_e(&self) -> SpecFunResult<f64>;
    /// Compute 1 / gamma(x)
    fn gammainv(&self) -> f64;
    /// Gamma-function with error estimate
    fn lngamma_e(&self) -> SpecFunResult<f64>;
    /// Gamma-function
    fn lngamma(&self) -> f64;
    // Compute the digamma function of a number along with the error.
    ///
    /// ```
    /// use cyphus::specfun::gamma::Gamma;
    ///
    /// let res = 2.0_f64.digamma_e();
    /// let abserr = (res.val + 0.5772156649).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn digamma_e(&self) -> SpecFunResult<f64>;
    /// Compute the digamma function of a number.
    ///
    /// ```
    /// use cyphus::specfun::gamma::Gamma;
    ///
    /// let res = 2.0_f64.digamma_e();
    /// let abserr = (res.val + 0.5772156649).abs();
    /// assert!(abserr < 1e-10);
    /// assert!(abserr < res.err);
    /// ```
    fn digamma(&self) -> f64;
    /// Compute the trigamma function of a number along with the error.
    ///
    /// ```
    /// use cyphus::specfun::gamma::Gamma;
    ///
    /// let res = 2.0_f64.trigamma_e();
    /// let abserr = (res.val - 1.644934067).abs();
    /// assert!(abserr < 1e-10);
    /// assert!(abserr < res.err);
    /// ```
    fn trigamma_e(&self) -> SpecFunResult<f64>;
    /// Compute the trigamma function of a number.
    ///
    /// ```
    /// use cyphus::specfun::gamma::Gamma;
    ///
    /// let res = 2.0_f64.trigamma();
    /// let abserr = (res - 1.644934067).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn trigamma(&self) -> f64;
    /// Compute the polygamma function of a number along with the error.
    ///
    /// ```
    /// use cyphus::specfun::gamma::Gamma;
    ///
    /// let res = 2.0_f64.polygamma_e(5);
    /// let abserr = (res.val - 122.0811674).abs();
    /// assert!(abserr < 1e-10);
    /// assert!(abserr < res.err);
    /// ```
    fn polygamma_e(&self, n: usize) -> SpecFunResult<f64>;
    /// Compute the polygamma function of a number.
    ///
    /// ```
    /// use cyphus::specfun::gamma::Gamma;
    ///
    /// let res = 2.0_f64.polygamma_e(5);
    /// let abserr = (res.val - 122.0811674).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn polygamma(&self, n: usize) -> f64;
    /// Compute the Pochhammer symbol (self)_x = gamma(self + x) / gamma(self)
    /// along with the associated error.
    /// ```
    /// use cyphus::specfun::gamma::Gamma;
    /// let a:f64 = 10.0;
    /// let x:f64 = -1.0;
    /// // Compute (10)_{-1} = 1/9
    /// assert!((a.poch_e(x).val - 1.0/9.0) < 1e-10);
    /// ```
    fn poch_e(&self, x: f64) -> SpecFunResult<f64>;
    /// Compute the Pochhammer symbol (self)_x = gamma(self + x) / gamma(self).
    /// ```
    /// use cyphus::specfun::gamma::Gamma;
    /// let a:f64 = 10.0;
    /// let x:f64 = -1.0;
    /// // Compute (10)_{-1} = 1/9
    /// assert!((a.poch(x) - 1.0/9.0) < 1e-10);
    /// ```
    fn poch(&self, x: f64) -> f64;
    /// Compute the relative Pochhammer symbol ((self)_x - 1) / x with the
    /// associated error.
    fn pochrel_e(&self, x: f64) -> SpecFunResult<f64>;
    /// Compute the relative Pochhammer symbol ((self)_x - 1) / x.
    fn pochrel(&self, x: f64) -> f64;
    /// Compute the nautral log of the Pochhammer symbol along with the
    /// associated error.
    fn lnpoch_e(&self, x: f64) -> SpecFunResult<f64>;
    /// Compute the nautral log of the Pochhammer symbol.
    fn lnpoch(&self, x: f64) -> f64;
}

macro_rules! impl_gamma_int {
    ($T:ty) => {
        impl Gamma for $T {
            fn gamma_e(&self) -> SpecFunResult<f64> {
                fact_e(*self as usize)
            }
            fn gamma(&self) -> f64 {
                fact_e(*self as usize).val
            }
            fn gammastar_e(&self) -> SpecFunResult<f64> {
                gammastar_e(*self as f64)
            }
            fn gammastar(&self) -> f64 {
                gammastar_e(*self as f64).val
            }
            fn gammainv_e(&self) -> SpecFunResult<f64> {
                gammainv_e(*self as f64)
            }
            fn gammainv(&self) -> f64 {
                gammainv_e(*self as f64).val
            }
            fn lngamma_e(&self) -> SpecFunResult<f64> {
                lnfact_e(*self as usize)
            }
            fn lngamma(&self) -> f64 {
                lnfact_e(*self as usize).val
            }
            fn digamma_e(&self) -> SpecFunResult<f64> {
                digamma_int_e(*self as usize)
            }
            fn digamma(&self) -> f64 {
                digamma_int_e(*self as usize).val
            }
            fn trigamma_e(&self) -> SpecFunResult<f64> {
                trigamma_int_e(*self as usize)
            }
            fn trigamma(&self) -> f64 {
                trigamma_int_e(*self as usize).val
            }
            fn polygamma_e(&self, n: usize) -> SpecFunResult<f64> {
                polygamma_e(n, *self as f64)
            }
            fn polygamma(&self, n: usize) -> f64 {
                polygamma_e(n, *self as f64).val
            }
            fn poch_e(&self, x: f64) -> SpecFunResult<f64> {
                poch_e(*self as f64, x as f64)
            }
            fn poch(&self, x: f64) -> f64 {
                poch_e(*self as f64, x as f64).val
            }
            fn pochrel_e(&self, x: f64) -> SpecFunResult<f64> {
                pochrel_e(*self as f64, x as f64)
            }
            fn pochrel(&self, x: f64) -> f64 {
                pochrel_e(*self as f64, x as f64).val
            }
            fn lnpoch_e(&self, x: f64) -> SpecFunResult<f64> {
                lnpoch_e(*self as f64, x as f64)
            }
            fn lnpoch(&self, x: f64) -> f64 {
                lnpoch_e(*self as f64, x as f64).val
            }
        }
    };
}

macro_rules! impl_gamma_float {
    ($T:ty) => {
        impl Gamma for $T {
            fn gamma_e(&self) -> SpecFunResult<f64> {
                gamma_e(*self as f64)
            }
            fn gamma(&self) -> f64 {
                gamma_e(*self as f64).val
            }
            fn gammastar_e(&self) -> SpecFunResult<f64> {
                gammastar_e(*self as f64)
            }
            fn gammastar(&self) -> f64 {
                gammastar_e(*self as f64).val
            }
            fn gammainv_e(&self) -> SpecFunResult<f64> {
                gammainv_e(*self as f64)
            }
            fn gammainv(&self) -> f64 {
                gammainv_e(*self as f64).val
            }
            fn lngamma_e(&self) -> SpecFunResult<f64> {
                lngamma_e(*self as f64)
            }
            fn lngamma(&self) -> f64 {
                lngamma_e(*self as f64).val
            }
            fn digamma_e(&self) -> SpecFunResult<f64> {
                digamma_e(*self as f64)
            }
            fn digamma(&self) -> f64 {
                digamma_e(*self as f64).val
            }
            fn trigamma_e(&self) -> SpecFunResult<f64> {
                trigamma_e(*self as f64)
            }
            fn trigamma(&self) -> f64 {
                trigamma_e(*self as f64).val
            }
            fn polygamma_e(&self, n: usize) -> SpecFunResult<f64> {
                polygamma_e(n, *self as f64)
            }
            fn polygamma(&self, n: usize) -> f64 {
                polygamma_e(n, *self as f64).val
            }
            fn poch_e(&self, x: f64) -> SpecFunResult<f64> {
                poch_e(*self as f64, x as f64)
            }
            fn poch(&self, x: f64) -> f64 {
                poch_e(*self as f64, x as f64).val
            }
            fn pochrel_e(&self, x: f64) -> SpecFunResult<f64> {
                pochrel_e(*self as f64, x as f64)
            }
            fn pochrel(&self, x: f64) -> f64 {
                pochrel_e(*self as f64, x as f64).val
            }
            fn lnpoch_e(&self, x: f64) -> SpecFunResult<f64> {
                lnpoch_e(*self as f64, x as f64)
            }
            fn lnpoch(&self, x: f64) -> f64 {
                lnpoch_e(*self as f64, x as f64).val
            }
        }
    };
}

impl_gamma_int!(i8);
impl_gamma_int!(i16);
impl_gamma_int!(i32);
impl_gamma_int!(i64);
impl_gamma_int!(i128);
impl_gamma_int!(u8);
impl_gamma_int!(u16);
impl_gamma_int!(u32);
impl_gamma_int!(usize);
impl_gamma_int!(u128);

impl_gamma_float!(f32);
impl_gamma_float!(f64);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_check_result_sgn_and_code;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL5: f64 = 131_072.0 * f64::EPSILON;

    #[test]
    fn test_ln_gamma_e() {
        test_check_result_and_code!(
            lngamma_e,
            (-0.1),
            2.368_961_332_728_788_6,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-1.0 / 256.0),
            5.547_444_766_967_471,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (1.0e-08),
            18.420_680_738_180_21,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (0.1),
            2.252_712_651_734_205,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (1.0 + 1.0 / 256.0),
            -0.002_242_222_659_961_150_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (2.0 + 1.0 / 256.0),
            0.001_656_417_755_696_172_8,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (100.0),
            359.134_205_369_575_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-1.0 - 1.0 / 65536.0),
            11.090_348_438_090_048,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-1.0 - 1.0 / 268_435_456.0),
            19.408_121_054_103_475,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-100.5),
            -364.900_968_309_427_36,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-100.0 - 1.0 / 65536.0),
            -352.649_091_011_709_8,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_lnfact_e() {
        test_check_result_and_code!(lnfact_e, (0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(lnfact_e, (1), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            lnfact_e,
            (7),
            8.525_161_361_065_415,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnfact_e,
            (33),
            85.054_467_017_581_52,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_gamma_e() {
        test_check_result_and_code!(
            gamma_e,
            (1.0 + 1.0 / 4096.0),
            0.999_859_137_145_940_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (1.0 + 1.0 / 32.0),
            0.982_901_099_283_626_9,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (2.0 + 1.0 / 256.0),
            1.001_657_790_373_358_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(gamma_e, (9.0), 40320.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gamma_e, (10.0), 362_880.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            gamma_e,
            (100.0),
            9.332_621_544_394_415e155,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (170.0),
            4.269_068_009_004_705e304,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (171.0),
            7.257_415_615_307_999e306,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (-10.5),
            -2.640_121_820_547_716e-7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (-11.25),
            6.027_393_816_261_932e-8,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (-1.0 + 1.0 / 65536.0),
            -65_536.422_805_878_19,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_gammastar_e() {
        test_check_result_and_code!(
            gammastar_e,
            (1.0e-08),
            3_989.423_555_759_891,
            2.0 * TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (1.0e-05),
            126.171_684_698_826_91,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (0.001),
            12.708_492_464_364_074,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (1.5),
            1.056_344_244_268_559_8,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (3.0),
            1.028_064_517_918_789_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (9.0),
            1.009_298_426_421_819,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (11.0),
            1.007_602_428_310_496_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (100.0),
            1.000_833_677_872_012_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (1.0e+05),
            1.000_000_833_333_680_5,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(gammastar_e, (1.0e+20), 1.0, TOL0, SpecFunCode::Success);
    }

    #[test]
    fn test_gammainv_e() {
        test_check_result_and_code!(gammainv_e, (1.0), 1.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gammainv_e, (2.0), 1.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gammainv_e, (3.0), 0.5, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gammainv_e, (4.0), 1.0 / 6.0, TOL0, SpecFunCode::Success);

        test_check_result_and_code!(
            gammainv_e,
            (10.0),
            1.0 / 362_880.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammainv_e,
            (100.0),
            1.071_510_288_125_467e-_156,
            TOL2,
            SpecFunCode::Success
        );

        test_check_result_and_code!(gammainv_e, (0.0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gammainv_e, (-1.0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gammainv_e, (-2.0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gammainv_e, (-3.0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gammainv_e, (-4.0), 0.0, TOL0, SpecFunCode::Success);

        test_check_result_and_code!(
            gammainv_e,
            (-10.5),
            -1.0 / 2.640_121_820_547_716e-7,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammainv_e,
            (-11.25),
            1.0 / 6.027_393_816_261_932e-8,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammainv_e,
            (-1.0 + 1.0 / 65536.0),
            -1.0 / 65_536.422_805_878_19,
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_lnpoch_e() {
        test_check_result_and_code!(lnpoch_e, (5.0, 0.0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            lnpoch_e,
            (5.0, 1.0 / 65536.0),
            0.000_022_981_557_571_259_39,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnpoch_e,
            (5.0, 1.0 / 256.0),
            0.005_884_960_217_985_189,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnpoch_e,
            (7.0, 3.0),
            6.222_576_268_071_369,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnpoch_e,
            (5.0, 2.0),
            3.401_197_381_662_155_5,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_lnpoch_sgn_e() {
        test_check_result_sgn_and_code!(
            lnpoch_sgn_e,
            (5.0, 0.0),
            0.0,
            TOL0,
            1.0,
            SpecFunCode::Success
        );
        test_check_result_sgn_and_code!(
            lnpoch_sgn_e,
            (-4.5, 0.25),
            0.743_011_647_511_992,
            TOL1,
            1.0,
            SpecFunCode::Success
        );
        test_check_result_sgn_and_code!(
            lnpoch_sgn_e,
            (-4.5, 1.25),
            2.189_930_630_448_317_5,
            TOL2,
            -1.0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_poch_e() {
        dbg!(f64::NEG_INFINITY.exp());
        test_check_result_and_code!(poch_e, (5.0, 0.0), 1.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (7.0, 3.0), 504.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (5.0, 2.0), 30.0, TOL1, SpecFunCode::Success);
        test_check_result_and_code!(
            poch_e,
            (5.0, 1.0 / 256.0),
            1.005_902_310_615_136_6,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            poch_e,
            (-9.0, -4.0),
            1.0 / 17160.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            poch_e,
            (-9.0, -3.0),
            -1.0 / 1320.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(poch_e, (-9.0, -3.5), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-9.0, 4.0), 3024.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-9.0, 3.0), -504.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-9.0, 3.5), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-9.0, 0.0), 1.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            poch_e,
            (-8.0, -4.0),
            1.0 / 11880.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            poch_e,
            (-8.0, -3.0),
            -1.0 / 990.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(poch_e, (-8.0, 4.0), 1680.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-8.0, 3.0), -336.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-3.0, 4.0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-3.0, 3.0), -6.0, TOL2, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-4.0, 4.0), 24.0, TOL2, SpecFunCode::Success);
        test_check_result_and_code!(poch_e, (-3.0, 100.0), 0.0, TOL0, SpecFunCode::Success);
    }

    #[test]
    fn test_pochrel_e() {
        test_check_result_and_code!(
            pochrel_e,
            (5.0, 0.0),
            1.506_117_668_431_800_5,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (7.0, 3.0),
            503.0 / 3.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (5.0, 2.0),
            29.0 / 2.0,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (5.0, 0.01),
            1.518_639_366_136_827_6,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (-5.5, 0.01),
            1.858_494_563_382_906_4,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (-5.5, -1.0 / 8.0),
            1.088_331_930_355_213_5,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (-5.5, -1.0 / 256.0),
            1.767_826_803_772_617_7,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (-5.5, -11.0),
            0.090_909_090_909_396_53,
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_taylorcoeff_e() {
        test_check_result_and_code!(
            taylorcoeff_e,
            (10, 1.0 / 1_048_576.0),
            1.714_896_185_477_607_2e-67,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            taylorcoeff_e,
            (10, 1.0 / 1024.0),
            2.173_889_178_849_79e-37,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            taylorcoeff_e,
            (10, 1.0),
            2.755_731_922_398_589e-7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            taylorcoeff_e,
            (10, 5.0),
            2.691_144_455_467_372_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            taylorcoeff_e,
            (10, 500.0),
            2.691_144_455_467_372e20,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            taylorcoeff_e,
            (100, 100.0),
            1.071_510_288_125_466_9e42,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            taylorcoeff_e,
            (1000, 200.0),
            2.662_879_055_815_474_7e-_267,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            taylorcoeff_e,
            (1000, 500.0),
            2.319_317_013_974_085_6e131,
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_fact_e() {
        test_check_result_and_code!(fact_e, (0), 1.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(fact_e, (1), 1.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(fact_e, (7), 5040.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            fact_e,
            (33),
            8.683_317_618_811_886e36,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_doublefact_e() {
        test_check_result_and_code!(doublefact_e, (0), 1.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(doublefact_e, (1), 1.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(doublefact_e, (7), 105.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            doublefact_e,
            (33),
            6.332_659_870_762_85e18,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_lndoublefact_e() {
        test_check_result_and_code!(lndoublefact_e, (0), 0.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            lndoublefact_e,
            (7),
            4.653_960_350_157_523,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lndoublefact_e,
            (33),
            43.292_252_022_541_72,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lndoublefact_e,
            (34),
            45.288_575_519_655_96,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lndoublefact_e,
            (1034),
            3_075.638_379_627_12,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lndoublefact_e,
            (1035),
            3_078.883_908_173_180_8,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_choose_e() {
        test_check_result_and_code!(choose_e, (7, 3), 35.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(choose_e, (7, 4), 35.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(choose_e, (5, 2), 10.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(choose_e, (5, 3), 10.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            choose_e,
            (500, 495),
            255_244_687_600.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            choose_e,
            (500, 5),
            255_244_687_600.0,
            TOL0,
            SpecFunCode::Success
        );
        //TODO: This test fails
        test_check_result_and_code!(
            choose_e,
            (500, 200),
            5.054_949_849_935_532_4e144,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            choose_e,
            (500, 300),
            5.054_949_849_935_532_4e144,
            TOL5,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_lnchoose_e() {
        test_check_result_and_code!(
            lnchoose_e,
            (7, 3),
            3.555_348_061_489_413_5,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnchoose_e,
            (5, 2),
            std::f64::consts::LN_10,
            TOL0,
            SpecFunCode::Success
        );
    }
}
