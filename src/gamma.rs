mod data;
pub(crate) mod mono;
pub(crate) mod poch;
pub(crate) mod psi;

use crate::result::SpecFunResult;

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
    /// use stercus::specfun::gamma::DiGamma;
    ///
    /// let res = 2.0_f64.digamma_e();
    /// let abserr = (res.val + 0.5772156649).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn digamma_e(&self) -> SpecFunResult<f64>;
    /// Compute the digamma function of a number.
    ///
    /// ```
    /// use stercus::specfun::gamma::DiGamma;
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
    /// use stercus::specfun::gamma::TriGamma;
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
    /// use stercus::specfun::gamma::TriGamma;
    ///
    /// let res = 2.0_f64.trigamma();
    /// let abserr = (res - 1.644934067).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn trigamma(&self) -> f64;
    /// Compute the polygamma function of a number along with the error.
    ///
    /// ```
    /// use stercus::specfun::gamma::PolyGamma;
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
    /// use stercus::specfun::gamma::PolyGamma;
    ///
    /// let res = 2.0_f64.polygamma_e(5);
    /// let abserr = (res.val - 122.0811674).abs();
    /// assert!(abserr < 1e-10);
    /// ```
    fn polygamma(&self, n: usize) -> f64;
}

macro_rules! impl_gamma_int {
    ($T:ty) => {
        impl Gamma for $T {
            fn gamma_e(&self) -> SpecFunResult<f64> {
                factorial(*self as usize)
            }
            fn gamma(&self) -> f64 {
                factorial(*self as usize).val
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
    use crate::consts::SQRT_DLB_EPS;
    #[macro_use]
    use crate::test_utils::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_check_result_sgn_and_code;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131072.0 * f64::EPSILON;

    #[test]
    fn test_ln_gamma_e() {
        test_check_result_and_code!(
            lngamma_e,
            (-0.1),
            2.368961332728788655,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-1.0 / 256.0),
            5.547444766967471595,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (1.0e-08),
            18.420680738180208905,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (0.1),
            2.252712651734205,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (1.0 + 1.0 / 256.0),
            -0.0022422226599611501448,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (2.0 + 1.0 / 256.0),
            0.0016564177556961728692,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (100.0),
            359.1342053695753,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-1.0 - 1.0 / 65536.0),
            11.090348438090047844,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-1.0 - 1.0 / 268435456.0),
            19.408121054103474300,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-100.5),
            -364.9009683094273518,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lngamma_e,
            (-100.0 - 1.0 / 65536.0),
            -352.6490910117097874,
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
            8.525161361065414300,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnfact_e,
            (33),
            85.05446701758151741,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_gamma_e() {
        test_check_result_and_code!(
            gamma_e,
            (1.0 + 1.0 / 4096.0),
            0.9998591371459403421,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (1.0 + 1.0 / 32.0),
            0.9829010992836269148,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (2.0 + 1.0 / 256.0),
            1.0016577903733583299,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(gamma_e, (9.0), 40320.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(gamma_e, (10.0), 362880.0, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(
            gamma_e,
            (100.0),
            9.332621544394415268e+155,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (170.0),
            4.269068009004705275e+304,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (171.0),
            7.257415615307998967e+306,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (-10.5),
            -2.640121820547716316e-07,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (-11.25),
            6.027393816261931672e-08,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gamma_e,
            (-1.0 + 1.0 / 65536.0),
            -65536.42280587818970,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_gammastar_e() {
        test_check_result_and_code!(
            gammastar_e,
            (1.0e-08),
            3989.423555759890865,
            2.0 * TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (1.0e-05),
            126.17168469882690233,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (0.001),
            12.708492464364073506,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (1.5),
            1.0563442442685598666,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (3.0),
            1.0280645179187893045,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (9.0),
            1.0092984264218189715,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (11.0),
            1.0076024283104962850,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (100.0),
            1.0008336778720121418,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammastar_e,
            (1.0e+05),
            1.0000008333336805529,
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
            1.0 / 362880.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammainv_e,
            (100.0),
            1.0715102881254669232e-156,
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
            -1.0 / 2.640121820547716316e-07,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammainv_e,
            (-11.25),
            1.0 / 6.027393816261931672e-08,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            gammainv_e,
            (-1.0 + 1.0 / 65536.0),
            -1.0 / 65536.42280587818970,
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
            0.000022981557571259389129,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnpoch_e,
            (5.0, 1.0 / 256.0),
            0.005884960217985189004,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnpoch_e,
            (7.0, 3.0),
            6.222576268071368616,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            lnpoch_e,
            (5.0, 2.0),
            3.401197381662155375,
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
            0.7430116475119920117,
            TOL1,
            1.0,
            SpecFunCode::Success
        );
        test_check_result_sgn_and_code!(
            lnpoch_sgn_e,
            (-4.5, 1.25),
            2.1899306304483174731,
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
            1.0059023106151364982,
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
            1.506117668431800472,
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
            1.5186393661368275330,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (-5.5, 0.01),
            1.8584945633829063516,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (-5.5, -1.0 / 8.0),
            1.0883319303552135488,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (-5.5, -1.0 / 256.0),
            1.7678268037726177453,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            pochrel_e,
            (-5.5, -11.0),
            0.09090909090939652475,
            TOL1,
            SpecFunCode::Success
        );
    }
}
