mod airy_core;
mod airy_data;
mod airy_deriv;
mod airy_zeros;
use crate::result::SpecFunResult;

use airy_core::*;
use airy_deriv::*;
use airy_zeros::*;

pub trait Airy {
    /// Compute the Airy function Ai(x) along with an error estimate.
    fn airy_ai_e(&self) -> SpecFunResult<f64>;
    /// Compute the Airy function Ai(x) along with an error estimate,
    /// scaled by S(x) = exp(2x^{3/2}/3) for x > 0 and 1 for x < 0.
    fn airy_ai_scaled_e(&self) -> SpecFunResult<f64>;
    /// Compute the Airy function Ai(x).
    fn airy_ai(&self) -> f64;
    /// Compute the Airy function Ai(x) scaled by S(x) = exp(2x^{3/2}/3)
    /// for x > 0 and 1 for x < 0.
    fn airy_ai_scaled(&self) -> f64;
    /// Compute the Airy function Bi(x) along with an error estimate.
    fn airy_bi_e(&self) -> SpecFunResult<f64>;
    /// Compute the Airy function Bi(x) along with an error estimate,
    /// scaled by S(x) = exp(-2x^{3/2}/3) for x > 0 and 1 for x < 0.
    fn airy_bi_scaled_e(&self) -> SpecFunResult<f64>;
    /// Compute the Airy function Bi(x).
    fn airy_bi(&self) -> f64;
    /// Compute the Airy function Bi(x) scaled by S(x) = exp(-2x^{3/2}/3)
    /// for x > 0 and 1 for x < 0.
    fn airy_bi_scaled(&self) -> f64;
    /// Compute first derivative of the Airy function Ai(x) along with an
    /// error estimate.
    fn airy_ai_deriv_e(&self) -> SpecFunResult<f64>;
    /// Compute first derivative of the Airy function Ai(x).
    fn airy_ai_deriv(&self) -> f64;
    /// Compute the first derivative of the Airy function Ai(x) along with an
    /// error estimate, scaled by S(x) = exp(2x^{3/2}/3) for x > 0 and 1 for
    /// x < 0.
    fn airy_ai_deriv_scaled_e(&self) -> SpecFunResult<f64>;
    /// Compute the first derivative of the Airy function Ai(x) scaled by
    /// S(x) = exp(2x^{3/2}/3) for x > 0 and 1 for x < 0.
    fn airy_ai_deriv_scaled(&self) -> f64;
    /// Compute first derivative of the Airy function Bi(x) along with an
    /// error estimate.
    fn airy_bi_deriv_e(&self) -> SpecFunResult<f64>;
    /// Compute first derivative of the Airy function Bi(x).
    fn airy_bi_deriv(&self) -> f64;
    /// Compute the first derivative of the Airy function Bi(x) along with an
    /// error estimate, scaled by S(x) = exp(-2x^{3/2}/3) for x > 0 and 1 for
    /// x < 0.
    fn airy_bi_deriv_scaled_e(&self) -> SpecFunResult<f64>;
    /// Compute the first derivative of the Airy function Bi(x) scaled by
    /// S(x) = exp(-2x^{3/2}/3) for x > 0 and 1 for x < 0.
    fn airy_bi_deriv_scaled(&self) -> f64;
}

macro_rules! impl_airy_float {
    ($T:ty) => {
        impl Airy for $T {
            fn airy_ai_e(&self) -> SpecFunResult<f64> {
                airy_ai_e(*self as f64)
            }
            fn airy_ai_scaled_e(&self) -> SpecFunResult<f64> {
                airy_ai_scaled_e(*self as f64)
            }
            fn airy_ai(&self) -> f64 {
                airy_ai_e(*self as f64).val
            }
            fn airy_ai_scaled(&self) -> f64 {
                airy_ai_scaled_e(*self as f64).val
            }
            fn airy_bi_e(&self) -> SpecFunResult<f64> {
                airy_bi_e(*self as f64)
            }
            fn airy_bi_scaled_e(&self) -> SpecFunResult<f64> {
                airy_bi_scaled_e(*self as f64)
            }
            fn airy_bi(&self) -> f64 {
                airy_bi_e(*self as f64).val
            }
            fn airy_bi_scaled(&self) -> f64 {
                airy_bi_scaled_e(*self as f64).val
            }

            fn airy_ai_deriv_e(&self) -> SpecFunResult<f64> {
                airy_ai_deriv_e(*self as f64)
            }
            fn airy_ai_deriv(&self) -> f64 {
                airy_ai_deriv_e(*self as f64).val
            }
            fn airy_ai_deriv_scaled_e(&self) -> SpecFunResult<f64> {
                airy_ai_deriv_scaled_e(*self as f64)
            }
            fn airy_ai_deriv_scaled(&self) -> f64 {
                airy_ai_deriv_scaled_e(*self as f64).val
            }

            fn airy_bi_deriv_e(&self) -> SpecFunResult<f64> {
                airy_bi_deriv_e(*self as f64)
            }
            fn airy_bi_deriv(&self) -> f64 {
                airy_bi_deriv_e(*self as f64).val
            }
            fn airy_bi_deriv_scaled_e(&self) -> SpecFunResult<f64> {
                airy_bi_deriv_scaled_e(*self as f64)
            }
            fn airy_bi_deriv_scaled(&self) -> f64 {
                airy_bi_deriv_scaled_e(*self as f64).val
            }
        }
    };
}

impl_airy_float!(f32);
impl_airy_float!(f64);

trait AiryZero {
    /// Compute a zero of the Airy function Ai(x) along with an error estimate.
    fn airy_zero_ai_e(&self) -> SpecFunResult<f64>;
    /// Compute a zero of the Airy function Ai(x).
    fn airy_zero_ai(&self) -> f64;
    /// Compute a zero of the first derivative of the Airy function Ai(x) along
    /// with an error estimate.
    fn airy_zero_ai_deriv_e(&self) -> SpecFunResult<f64>;
    /// Compute a zero of the first derivative of the Airy function Ai(x).
    fn airy_zero_ai_deriv(&self) -> f64;
    /// Compute a zero of the Airy function Bi(x) along with an error estimate.
    fn airy_zero_bi_e(&self) -> SpecFunResult<f64>;
    /// Compute a zero of the Airy function Bi(x).
    fn airy_zero_bi(&self) -> f64;
    /// Compute a zero of the first derivative of the Airy function Bi(x) along
    /// with an error estimate.
    fn airy_zero_bi_deriv_e(&self) -> SpecFunResult<f64>;
    /// Compute a zero of the first derivative of the Airy function Bi(x).
    fn airy_zero_bi_deriv(&self) -> f64;
}

macro_rules! impl_airy_zero_int {
    ($T:ty) => {
        impl AiryZero for $T {
            fn airy_zero_ai_e(&self) -> SpecFunResult<f64> {
                airy_zero_ai_e(*self as usize)
            }
            fn airy_zero_ai(&self) -> f64 {
                airy_zero_ai_e(*self as usize).val
            }
            fn airy_zero_ai_deriv_e(&self) -> SpecFunResult<f64> {
                airy_zero_ai_deriv_e(*self as usize)
            }
            fn airy_zero_ai_deriv(&self) -> f64 {
                airy_zero_ai_deriv_e(*self as usize).val
            }

            fn airy_zero_bi_e(&self) -> SpecFunResult<f64> {
                airy_zero_bi_e(*self as usize)
            }
            fn airy_zero_bi(&self) -> f64 {
                airy_zero_bi_e(*self as usize).val
            }
            fn airy_zero_bi_deriv_e(&self) -> SpecFunResult<f64> {
                airy_zero_bi_deriv_e(*self as usize)
            }
            fn airy_zero_bi_deriv(&self) -> f64 {
                airy_zero_bi_deriv_e(*self as usize).val
            }
        }
    };
}

impl_airy_zero_int!(u8);
impl_airy_zero_int!(u16);
impl_airy_zero_int!(u32);
impl_airy_zero_int!(u64);
impl_airy_zero_int!(usize);

#[cfg(test)]
mod test {
    use super::airy_deriv::*;
    use super::airy_zeros::*;
    use super::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_utils::*;
    const TOL0: f64 = 2.0 * f64::EPSILON;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;

    #[test]
    fn test_airy_ai_e() {
        test_check_result_and_code!(
            airy_ai_e,
            (-500.0),
            0.072_590_120_104_041_14,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (-5.0),
            0.350_761_009_024_114_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (-0.300_000_000_000_009_4),
            0.430_903_095_285_583_1,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (0.699_999_999_999_990_7),
            0.189_162_400_398_151_9,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            airy_ai_e,
            (1.649_999_999_999_991),
            0.058_310_586_187_208_85,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            airy_ai_e,
            (2.549_999_999_999_99),
            0.014_461_495_132_954_28,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (3.499_999_999_999_987),
            0.002_584_098_786_989_702,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (5.399_999_999_999_98),
            4.272_986_169_411_866e-5,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_ai_scaled_e() {
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (-5.0),
            0.350_761_009_024_114_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (0.699_999_999_999_990_7),
            0.279_512_566_768_121_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (1.649_999_999_999_991),
            0.239_549_300_144_274_1,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (2.549_999_999_999_99),
            0.218_365_859_589_938_8,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (3.499_999_999_999_987),
            0.203_292_080_816_351_9,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (5.399_999_999_999_98),
            0.183_605_009_328_222_9,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_airy_bi_e() {
        test_check_result_and_code!(
            airy_bi_e,
            (-500.0),
            -0.094_688_570_132_991_03,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (-5.0),
            -0.138_369_134_901_600_5,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (0.699_999_999_999_990_7),
            0.973_328_655_878_159_9,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (1.649_999_999_999_991),
            2.196_407_956_850_028,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (2.549_999_999_999_99),
            6.973_628_612_493_443,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (3.499_999_999_999_987),
            33.055_506_754_610_69,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (5.399_999_999_999_98),
            1_604.476_078_241_272,
            TOL1,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_airy_bi_scaled_e() {
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (-5.0),
            -0.138_369_134_901_600_5,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (0.699_999_999_999_990_7),
            0.658_708_075_458_230_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (1.649_999_999_999_991),
            0.534_644_999_559_753_9,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (2.549_999_999_999_99),
            0.461_835_455_542_297,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (3.499_999_999_999_987),
            0.420_177_188_235_306_1,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (5.399_999_999_999_98),
            0.373_405_067_572_047_3,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_ai_deriv_e() {
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (-5.0),
            0.327_192_818_554_443_5,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (-0.550_000_000_000_009_4),
            -0.191_460_498_714_362_9,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (0.499_999_999_999_990_6),
            -0.224_910_532_664_685,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (1.899_999_999_999_992),
            -0.060_436_781_785_757_18,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (3.249_999_999_999_988),
            -0.007_792_687_926_790_889,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (5.199_999_999_999_981),
            -0.000_158_943_452_645_954_3,
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_ai_deriv_scaled_e() {
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (-5.0),
            0.327_192_818_554_443_5,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (0.549_999_999_999_990_6),
            -0.287_405_727_917_016_6,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (1.499_999_999_999_991),
            -0.331_419_979_686_363_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (2.499_999_999_999_99),
            -0.366_108_938_475_162,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (3.649_999_999_999_986),
            -0.397_403_383_145_396_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (6.299_999_999_999_977),
            -0.450_879_918_958_594_7,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_bi_deriv_e() {
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (-5.0),
            0.778_411_773_001_899,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (-0.550_000_000_000_009_4),
            0.515_578_535_876_501_4,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (0.499_999_999_999_990_6),
            0.544_572_564_140_588_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (1.899_999_999_999_992),
            3.495_165_862_891_568,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (3.249_999_999_999_988),
            36.554_851_492_503_38,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (5.199_999_999_999_981),
            2_279.748_293_583_233,
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_bi_deriv_scaled_e() {
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (-5.0),
            0.778_411_773_001_899,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (0.549_999_999_999_990_6),
            0.432_281_128_181_756_6,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (1.499_999_999_999_991),
            0.554_230_756_391_803_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (2.499_999_999_999_99),
            0.675_538_444_164_498_5,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (3.649_999_999_999_986),
            0.761_395_937_300_022_8,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (6.299_999_999_999_977),
            0.885_206_413_973_757_1,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_zero_ai_e() {
        test_check_result_and_code!(
            airy_zero_ai_e,
            (2),
            -4.087_949_444_130_97,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_e,
            (50),
            -38.021_008_677_255_25,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_e,
            (100),
            -60.455_557_274_116_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_e,
            (110),
            -64.431_356_709_913_25,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_zero_bi_e() {
        test_check_result_and_code!(
            airy_zero_bi_e,
            (2),
            -3.271_093_302_836_353,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (50),
            -37.765_834_381_651_8,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (100),
            -60.253_364_825_808_37,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (110),
            -64.235_516_760_656_15,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (111),
            -64.626_899_481_951_94,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (200),
            -95.886_991_473_566_82,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_zero_ai_deriv_e() {
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (2),
            -3.248_197_582_179_836_6,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (50),
            -37.765_659_100_538_87,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (100),
            -60.253_295_964_424_794,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (110),
            -64.235_456_172_435_46,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (1000),
            -280.937_808_035_893_5,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_zero_bi_deriv_e() {
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (2),
            -4.073_155_089_071_828,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (50),
            -38.020_835_740_957_885,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (100),
            -60.455_488_872_571_41,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (110),
            -64.431_296_489_448_46,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (111),
            -64.822_087_375_842_06,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (200),
            -96.047_310_503_103_25,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (1000),
            -281.031_516_447_111_87,
            TOL0,
            SpecFunCode::Success
        );
    }
}
