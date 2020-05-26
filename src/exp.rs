pub(crate) mod core;

#[cfg(test)]
mod test {
    use super::core::*;
    use crate::consts::{LN_DBL_MAX, SQRT_DLB_EPS};
    use crate::test_utils::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_check_result_and_code_e10;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131_072.0 * f64::EPSILON;

    #[test]
    fn test_exp_e() {
        test_check_result_and_code!(
            exp_e,
            (-10.0),
            (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(exp_e, (10.0), (10.0_f64).exp(), TOL0, SpecFunCode::Success);
    }

    #[test]
    fn test_exp_e10_e() {
        test_check_result_and_code_e10!(
            exp_e10_e,
            (1.0),
            std::f64::consts::E,
            0,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (2000.0),
            3.881_180_194_283_637_3,
            868,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (100.0),
            2.688_117_141_816_135_6e43,
            0,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (1000.0),
            1.970_071_114_017_047,
            434,
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (-100.0),
            3.720_075_976_020_836e-44,
            0,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_e10_e,
            (-1000.0),
            5.075_958_897_549_457,
            -435,
            TOL3,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_err_e() {
        test_check_result_and_code!(
            exp_err_e,
            (-10.0, TOL1),
            (-10.0_f64).exp(),
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_err_e,
            (10.0, TOL1),
            (10.0_f64).exp(),
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_err_e10_e() {
        test_check_result_and_code_e10!(
            exp_err_e10_e,
            (1.0, SQRT_TOL0),
            std::f64::consts::E,
            0,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_err_e10_e,
            (2000.0, 1e-10),
            3.881_180_194_283_637_3,
            868,
            1e-7,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_mult_e() {
        let x = 0.8 * LN_DBL_MAX;
        test_check_result_and_code!(
            exp_mult_e,
            (-10.0, 1.0e-06),
            1.0e-06 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (-10.0, 2.0),
            2.0 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (-10.0, -2.0),
            -2.0 * (-10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (10.0, 1.0e-06),
            1.0e-06 * (10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (10.0, -2.0),
            -2.0 * (10.0_f64).exp(),
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 1.00001),
            1.00001 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 1.000_001),
            1.000_001 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 1.000_000_001),
            1.000_000_001 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 100.0),
            100.0 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, 1.0e+20),
            1.0e+20 * x.exp(),
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_e,
            (x, (-x).exp() * std::f64::consts::LN_2.exp()),
            2.0,
            TOL4,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_mult_err_e() {
        let x = 0.8 * LN_DBL_MAX;
        test_check_result_and_code!(
            exp_mult_err_e,
            (-10.0, SQRT_TOL0, 2.0, SQRT_TOL0),
            2.0 * (-10.0_f64).exp(),
            SQRT_TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exp_mult_err_e,
            (
                x,
                SQRT_TOL0 * x,
                (-x).exp() * std::f64::consts::LN_2.exp(),
                SQRT_TOL0 * (-x).exp() * std::f64::consts::LN_2.exp()
            ),
            2.0,
            SQRT_TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_mult_e10_e() {
        test_check_result_and_code_e10!(
            exp_mult_e10_e,
            (1.0, 1.0),
            std::f64::consts::E,
            0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_mult_e10_e,
            (1000.0, 1e200),
            1.970_071_114_017_047,
            634,
            TOL3,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_mult_e10_e,
            (10000.0, 1.0),
            8.806_818_225_662_921,
            4342,
            TOL5,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_mult_e10_e,
            (100.0, 1.0),
            2.688_117_141_816_135_6e43,
            0,
            TOL2,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exp_mult_err_e10_e() {
        test_check_result_and_code_e10!(
            exp_mult_err_e10_e,
            (1.0, TOL0, 1.0, TOL0),
            std::f64::consts::E,
            0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code_e10!(
            exp_mult_err_e10_e,
            (1000.0, 1e-12, 1e200, 1e190),
            1.970_071_114_016_566,
            634,
            TOL3,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_expm1_e() {
        test_check_result_and_code!(
            expm1_e,
            (-10.0),
            (-10.0_f64).exp() - 1.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (-0.001),
            -0.000_999_500_166_625_008_5,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (-1.0e-8),
            -1.0e-08 + 0.5e-16,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (1.0e-8),
            1.0e-08 + 0.5e-16,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (0.001),
            0.001_000_500_166_708_341_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            expm1_e,
            (10.0),
            (10.0_f64).exp() - 1.0,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exprel_e() {
        test_check_result_and_code!(
            exprel_e,
            (-10.0),
            0.099_995_460_007_023_75,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (-0.001),
            0.999_500_166_625_008_4,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (-1.0e-8),
            1.0 - 0.5e-08,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (1.0e-8),
            1.0 + 0.5e-08,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (0.001),
            1.000_500_166_708_341_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_e,
            (10.0),
            2_202.546_579_480_671_4,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exprel_2_e() {
        test_check_result_and_code!(
            exprel_2_e,
            (-10.0),
            0.180_000_907_998_595_25,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (-0.001),
            0.999_666_749_983_336_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (-1.0e-8),
            0.999_999_996_666_666_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (1.0e-8),
            1.000_000_003_333_333_4,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (0.001),
            1.000_333_416_683_336_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_2_e,
            (10.0),
            440.309_315_896_134_3,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_exprel_n_e() {
        test_check_result_and_code!(
            exprel_n_e,
            (3, -1000.0),
            0.002_994_006,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -100.0),
            0.029_406,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -10.0),
            0.245_999_727_600_421_44,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -3.0),
            0.544_491_762_584_919_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -0.001),
            0.999_750_049_991_667_8,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, -1.0e-8),
            0.999_999_997_5,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 1.0e-8),
            1.000_000_002_5,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 0.001),
            1.000_250_050_008_334_4,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 3.0),
            2.574_563_760_708_370_5,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 3.1),
            2.677_241_706_846_021,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 10.0),
            131.792_794_768_840_3,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (3, 100.0),
            1.612_870_285_089_681_2e38,
            TOL2,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            exprel_n_e,
            (50, -1000.0),
            0.047_662_316_092_539_76,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -100.0),
            0.334_824_757_234_588_93,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -10.0),
            0.835_628_705_185_328_6,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -3.0),
            0.944_388_160_915_216_4,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -1.0),
            0.980_762_245_565_660_6,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, -1.0e-8),
            1.0 - 1.0e-8 / 51.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 1.0e-8),
            1.0 + 1.0e-8 / 51.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 1.0),
            1.019_992_165_836_668,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 3.0),
            1.062_420_575_746_036_8,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 48.0),
            7.499_573_876_877_195,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 50.1),
            9.311_803_306_230_992,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 100.0),
            8.175_664_432_485_807e7,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (50, 500.0),
            4.806_352_370_663_186e146,
            TOL3,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            exprel_n_e,
            (500, -1000.0),
            0.333_481_580_312_761_9,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -100.0),
            0.833_564_621_753_618_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -10.0),
            0.980_429_780_313_182_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -3.0),
            0.994_047_548_885_067_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -1.0),
            0.998_007_960_238_348_9,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, -1.0e-8),
            1.0 - 1.0e-8 / 501.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 1.0e-8),
            1.0 + 1.0e-8 / 501.0,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 1.0),
            1.001_999_992_016_063_4,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 3.0),
            1.006_024_023_663_244_5,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 48.0),
            1.105_935_551_798_127_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 100.0),
            1.249_222_146_487_828_8,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 500.0),
            28.363_019_877_927_63,
            TOL2,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 1000.0),
            2.403_756_316_033_53e68,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            exprel_n_e,
            (500, 1600.0),
            7.899_293_535_320_607e226,
            TOL4,
            SpecFunCode::Success
        );
    }
}
