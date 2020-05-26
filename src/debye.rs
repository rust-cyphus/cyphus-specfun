use crate::cheb::ChebSeries;
use crate::result::{SpecFunCode, SpecFunResult};
use lazy_static::lazy_static;

lazy_static! {
    static ref ADEB1_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.400_659_719_038_141,
            0.193_721_304_218_936_02,
            -6.232_912_455_489_577e-3,
            3.511_174_770_206_48e-4,
            -0.228_222_466_701_231e-04,
            0.158_054_678_750_30e-05,
            -0.113_537_819_707_2e-06,
            0.835_833_611_88e-08,
            -0.626_442_478_7e-09,
            0.476_033_489e-10,
            -0.365_741_54e-11,
            0.283_543_1e-12,
            -0.221_473e-13,
            0.17409e-14,
            -0.1376e-15,
            0.109e-16,
            -0.9e-18,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB2_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.594_381_023_257_077,
            0.286_335_720_453_072,
            -1.020_626_561_580_467_2e-2,
            6.049_109_775_346_844e-4,
            -0.405_257_658_950_210e-04,
            0.286_338_263_288_11e-05,
            -0.208_639_430_306_5e-06,
            0.155_237_875_826e-07,
            -0.117_312_800_87e-08,
            0.897_358_589e-10,
            -0.693_176_14e-11,
            0.539_805_7e-12,
            -0.423_241e-13,
            0.33378e-14,
            -0.2645e-15,
            0.211e-16,
            -0.17e-17,
            0.1e-18,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB3_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.707_737_068_327_441,
            0.340_068_135_211_091_75,
            -1.294_515_018_444_087e-2,
            7.963_755_380_173_816e-4,
            -0.546_360_009_590_824e-04,
            0.392_430_195_988_05e-05,
            -0.289_403_282_353_9e-06,
            0.217_317_613_962e-07,
            -0.165_420_999_50e-08,
            0.127_279_618_9e-09,
            -0.987_963_460e-11,
            0.772_507_4e-12,
            -0.607_797e-13,
            0.48076e-14,
            -0.3820e-15,
            0.305e-16,
            -0.24e-17,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB4_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.781_869_415_020_523_7,
            0.374_976_783_526_892_84,
            -1.494_090_739_903_158_4e-2,
            0.945_679_811_437_042e-03,
            -0.661_329_161_389_33e-04,
            0.481_563_298_214_4e-05,
            -0.358_808_395_875_9e-06,
            0.271_601_187_416e-07,
            -0.208_070_991_22e-08,
            0.160_938_386_9e-09,
            -0.125_470_979e-10,
            0.984_726_5e-12,
            -0.777_237e-13,
            0.61648e-14,
            -0.4911e-15,
            0.393e-16,
            -0.32e-17,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB5_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.834_026_954_683_453,
            0.399_409_885_710_626_63,
            -1.645_667_647_730_996_6e-2,
            1.065_213_834_066_454e-3,
            -0.756_730_374_875_418e-4,
            0.557_459_852_402_73e-5,
            -0.419_069_233_091_8e-6,
            0.319_456_143_678e-7,
            -0.246_133_181_71e-8,
            0.191_280_163_3e-9,
            -0.149_720_049e-10,
            0.117_903_12e-11,
            -0.933_329e-13,
            0.74218e-14,
            -0.5925e-15,
            0.475e-16,
            -0.39e-17,
        ],
        a: -1.0,
        b: 1.0,
    };
    static ref ADEB6_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            2.872_672_713_413_012_3,
            0.417_437_535_233_902_76,
            -1.764_538_493_540_678_7e-2,
            1.162_985_273_349_455_5e-3,
            -0.837_118_027_357_117e-4,
            0.622_836_115_961_89e-5,
            -0.471_864_446_563_6e-6,
            0.361_950_397_806e-7,
            -0.280_303_680_10e-8,
            0.218_768_198_3e-9,
            -0.171_857_387e-10,
            0.135_758_09e-11,
            -0.107_758_0e-12,
            0.85893e-14,
            -0.6872e-15,
            0.552e-16,
            -0.44e-17,
        ],
        a: -1.0,
        b: 1.0,
    };
}

pub fn debye1_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 1.644_934_066_848_226_4;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 0.25 * x + x * x / 36.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB1_CHEB).eval(t);
        result.val = c.val - 0.25 * x;
        result.err = c.err + 0.25 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for _i in (1..(nexp + 1)).rev() {
            sum *= ex;
            sum += 1f64.mul_add(xk.recip(), 1.0) / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / x - sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        result.val = (val_infinity - (-x).exp() * (x + 1.0)) / x;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = val_infinity / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye2_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 4.808_227_612_638_377;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - x / 3.0 + x * x / 24.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB2_CHEB).eval(t);
        result.val = c.val - x / 3.0;
        result.err = c.err + x / 3.0 * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for _i in (1..(nexp + 1)).rev() {
            let xk_inv = xk.recip();
            sum *= ex;
            sum += 2f64.mul_add(xk_inv, 2.0).mul_add(xk_inv, 1.0) / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x) - 2.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x2 = x * x;
        let sum = 2.0 + 2.0 * x + x2;
        result.val = (val_infinity - 2.0 * sum * (-x).exp()) / x2;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = val_infinity / x / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye3_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 19.481_818_206_800_487;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 3.0 * x / 8.0 + x * x / 20.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB3_CHEB).eval(t);
        result.val = c.val - 0.375 * x;
        result.err = c.err + 0.375 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for _i in (1..(nexp + 1)).rev() {
            let xk_inv = 1.0 / xk;
            sum *= ex;
            sum += 6f64
                .mul_add(xk_inv, 6.0)
                .mul_add(xk_inv, 3.0)
                .mul_add(xk_inv, 1.0)
                / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x * x) - 3.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x3 = x * x * x;
        let sum = 6.0 + 6.0 * x + 3.0 * x * x + x3;
        result.val = (val_infinity - 3.0 * sum * (-x).exp()) / x3;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = val_infinity / x / x / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye4_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 99.545_064_493_763_52;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * std::f64::consts::SQRT_2 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 2.0 * x / 5.0 + x * x / 18.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB4_CHEB).eval(t);
        result.val = c.val - 2.0 / 5.0 * x;
        result.err = c.err + 2.0 / 5.0 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for _i in (1..(nexp + 1)).rev() {
            let xk_inv = 1.0 / xk;
            sum *= ex;
            sum += 24f64
                .mul_add(xk_inv, 24.0)
                .mul_add(xk_inv, 12.0)
                .mul_add(xk_inv, 4.0)
                .mul_add(xk_inv, 1.0)
                / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x * x * x) - 4.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x2 = x * x;
        let x4 = x2 * x2;
        let sum = 24.0 + 24.0 * x + 12.0 * x2 + 4.0 * x2 * x + x4;
        result.val = (val_infinity - 4.0 * sum * (-x).exp()) / x4;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = val_infinity / x / x / x / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye5_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 610.405_837_190_669_4;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 5.0 * x / 12.0 + 5.0 * x * x / 84.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB5_CHEB).eval(t);
        result.val = c.val - 5.0 / 12.0 * x;
        result.err = c.err + 5.0 / 12.0 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for _i in (1..(nexp + 1)).rev() {
            let xk_inv = 1.0 / xk;
            sum *= ex;
            sum += 120f64
                .mul_add(xk_inv, 120.0)
                .mul_add(xk_inv, 60.0)
                .mul_add(xk_inv, 20.0)
                .mul_add(xk_inv, 5.0)
                .mul_add(xk_inv, 1.0)
                / rk;

            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x * x * x * x) - 5.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x2 = x * x;
        let x4 = x2 * x2;
        let x5 = x4 * x;
        let sum = 120.0 + 120.0 * x + 60.0 * x2 + 20.0 * x2 * x + 5.0 * x4 + x5;
        result.val = (val_infinity - 5.0 * sum * (-x).exp()) / x5;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = ((((val_infinity / x) / x) / x) / x) / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

pub fn debye6_e(x: f64) -> SpecFunResult<f64> {
    let val_infinity = 4_356.068_878_289_907;
    let xcut = -crate::consts::LN_DBL_MIN;

    let mut result = SpecFunResult::<f64>::default();

    if x < 0.0 {
        result.code = SpecFunCode::DomainErr;
        result.val = f64::NAN;
        result.err = f64::NAN;
    } else if x < 2.0 * crate::consts::SQRT_DLB_EPS {
        result.val = 1.0 - 3.0 * x / 7.0 + x * x / 16.0;
        result.err = f64::EPSILON * result.val.abs();
    } else if x <= 4.0 {
        let t = x * x / 8.0 - 1.0;
        let c = (*ADEB6_CHEB).eval(t);
        result.val = c.val - 3.0 / 7.0 * x;
        result.err = c.err + 3.0 / 7.0 * x * f64::EPSILON;
    } else if x < -(std::f64::consts::LN_2 + crate::consts::LN_DBL_EPS) {
        let nexp = (xcut / x).floor() as usize;
        let ex = (-x).exp();
        let mut sum = 0.0;
        let mut xk = x * nexp as f64;
        let mut rk = nexp as f64;
        for _i in (1..(nexp + 1)).rev() {
            let xk_inv = 1.0 / xk;
            sum *= ex;
            sum += 720f64
                .mul_add(xk_inv, 720.0)
                .mul_add(xk_inv, 360.0)
                .mul_add(xk_inv, 120.0)
                .mul_add(xk_inv, 30.0)
                .mul_add(xk_inv, 6.0)
                .mul_add(xk_inv, 1.0)
                / rk;
            rk -= 1.0;
            xk -= x;
        }
        result.val = val_infinity / (x * x * x * x * x * x) - 6.0 * sum * ex;
        result.err = f64::EPSILON * result.val.abs();
    } else if x < xcut {
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x4 * x2;
        let sum = 720.0 + 720.0 * x + 360.0 * x2 + 120.0 * x2 * x + 30.0 * x4 + 6.0 * x4 * x + x6;
        result.val = (val_infinity - 6.0 * sum * (-x).exp()) / x6;
        result.err = f64::EPSILON * result.val.abs();
    } else {
        result.val = (((((val_infinity / x) / x) / x) / x) / x) / x;
        result.err = f64::EPSILON * result.val.abs();
    }
    result
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_utils::*;

    const TOL0: f64 = 2.0 * f64::EPSILON;

    #[test]
    fn test_debye1_e() {
        test_check_result_and_code!(
            debye1_e,
            (0.1),
            0.975_277_750_004_723_3,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye1_e,
            (1.0),
            0.777_504_634_112_248_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye1_e,
            (10.0),
            0.164_443_465_679_946_03,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_debye2_e() {
        test_check_result_and_code!(
            debye2_e,
            (0.1),
            0.967_083_287_045_302_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye2_e,
            (1.0),
            0.707_878_475_627_829_2,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye2_e,
            (10.0),
            0.047_971_498_020_121_87,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_debye3_e() {
        test_check_result_and_code!(
            debye3_e,
            (0.1),
            0.962_999_940_487_211,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye3_e,
            (1.0),
            0.674_415_564_077_814_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye3_e,
            (10.0),
            0.019_295_765_690_345_49,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_debye4_e() {
        test_check_result_and_code!(
            debye4_e,
            (0.1),
            0.960_555_486_124_335_9,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye4_e,
            (1.0),
            0.654_874_068_886_737,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye4_e,
            (10.0),
            0.009_673_675_560_271_159,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_debye5_e() {
        test_check_result_and_code!(
            debye5_e,
            (0.1),
            0.958_928_494_283_105_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye5_e,
            (1.0),
            0.642_100_258_021_779_1,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye5_e,
            (10.0),
            0.005_701_535_852_992_909,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_debye6_e() {
        test_check_result_and_code!(
            debye6_e,
            (0.1),
            0.957_767_773_826_054_7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye6_e,
            (1.0),
            0.633_111_425_834_951_1,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            debye6_e,
            (10.0),
            3.793_849_329_461_595_3e-3,
            TOL0,
            SpecFunCode::Success
        );
    }
}
