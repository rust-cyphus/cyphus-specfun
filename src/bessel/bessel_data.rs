use crate::cheb::ChebSeries;
use lazy_static::lazy_static;

// based on SLATEC besj0, 1977 version, w. fullerton
// Chebyshev expansions for Bessel functions
// series for bj0 on the interval 0.0 to 16.0
// with:
// - weighted error: 7.47e-18
// - log weighted error: 17.13
// - significant figures required: 16.98
// - decimal places required: 17.68
pub(super) const BJ0_DATA: [f64; 13] = [
    0.100_254_161_968_939_13,
    -0.665_223_007_764_405_1,
    0.248_983_703_498_281_3,
    -0.033_252_723_170_035_77,
    0.002_311_417_930_469_401_7,
    -0.000_099_112_774_199_508_0,
    0.000_002_891_670_864_399_8,
    -0.000_000_061_210_858_663_0,
    0.000_000_000_983_865_079_3,
    -0.000_000_000_012_423_551_5,
    0.000_000_000_000_126_543_3,
    -0.000_000_000_000_001_061_9,
    0.000_000_000_000_000_007_4,
];

// based on SLATEC besj1, 1977 version, w. fullerton
// Chebyshev expansions for Bessel functions
// series for bj0 on the interval 0.0 to 16.0
// with:
// - weighted error: 4.48e-17
// - log weighted error: 16.35
// - significant figures required: 15.77
// - decimal places required: 16.89
pub(super) const BJ1_DATA: [f64; 12] = [
    -0.117_261_415_133_327_87,
    -0.253_615_218_307_906_4,
    0.050_127_080_984_469_566,
    -0.004_631_514_809_625_081,
    0.000_247_996_229_415_914,
    -0.000_008_678_948_686_278,
    0.000_000_214_293_917_143,
    -0.000_000_003_936_093_079,
    0.000_000_000_055_911_823,
    -0.000_000_000_000_632_761,
    0.000_000_000_000_005_840,
    -0.000_000_000_000_000_044,
];

// chebyshev expansions for amplitude and phase
// functions used in bessel evaluations
// These are the same for J0,Y0 and for J1,Y1, so
// they sit outside those functions.

pub(super) const BM0_DATA: [f64; 21] = [
    0.092_849_616_373_816_45,
    -0.001_429_877_074_034_84,
    0.000_028_305_792_712_57,
    -0.000_001_433_006_114_24,
    0.000_000_120_286_280_46,
    -0.000_000_013_971_130_13,
    0.000_000_002_040_761_88,
    -0.000_000_000_353_996_69,
    0.000_000_000_070_247_59,
    -0.000_000_000_015_541_07,
    0.000_000_000_003_762_26,
    -0.000_000_000_000_982_82,
    0.000_000_000_000_274_08,
    -0.000_000_000_000_080_91,
    0.000_000_000_000_025_11,
    -0.000_000_000_000_008_14,
    0.000_000_000_000_002_75,
    -0.000_000_000_000_000_96,
    0.000_000_000_000_000_34,
    -0.000_000_000_000_000_12,
    0.000_000_000_000_000_04,
];

pub(super) const BTH0_DATA: [f64; 24] = [
    -0.246_391_637_743_001_19,
    0.001_737_098_307_508_963,
    -0.000_062_183_633_402_968,
    0.000_004_368_050_165_742,
    -0.000_000_456_093_019_869,
    0.000_000_062_197_400_101,
    -0.000_000_010_300_442_889,
    0.000_000_001_979_526_776,
    -0.000_000_000_428_198_396,
    0.000_000_000_102_035_840,
    -0.000_000_000_026_363_898,
    0.000_000_000_007_297_935,
    -0.000_000_000_002_144_188,
    0.000_000_000_000_663_693,
    -0.000_000_000_000_215_126,
    0.000_000_000_000_072_659,
    -0.000_000_000_000_025_465,
    0.000_000_000_000_009_229,
    -0.000_000_000_000_003_448,
    0.000_000_000_000_001_325,
    -0.000_000_000_000_000_522,
    0.000_000_000_000_000_210,
    -0.000_000_000_000_000_087,
    0.000_000_000_000_000_036,
];

pub(super) const BM1_DATA: [f64; 21] = [
    0.104_736_251_093_128_5,
    0.004_424_438_937_023_45,
    -0.000_056_616_395_040_35,
    0.000_002_313_494_173_39,
    -0.000_000_173_771_820_07,
    0.000_000_018_932_099_30,
    -0.000_000_002_654_160_23,
    0.000_000_000_447_402_09,
    -0.000_000_000_086_917_95,
    0.000_000_000_018_914_92,
    -0.000_000_000_004_518_84,
    0.000_000_000_001_167_65,
    -0.000_000_000_000_322_65,
    0.000_000_000_000_094_50,
    -0.000_000_000_000_029_13,
    0.000_000_000_000_009_39,
    -0.000_000_000_000_003_15,
    0.000_000_000_000_001_09,
    -0.000_000_000_000_000_39,
    0.000_000_000_000_000_14,
    -0.000_000_000_000_000_05,
];

pub(super) const BTH1_DATA: [f64; 24] = [
    0.740_601_410_263_138_5,
    -0.004_571_755_659_637_69,
    0.000_119_818_510_964_326,
    -0.000_006_964_561_891_648,
    0.000_000_655_495_621_447,
    -0.000_000_084_066_228_945,
    0.000_000_013_376_886_564,
    -0.000_000_002_499_565_654,
    0.000_000_000_529_495_100,
    -0.000_000_000_124_135_944,
    0.000_000_000_031_656_485,
    -0.000_000_000_008_668_640,
    0.000_000_000_002_523_758,
    -0.000_000_000_000_775_085,
    0.000_000_000_000_249_527,
    -0.000_000_000_000_083_773,
    0.000_000_000_000_029_205,
    -0.000_000_000_000_010_534,
    0.000_000_000_000_003_919,
    -0.000_000_000_000_001_500,
    0.000_000_000_000_000_589,
    -0.000_000_000_000_000_237,
    0.000_000_000_000_000_097,
    -0.000_000_000_000_000_040,
];

/* chebyshev expansions for amplitude and phase
   functions used in bessel evaluations

   These are the same for J0,Y0 and for J1,Y1, so
   they sit outside those functions.
*/

lazy_static! {
    // based on SLATEC besy0, 1980 version, w. fullerton */
    // chebyshev expansions
    // series for by0 on the interval  0. ->  1.60000d+01
    //      with weighted error:         1.20e-17
    //      log weighted error:             16.92
    //      significant figures required:   16.15
    //      decimal places required:        17.48
    //
    pub(super) static ref BY0_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            -0.011_277_839_392_865_573,
            -0.128_345_237_560_420_35,
            -0.104_378_847_997_942_49,
            0.023_662_749_183_969_694,
            -0.002_090_391_647_700_486,
            0.000_103_975_453_939_057,
            -0.000_003_369_747_162_423,
            0.000_000_077_293_842_676,
            -0.000_000_001_324_976_772,
            0.000_000_000_017_648_232,
            -0.000_000_000_000_188_105,
            0.000_000_000_000_001_641,
            -0.000_000_000_000_000_011,
        ],
        a: -1.0,
        b: 1.0,
    };
    pub(super) static ref BY1_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            0.032_080_471_006_119_084,
            1.262_707_897_433_500_4,
            0.006_499_961_899_923_175,
            -0.089_361_645_288_605_05,
            0.013_250_881_221_757_096,
            -0.000_897_905_911_964_835_2,
            0.000_036_473_614_879_583_06,
            -0.000_001_001_374_381_666_00,
            0.000_000_019_945_396_573_90,
            -0.000_000_000_302_306_560_18,
            0.000_000_000_003_609_878_15,
            -0.000_000_000_000_034_874_88,
            0.000_000_000_000_000_278_38,
            -0.000_000_000_000_000_001_86,
        ],
        a: -1.0,
        b: 1.0,
    };
    pub(super) static ref BESSEL_AMP_PHASE_BM0_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            0.092_849_616_373_816_45,
            -0.001_429_877_074_034_84,
            0.000_028_305_792_712_57,
            -0.000_001_433_006_114_24,
            0.000_000_120_286_280_46,
            -0.000_000_013_971_130_13,
            0.000_000_002_040_761_88,
            -0.000_000_000_353_996_69,
            0.000_000_000_070_247_59,
            -0.000_000_000_015_541_07,
            0.000_000_000_003_762_26,
            -0.000_000_000_000_982_82,
            0.000_000_000_000_274_08,
            -0.000_000_000_000_080_91,
            0.000_000_000_000_025_11,
            -0.000_000_000_000_008_14,
            0.000_000_000_000_002_75,
            -0.000_000_000_000_000_96,
            0.000_000_000_000_000_34,
            -0.000_000_000_000_000_12,
            0.000_000_000_000_000_04,
        ],
        a: -1.0,
        b: 1.0,
    };
    pub(super) static ref BESSEL_AMP_PHASE_BTH0_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            -0.246_391_637_743_001_19,
            0.001_737_098_307_508_963,
            -0.000_062_183_633_402_968,
            0.000_004_368_050_165_742,
            -0.000_000_456_093_019_869,
            0.000_000_062_197_400_101,
            -0.000_000_010_300_442_889,
            0.000_000_001_979_526_776,
            -0.000_000_000_428_198_396,
            0.000_000_000_102_035_840,
            -0.000_000_000_026_363_898,
            0.000_000_000_007_297_935,
            -0.000_000_000_002_144_188,
            0.000_000_000_000_663_693,
            -0.000_000_000_000_215_126,
            0.000_000_000_000_072_659,
            -0.000_000_000_000_025_465,
            0.000_000_000_000_009_229,
            -0.000_000_000_000_003_448,
            0.000_000_000_000_001_325,
            -0.000_000_000_000_000_522,
            0.000_000_000_000_000_210,
            -0.000_000_000_000_000_087,
            0.000_000_000_000_000_036,
        ],
        a: -1.0,
        b: 1.0,
    };
    pub(super) static ref BESSEL_AMP_PHASE_BM1_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            0.104_736_251_093_128_5,
            0.004_424_438_937_023_45,
            -0.000_056_616_395_040_35,
            0.000_002_313_494_173_39,
            -0.000_000_173_771_820_07,
            0.000_000_018_932_099_30,
            -0.000_000_002_654_160_23,
            0.000_000_000_447_402_09,
            -0.000_000_000_086_917_95,
            0.000_000_000_018_914_92,
            -0.000_000_000_004_518_84,
            0.000_000_000_001_167_65,
            -0.000_000_000_000_322_65,
            0.000_000_000_000_094_50,
            -0.000_000_000_000_029_13,
            0.000_000_000_000_009_39,
            -0.000_000_000_000_003_15,
            0.000_000_000_000_001_09,
            -0.000_000_000_000_000_39,
            0.000_000_000_000_000_14,
            -0.000_000_000_000_000_05,
        ],
        a: -1.0,
        b: 1.0,
    };
    pub(super) static ref BESSEL_AMP_PHASE_BTH1_CHEB: ChebSeries<f64> = ChebSeries {
        coeffs: vec![
            0.740_601_410_263_138_5,
            -0.004_571_755_659_637_69,
            0.000_119_818_510_964_326,
            -0.000_006_964_561_891_648,
            0.000_000_655_495_621_447,
            -0.000_000_084_066_228_945,
            0.000_000_013_376_886_564,
            -0.000_000_002_499_565_654,
            0.000_000_000_529_495_100,
            -0.000_000_000_124_135_944,
            0.000_000_000_031_656_485,
            -0.000_000_000_008_668_640,
            0.000_000_000_002_523_758,
            -0.000_000_000_000_775_085,
            0.000_000_000_000_249_527,
            -0.000_000_000_000_083_773,
            0.000_000_000_000_029_205,
            -0.000_000_000_000_010_534,
            0.000_000_000_000_003_919,
            -0.000_000_000_000_001_500,
            0.000_000_000_000_000_589,
            -0.000_000_000_000_000_237,
            0.000_000_000_000_000_097,
            -0.000_000_000_000_000_040,
        ],
        a: -1.0,
        b: 1.0,
    };
}