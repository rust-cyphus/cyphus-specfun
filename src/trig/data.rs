// Chebyshev expansion for f(t) = sinc((t+1)/2), -1 < t < 1
pub const SINC_DATA: [f64; 17] = [
    1.133_648_177_811_747_f64,
    -0.532_677_564_732_557_f64,
    -0.068_293_048_346_633_f64,
    0.033_403_684_226_353_f64,
    0.001_485_679_893_925_f64,
    -0.000_734_421_305_768_f64,
    -0.000_016_837_282_388_f64,
    0.000_008_359_950_146_f64,
    0.000_000_117_382_095_f64,
    -0.000_000_058_413_665_f64,
    -0.000_000_000_554_763_f64,
    0.000_000_000_276_434_f64,
    0.000_000_000_001_895_f64,
    -0.000_000_000_000_945_f64,
    -0.000_000_000_000_004_f64,
    0.000_000_000_000_002_f64,
    0.000_000_000_000_000_f64,
];

// Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
// g(x) = (sin(x)/x - 1)/(x*x)
pub const SIN_DATA: [f64; 12] = [
    -0.329_519_016_066_351_f64,
    0.002_537_428_467_166_f64,
    0.000_626_192_878_264_f64,
    -4.649_554_752_185_404e-6_f64,
    -5.691_753_154_937_97e-7_f64,
    3.728_333_514_097_38e-9_f64,
    3.026_737_648_474_747e-10_f64,
    -1.740_087_501_643_662e-12_f64,
    -1.055_467_830_579_084e-13_f64,
    5.370_198_140_913_241e-16_f64,
    2.598_413_798_309_902e-17_f64,
    -1.182_155_525_536_483e-19_f64,
];

// Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
// g(x) = (2(cos(x) - 1)/(x^2) + 1) / x^2
pub const COS_DATA: [f64; 11] = [
    0.165_391_825_637_921_f64,
    -0.000_848_528_838_450_f64,
    -0.000_210_086_507_222_f64,
    1.165_822_696_197_602e-6_f64,
    1.433_193_758_562_598e-7_f64,
    -7.477_088_342_900_714e-10_f64,
    -6.096_999_494_458_425e-11_f64,
    2.907_482_492_019_093e-13_f64,
    1.771_267_398_762_614e-14_f64,
    -7.689_642_150_281_557e-17_f64,
    -3.736_312_113_307_941e-18_f64,
];

/// sinh(x) series
///
/// double-precision for |x| < 1.0
pub fn sinh_series(x: f64) -> f64 {
    let y = x * x;
    let c0: f64 = 1.0 / 6.0;
    let c1: f64 = 1.0 / 120.0;
    let c2: f64 = 1.0 / 5040.0;
    let c3: f64 = 1.0 / 362_880.0;
    let c4: f64 = 1.0 / 39_916_800.0;
    let c5: f64 = 1.0 / 6_227_020_800.0;
    let c6: f64 = 1.0 / 1_307_674_368_000.0;
    let c7: f64 = 1.0 / 355_687_428_096_000.0;

    c7.mul_add(y, c6)
        .mul_add(y, c5)
        .mul_add(y, c4)
        .mul_add(y, c3)
        .mul_add(y, c2)
        .mul_add(y, c1)
        .mul_add(y, c0)
        .mul_add(y, 1.0)
        * x
}

/// cosh(x) - 1 series
///
/// double-precision for |x| < 1.0
pub fn cosh_m1_series(x: f64) -> f64 {
    let y = x * x;
    let c0: f64 = 0.5;
    let c1: f64 = 1.0 / 24.0;
    let c2: f64 = 1.0 / 720.0;
    let c3: f64 = 1.0 / 40320.0;
    let c4: f64 = 1.0 / 3_628_800.0;
    let c5: f64 = 1.0 / 479_001_600.0;
    let c6: f64 = 1.0 / 87_178_291_200.0;
    let c7: f64 = 1.0 / 20_922_789_888_000.0;
    let c8: f64 = 1.0 / 6_402_373_705_728_000.0;

    c8.mul_add(y, c7)
        .mul_add(y, c6)
        .mul_add(y, c5)
        .mul_add(y, c4)
        .mul_add(y, c3)
        .mul_add(y, c2)
        .mul_add(y, c1)
        .mul_add(y, c0)
        * y
}
