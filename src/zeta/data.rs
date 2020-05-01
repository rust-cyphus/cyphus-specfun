// Chebyshev fit for (s[t]-1) * Zeta[s(t)], where
// s[t] = (t+1)/2, -1 < t < 1
pub(crate) const ZETA_XLT1_DATA: [f64; 14] = [
    1.480_186_771_569_315_61,
    0.250_120_625_398_894_26,
    0.009_911_375_021_353_60,
    -0.000_120_847_596_566_76,
    -4.758_586_636_766_255_65e-06,
    2.222_994_669_446_639_18e-07,
    -2.223_749_649_803_025_71e-09,
    -1.017_322_651_322_902_83e-10,
    4.375_664_345_042_455_82e-12,
    -6.222_963_259_310_055_14e-14,
    -6.611_620_100_327_220_71e-16,
    4.947_727_953_337_391_23e-17,
    -1.042_981_909_345_618_97e-18,
    6.992_521_616_658_002_10e-21,
];

// Chebyshev fit for (s[t] - 1) * Zeta[s(t)], where
// s[t] = (19t+21)/2, -1 < t < 1
pub(crate) const ZETA_XGT1_DATA: [f64; 30] = [
    19.391_851_572_672_411_94,
    9.152_532_969_251_075_61,
    0.242_789_765_886_737_99,
    -0.133_900_068_826_202_73,
    0.057_782_706_406_502_85,
    -0.018_762_598_375_400_22,
    0.003_940_301_425_832_03,
    -0.000_058_150_827_315_81,
    -0.000_375_614_890_721_48,
    0.000_189_253_054_810_92,
    -0.000_054_903_219_969_55,
    8.708_648_400_893_903_86e-6,
    6.460_947_792_481_188_90e-7,
    -9.674_977_391_505_908_92e-7,
    3.658_540_076_676_725_77e-7,
    -8.459_251_642_727_516_43e-8,
    9.995_678_614_449_793_65e-9,
    1.426_003_642_095_111_81e-9,
    -1.176_196_882_338_287_91e-9,
    3.711_457_589_978_520_46e-10,
    -7.475_685_519_421_096_16e-11,
    7.853_693_420_918_370_04e-12,
    9.982_718_225_968_553_96e-13,
    -7.527_668_703_019_222_15e-13,
    2.195_502_639_396_427_99e-13,
    -4.193_485_985_283_464_74e-14,
    4.634_114_963_593_355_07e-15,
    2.374_248_850_904_834_01e-16,
    -2.727_651_638_812_478_61e-16,
    7.847_357_013_463_604_47e-17,
];

// Chebyshev fit for Ln[Zeta[s[t]]] - 1 - 2^(-s[t]), where
// s[t] = 10+5t -1 <= t <= 1; 5 <= s <= 15
pub(crate) const ZETAM1_INTER_DATA: [f64; 24] = [
    -21.750_943_565_308_848_34,
    -5.630_368_776_981_217_82,
    0.052_804_135_868_422_943_5,
    -0.015_638_180_917_967_076_2,
    0.004_082_184_743_723_559_27,
    -0.001_026_486_734_947_488_2,
    0.000_260_469_880_409_882_387,
    -0.000_067_617_584_720_991_944_7,
    0.000_017_928_447_258_787_112_4,
    -4.832_386_513_185_561_88e-6,
    1.319_137_889_649_992_88e-6,
    -3.637_605_006_563_299_72e-7,
    1.011_468_475_131_947_44e-7,
    -2.832_152_251_418_065_01e-8,
    7.977_337_102_520_214_23e-9,
    -2.258_501_685_539_568_86e-9,
    6.422_693_929_501_643_06e-10,
    -1.833_638_618_461_272_84e-10,
    5.253_097_638_952_831_79e-11,
    -1.509_586_870_425_898_21e-11,
    4.349_975_455_160_492_44e-12,
    -1.255_977_827_481_904_16e-12,
    3.612_807_400_722_226_50e-13,
    -9.664_372_392_057_452_07e-14,
];

// Zeta[n] - 1
pub(crate) const ZETA_POS_TABLE_NMAX: usize = 100;

pub(crate) const ZETAM1_POS_INT_TABLE: [f64; ZETA_POS_TABLE_NMAX + 1] = [
    -1.5, // zeta(0) */
    0.0,  // FIXME: Infinity */
    // zeta(1) - 1 */
    0.644_934_066_848_226_436, // zeta(2) - 1 */
    0.202_056_903_159_594_281,
    0.082_323_233_711_138_191,
    0.036_927_755_143_369_927,
    0.017_343_061_984_449_130,
    0.008_349_277_381_922_829,
    0.004_077_356_197_944_338,
    0.002_008_392_826_082_212,
    0.000_994_575_127_818_080,
    0.000_494_188_604_119_466,
    0.000_246_086_553_308_047,
    0.000_122_713_347_578_486,
    0.000_061_248_135_058_705,
    0.000_030_588_236_307_020,
    0.000_015_282_259_408_657,
    7.637_197_637_899_762_27e-6,
    3.817_293_264_999_839_85e-6,
    1.908_212_716_553_938_92e-6,
    9.539_620_338_727_961_13e-7,
    4.769_329_867_878_064_63e-7,
    2.384_505_027_277_329_90e-7,
    1.192_199_259_653_110_73e-7,
    5.960_818_905_125_947_96e-8,
    2.980_350_351_465_228_01e-8,
    1.490_155_482_836_504_12e-8,
    7.450_711_789_835_429_49e-9,
    3.725_334_024_788_457_05e-9,
    1.862_659_723_513_049_00e-9,
    9.313_274_324_196_681_82e-10,
    4.656_629_065_033_784_07e-10,
    2.328_311_833_676_505_49e-10,
    1.164_155_017_270_051_97e-10,
    5.820_772_087_902_700_88e-11,
    2.910_385_044_497_099_68e-11,
    1.455_192_189_104_198_42e-11,
    7.275_959_835_057_481_01e-12,
    3.637_979_547_378_651_19e-12,
    1.818_989_650_307_065_94e-12,
    9.094_947_840_263_889_28e-13,
    4.547_473_783_042_154_02e-13,
    2.273_736_845_824_652_51e-13,
    1.136_868_407_680_227_84e-13,
    5.684_341_987_627_585_60e-14,
    2.842_170_976_889_301_85e-14,
    1.421_085_482_803_160_67e-14,
    7.105_427_395_210_852_71e-15,
    3.552_713_691_337_113_67e-15,
    1.776_356_843_579_120_32e-15,
    8.881_784_210_930_815_90e-16,
    4.440_892_103_143_813_36e-16,
    2.220_446_050_798_041_98e-16,
    1.110_223_025_141_066_13e-16,
    5.551_115_124_845_481_24e-17,
    2.775_557_562_136_124_17e-17,
    1.387_778_780_972_523_27e-17,
    6.938_893_904_544_153_69e-18,
    3.469_446_952_165_922_62e-18,
    1.734_723_476_047_576_57e-18,
    8.673_617_380_119_933_72e-19,
    4.336_808_690_020_650_48e-19,
    2.168_404_344_997_219_78e-19,
    1.084_202_172_494_241_40e-19,
    5.421_010_862_456_645_41e-20,
    2.710_505_431_223_468_83e-20,
    1.355_252_715_610_116_45e-20,
    6.776_263_578_045_189_09e-21,
    3.388_131_789_020_796_81e-21,
    1.694_065_894_509_799_16e-21,
    8.470_329_472_546_998_34e-22,
    4.235_164_736_272_833_34e-22,
    2.117_582_368_136_194_73e-22,
    1.058_791_184_068_023_38e-22,
    5.293_955_920_339_870_32e-23,
    2.646_977_960_169_852_96e-23,
    1.323_488_980_084_899_08e-23,
    6.617_444_900_424_404_06e-24,
    3.308_722_450_212_171_58e-24,
    1.654_361_225_106_075_64e-24,
    8.271_806_125_530_344_40e-25,
    4.135_903_062_765_160_92e-25,
    2.067_951_531_382_576_70e-25,
    1.033_975_765_691_287_09e-25,
    5.169_878_828_456_431_32e-26,
    2.584_939_414_228_214_26e-26,
    1.292_469_707_114_106_67e-26,
    6.462_348_535_570_531_80e-27,
    3.231_174_267_785_265_38e-27,
    1.615_587_133_892_632_52e-27,
    8.077_935_669_463_162_03e-28,
    4.038_967_834_731_580_82e-28,
    2.019_483_917_365_790_34e-28,
    1.009_741_958_682_895_15e-28,
    5.048_709_793_414_475_69e-29,
    2.524_354_896_707_237_82e-29,
    1.262_177_448_353_618_90e-29,
    6.310_887_241_768_094_49e-30,
    3.155_443_620_884_047_23e-30,
    1.577_721_810_442_023_61e-30,
    7.888_609_052_210_118_07e-31,
];

pub(crate) const ZETA_NEG_TABLE_NMAX: usize = 99;
pub(crate) const ZETA_NEG_TABLE_SIZE: usize = 50;

pub(crate) const ZETA_NEG_INT_TABLE: [f64; ZETA_NEG_TABLE_SIZE] = [
    -0.083_333_333_333_333_333, /* zeta(-1) */
    0.008_333_333_333_333_333,  /* zeta(-3) */
    -0.003_968_253_968_253_968, /* ...      */
    0.004_166_666_666_666_667,
    -0.007_575_757_575_757_576,
    0.021_092_796_092_796_093,
    -0.083_333_333_333_333_333,
    0.443_259_803_921_568_62,
    -3.053_954_330_270_119_74,
    26.456_212_121_212_121_21,
    -281.460_144_927_536_231_88,
    3_607.510_546_398_046_398_04,
    -54_827.583_333_333_333_333_33,
    974_936.823_850_574_712_643_678_160_92,
    -2.005_269_579_668_807_89e+07,
    4.723_848_677_216_299_01e+08,
    -1.263_572_479_591_666_66e+10,
    3.808_793_112_524_536_88e+11,
    -1.285_085_049_930_508_33e+13,
    4.824_144_835_485_017_03e+14,
    -2.004_031_065_651_625_27e+16,
    9.167_743_603_195_330_77e+17,
    -4.597_988_834_365_650_34e+19,
    2.518_047_192_145_109_56e+21,
    -1.500_173_349_215_392_87e+23,
    9.689_957_887_463_594_06e+24,
    -6.764_588_237_929_282_09e+26,
    5.089_065_946_866_228_96e+28,
    -4.114_728_879_255_797_86e+30,
    3.566_658_209_537_555_61e+32,
    -3.306_608_987_657_757_67e+34,
    3.271_563_423_647_871_62e+36,
    -3.447_378_255_827_805_38e+38,
    3.861_427_983_270_525_88e+40,
    -4.589_297_443_245_433_21e+42,
    5.777_538_634_277_043_18e+44,
    -7.691_985_875_950_713_51e+46,
    1.081_363_544_997_165_46e+49,
    -1.602_936_452_200_896_54e+51,
    2.501_947_904_156_046_28e+53,
    -4.106_705_233_581_021_24e+55,
    7.079_877_440_849_458_06e+57,
    -1.280_454_688_793_950_87e+60,
    2.426_734_039_233_352_40e+62,
    -4.814_321_887_404_576_93e+64,
    9.987_557_417_572_753_06e+66,
    -2.164_563_486_843_518_56e+69,
    4.896_232_703_962_055_32e+71,  /* ...        */
    -1.154_902_392_396_351_96e+74, /* zeta(-97)  */
    2.838_224_957_069_370_69e+76,  /* zeta(-99)  */
];

// coefficients for Maclaurin summation in hzeta()
pub(crate) const HZETA_C: [f64; 15] = [
    1.000_000_000_000_000_00,
    0.083_333_333_333_333_333,
    -0.001_388_888_888_888_888_89,
    0.000_033_068_783_068_783_069,
    -8.267_195_767_195_767_19e-07,
    2.087_675_698_786_809_89e-08,
    -5.284_190_138_687_493_18e-10,
    1.338_253_653_068_467_88e-11,
    -3.389_680_296_322_582_86e-13,
    8.586_062_056_277_844_56e-15,
    -2.174_868_698_558_061_87e-16,
    5.509_002_828_360_229_51e-18,
    -1.395_446_468_581_252_33e-19,
    3.534_707_039_629_467_47e-21,
    -8.953_517_427_037_546_85e-23,
];

pub(crate) const ETA_POS_TABLE_NMAX: usize = 100;

pub(crate) const ETA_POS_INT_TABLE: [f64; ETA_POS_TABLE_NMAX + 1] = [
    0.500_000_000_000_000_00, /* eta(0) */
    std::f64::consts::LN_2,   /* eta(1) */
    0.822_467_033_424_113_21, /* ...    */
    0.901_542_677_369_695_71,
    0.947_032_829_497_245_91,
    0.972_119_770_446_909_30,
    0.985_551_091_297_435_10,
    0.992_593_819_922_830_28,
    0.996_233_001_852_647_89,
    0.998_094_297_541_605_33,
    0.999_039_507_598_271_56,
    0.999_517_143_498_060_75,
    0.999_757_685_143_858_19,
    0.999_878_542_763_265_11,
    0.999_939_170_345_979_71,
    0.999_969_551_213_099_23,
    0.999_984_764_214_906_10,
    0.999_992_378_292_041_01,
    0.999_996_187_869_610_11,
    0.999_998_093_508_171_67,
    0.999_999_046_611_581_52,
    0.999_999_523_258_215_54,
    0.999_999_761_613_230_82,
    0.999_999_880_801_318_43,
    0.999_999_940_398_892_39,
    0.999_999_970_198_856_96,
    0.999_999_985_099_231_99,
    0.999_999_992_549_550_48,
    0.999_999_996_274_753_40,
    0.999_999_998_137_369_41,
    0.999_999_999_068_682_28,
    0.999_999_999_534_340_33,
    0.999_999_999_767_169_89,
    0.999_999_999_883_584_85,
    0.999_999_999_941_792_39,
    0.999_999_999_970_896_18,
    0.999_999_999_985_448_09,
    0.999_999_999_992_724_04,
    0.999_999_999_996_362_02,
    0.999_999_999_998_181_01,
    0.999_999_999_999_090_50,
    0.999_999_999_999_545_25,
    0.999_999_999_999_772_62,
    0.999_999_999_999_886_31,
    0.999_999_999_999_943_15,
    0.999_999_999_999_971_57,
    0.999_999_999_999_985_78,
    0.999_999_999_999_992_89,
    0.999_999_999_999_996_44,
    0.999_999_999_999_998_22,
    0.999_999_999_999_999_11,
    0.999_999_999_999_999_55,
    0.999_999_999_999_999_77,
    0.999_999_999_999_999_88,
    0.999_999_999_999_999_94,
    0.999_999_999_999_999_97,
    0.999_999_999_999_999_98,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    0.999_999_999_999_999_99,
    1.000_000_000_000_000_00,
    1.000_000_000_000_000_00,
    1.000_000_000_000_000_00,
];

pub(crate) const ETA_NEG_TABLE_NMAX: usize = 99;
pub(crate) const ETA_NEG_TABLE_SIZE: usize = 50;

// 0.250_000_000_000_000_00
// (\d.\d\d\d)(\d\d\d)(\d\d\d)(\d\d\d)(\d\d\d)(\d\d)(\d\d\d\d\d\d\d\d\d\d\d\d)

pub(crate) const ETA_NEG_INT_TABLE: [f64; ETA_NEG_TABLE_SIZE] = [
    0.250_000_000_000_000_00,  /* eta(-1) */
    -0.125_000_000_000_000_00, /* eta(-3) */
    0.250_000_000_000_000_00,  /* ...      */
    -1.062_500_000_000_000_00,
    7.750_000_000_000_000_00,
    -86.375_000_000_000_000_00,
    1_365.250_000_000_000_000_00,
    -29_049.031_250_000_000_000_00,
    800_572.750_000_000_000_000_00,
    -2.774_132_262_500_000_00e+7,
    1.180_529_130_250_000_00e+9,
    -6.052_398_005_168_750_00e+10,
    3.679_416_778_537_750_00e+12,
    -2.617_076_099_065_838_75e+14,
    2.153_141_814_080_029_52e+16,
    -2.028_877_557_517_301_59e+18,
    2.170_800_990_262_377_05e+20,
    -2.617_382_696_845_581_49e+22,
    3.532_414_887_686_387_78e+24,
    -5.304_203_340_686_490_66e+26,
    8.813_821_836_431_157_67e+28,
    -1.612_806_510_749_077_85e+31,
    3.235_547_000_172_273_42e+33,
    -7.087_672_747_653_749_31e+35,
    1.689_045_034_129_396_57e+38,
    -4.363_969_073_121_683_11e+40,
    1.218_599_882_706_126_13e+43,
    -3.667_058_480_315_300_61e+45,
    1.185_989_852_630_209_91e+48,
    -4.112_076_949_358_401_50e+50,
    1.524_904_243_678_762_03e+53,
    -6.034_969_319_694_130_70e+55,
    2.543_716_176_421_069_58e+58,
    -1.139_692_380_263_228_78e+61,
    5.418_086_106_475_397_91e+63,
    -2.728_365_479_999_437_38e+66,
    1.452_975_051_491_854_32e+69,
    -8.170_551_937_106_745_00e+71,
    4.844_578_160_667_836_77e+74,
    -3.024_669_420_664_951_93e+77,
    1.985_880_796_169_049_30e+80,
    -1.369_447_462_072_008_69e+83,
    9.907_038_298_429_580_78e+85,
    -7.510_378_079_659_264_59e+88,
    5.959_841_826_426_088_08e+91,
    -4.945_598_888_750_002_03e+94,
    4.287_359_692_702_024_12e+97,
    -3.879_195_203_771_616_29e+100,
    3.660_031_777_315_634_22e+103,
    -3.597_877_570_411_728_38e+106, /* eta(-99)  */
];
