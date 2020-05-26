use crate::consts::EUL_GAMMA;
use lazy_static::lazy_static;

/// The maximum x such that gamma(x) is not considered an overflow.
pub(crate) const GAMMA_XMAX: f64 = 171.0;

lazy_static! {

/// Table of all factorials up to 170! (where overflow occurs)
#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref FACT_TABLE: [f64; 171] = [
    1.0,
    1.0,
    2.0,
    6.0,
    24.0,
    120.0,
    720.0,
    5040.0,
    40320.0,
    362_880.0,
    3_628_800.0,
    39_916_800.0,
    479_001_600.0,
    6_227_020_800.0,
    87_178_291_200.0,
    1_307_674_368_000.0,
    20_922_789_888_000.0,
    355_687_428_096_000.0,
    6_402_373_705_728_000.0,
    121_645_100_408_832_000.0,
    2_432_902_008_176_640_000.0,
    51_090_942_171_709_440_000.0,
    1_124_000_727_777_607_680_000.0,
    25_852_016_738_884_976_640_000.0,
    620_448_401_733_239_439_360_000.0,
    15_511_210_043_330_985_984_000_000.0,
    403_291_461_126_605_635_584_000_000.0,
    10_888_869_450_418_352_160_768_000_000.0,
    304_888_344_611_713_860_501_504_000_000.0,
    8_841_761_993_739_701_954_543_616_000_000.0,
    265_252_859_812_191_058_636_308_480_000_000.0,
    8_222_838_654_177_922_817_725_562_880_000_000.0,
    263_130_836_933_693_530_167_218_012_160_000_000.0,
    8_683_317_618_811_886_495_518_194_401_280_000_000.0,
    2.952_327_990_396_041_6e38,
    1.033_314_796_638_614_5e40,
    3.719_933_267_899_012_5e41,
    1.376_375_309_122_634_6e43,
    5.230_226_174_666_011e44,
    2.039_788_208_119_744_4e46,
    8.159_152_832_478_977e47,
    3.345_252_661_316_381e49,
    1.405_006_117_752_88e51,
    6.041_526_306_337_383e52,
    2.658_271_574_788_449e54,
    1.196_222_208_654_801_9e56,
    5.502_622_159_812_089e57,
    2.586_232_415_111_681_8e59,
    1.241_391_559_253_607_3e61,
    6.082_818_640_342_675e62,
    3.041_409_320_171_337_6e64,
    1.551_118_753_287_382_2e66,
    8.065_817_517_094_388e67,
    4.274_883_284_060_025_5e69,
    2.308_436_973_392_414e71,
    1.269_640_335_365_827_6e73,
    7.109_985_878_048_635e74,
    4.052_691_950_487_721_4e76,
    2.350_561_331_282_878_5e78,
    1.386_831_185_456_898_4e80,
    8.320_987_112_741_39e81,
    5.075_802_138_772_248e83,
    3.146_997_326_038_794e85,
    1.982_608_315_404_44e87,
    1.268_869_321_858_841_7e89,
    8.247_650_592_082_472e90,
    5.443_449_390_774_431e92,
    3.647_111_091_818_868e94,
    2.480_035_542_436_830_5e96,
    1.711_224_524_281_413e98,
    1.197_857_166_996_989_2e100,
    8.504_785_885_678_623e101,
    6.123_445_837_688_608_5e103,
    4.470_115_461_512_684_4e105,
    3.307_885_441_519_386e107,
    2.480_914_081_139_54e109,
    1.885_494_701_666_050_4e111,
    1.451_830_920_282_858_7e113,
    1.132_428_117_820_629_7e115,
    8.946_182_130_782_976e116,
    7.156_945_704_626_381e118,
    5.797_126_020_747_368e120,
    4.753_643_337_012_842e122,
    3.945_523_969_720_659e124,
    3.314_240_134_565_353e126,
    2.817_104_114_380_55e128,
    2.422_709_538_367_273_4e130,
    2.107_757_298_379_528e132,
    1.854_826_422_573_984_4e134,
    1.650_795_516_090_846e136,
    1.485_715_964_481_761_5e138,
    1.352_001_527_678_403e140,
    1.243_841_405_464_130_8e142,
    1.156_772_507_081_641_6e144,
    1.087_366_156_656_743e146,
    1.032_997_848_823_906e148,
    9.916_779_348_709_496e149,
    9.619_275_968_248_212e151,
    9.426_890_448_883_248e153,
    9.332_621_544_394_415e155,
    9.332_621_544_394_415e157,
    9.425_947_759_838_36e159,
    9.614_466_715_035_127e161,
    9.902_900_716_486_18e163,
    1.029_901_674_514_562_8e166,
    1.081_396_758_240_291e168,
    1.146_280_563_734_708_4e170,
    1.226_520_203_196_138e172,
    1.324_641_819_451_829e174,
    1.443_859_583_202_493_7e176,
    1.588_245_541_522_743e178,
    1.762_952_551_090_244_6e180,
    1.974_506_857_221_074e182,
    2.231_192_748_659_813_8e184,
    2.543_559_733_472_187_7e186,
    2.925_093_693_493_016e188,
    3.393_108_684_451_898e190,
    3.969_937_160_808_721e192,
    4.684_525_849_754_291e194,
    5.574_585_761_207_606e196,
    6.689_502_913_449_127e198,
    8.094_298_525_273_444e200,
    9.875_044_200_833_601e202,
    1.214_630_436_702_533e205,
    1.506_141_741_511_141e207,
    1.882_677_176_888_926e209,
    2.372_173_242_880_047e211,
    3.012_660_018_457_659_4e213,
    3.856_204_823_625_804e215,
    4.974_504_222_477_287e217,
    6.466_855_489_220_474e219,
    8.471_580_690_878_82e221,
    1.118_248_651_196_004_3e224,
    1.487_270_706_090_685_7e226,
    1.992_942_746_161_518_8e228,
    2.690_472_707_318_050_4e230,
    3.659_042_881_952_549e232,
    5.012_888_748_274_992e234,
    6.917_786_472_619_489e236,
    9.615_723_196_941_089e238,
    1.346_201_247_571_752_6e241,
    1.898_143_759_076_171e243,
    2.695_364_137_888_163e245,
    3.854_370_717_180_073e247,
    5.550_293_832_739_304_4e249,
    8.047_926_057_471_992e251,
    1.174_997_204_390_910_7e254,
    1.727_245_890_454_639e256,
    2.556_323_917_872_865_4e258,
    3.808_922_637_630_57e260,
    5.713_383_956_445_855e262,
    8.627_209_774_233_24e264,
    1.311_335_885_683_452_4e267,
    2.006_343_905_095_682_3e269,
    3.089_769_613_847_350_8e271,
    4.789_142_901_463_394e273,
    7.471_062_926_282_894e275,
    1.172_956_879_426_414_5e278,
    1.853_271_869_493_735e280,
    2.946_702_272_495_038_4e282,
    4.714_723_635_992_061_6e284,
    7.590_705_053_947_219e286,
    1.229_694_218_739_449_4e289,
    2.004_401_576_545_302_6e291,
    3.287_218_585_534_296e293,
    5.423_910_666_131_589e295,
    9.003_691_705_778_438e297,
    1.503_616_514_864_999e300,
    2.526_075_744_973_198_4e302,
    4.269_068_009_004_705e304,
    7.257_415_615_307_999e306,
];

/// Table of all double factorials up to 297!! (where overflow occurs)
#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix,
dead_code
)]
pub(crate) static ref DOUBLE_FACT_TABLE: [f64; 298] = [
    1.000_000_000_000_000_000_000_000_000,
    1.000_000_000_000_000_000_000_000_000,
    2.000_000_000_000_000_000_000_000_000,
    3.000_000_000_000_000_000_000_000_000,
    8.000_000_000_000_000_000_000_000_000,
    15.000_000_000_000_000_000_000_000_00,
    48.000_000_000_000_000_000_000_000_00,
    105.000_000_000_000_000_000_000_000_0,
    384.000_000_000_000_000_000_000_000_0,
    945.000_000_000_000_000_000_000_000_0,
    3_840.000_000_000_000_000_000_000_000,
    10_395.000_000_000_000_000_000_000_00,
    46_080.000_000_000_000_000_000_000_00,
    135_135.000_000_000_000_000_000_000_0,
    645_120.000_000_000_000_000_000_000_00,
    2.027_025e6,
    1.032_192e7,
    3.445_942_5e7,
    1.857_945_6e8,
    6.547_290_75e8,
    3.715_891_2e9,
    1.374_931_057_5e10,
    8.174_960_64e10,
    3.162_341_432_25e11,
    1.961_990_553_6e12,
    7.905_853_580_625e12,
    5.101_175_439_36e13,
    2.134_580_466_768_75e14,
    1.428_329_123_020_8e15,
    6.190_283_353_629_375e15,
    4.284_987_369_062_4e16,
    1.918_987_839_625_106_2e17,
    1.371_195_958_099_968e18,
    6.332_659_870_762_85e18,
    4.662_066_257_539_891e19,
    2.216_430_954_766_997_6e20,
    1.678_343_852_714_360_8e21,
    8.200_794_532_637_892e21,
    6.377_706_640_314_571e22,
    3.198_309_867_728_777_5e23,
    2.551_082_656_125_828_5e24,
    1.311_307_045_768_798_8e25,
    1.071_454_715_572_848e26,
    5.638_620_296_805_835e26,
    4.714_400_748_520_531e27,
    2.537_379_133_562_625_6e28,
    2.168_624_344_319_444_4e29,
    1.192_568_192_774_434_2e30,
    1.040_939_685_273_333_2e31,
    5.843_584_144_594_727e31,
    5.204_698_426_366_666e32,
    2.980_227_913_743_311e33,
    2.706_443_181_710_666_5e34,
    1.579_520_794_283_954_7e35,
    1.461_479_318_123_759_8e36,
    8.687_364_368_561_751e36,
    8.184_284_181_493_056e37,
    4.951_797_690_080_198e38,
    4.746_884_825_265_972e39,
    2.921_560_637_147_317e40,
    2.848_130_895_159_583_4e41,
    1.782_151_988_659_863_4e42,
    1.765_841_154_998_941_5e43,
    1.122_755_752_855_713_8e44,
    1.130_138_339_199_322_6e45,
    7.297_912_393_562_14e45,
    7.458_913_038_715_529e46,
    4.889_601_303_686_634e47,
    5.072_060_866_326_56e48,
    3.373_824_899_543_777_5e49,
    3.550_442_606_428_592e50,
    2.395_415_678_676_082e51,
    2.556_318_676_628_586_5e52,
    1.748_653_445_433_539_8e53,
    1.891_675_820_705_154e54,
    1.311_490_084_075_154_8e55,
    1.437_673_623_735_917e56,
    1.009_847_364_737_869_3e57,
    1.121_385_426_514_015_2e58,
    7.977_794_181_429_167e58,
    8.971_083_412_112_12e59,
    6.462_013_286_957_625e60,
    7.356_288_397_931_94e61,
    5.363_471_028_174_829e62,
    6.179_282_254_262_829_5e63,
    4.558_950_373_948_605e64,
    5.314_182_738_666_033e65,
    3.966_286_825_335_286_5e66,
    4.676_480_810_026_109_3e67,
    3.529_995_274_548_405e68,
    4.208_832_729_023_498e69,
    3.212_295_699_839_048e70,
    3.872_126_110_701_618_5e71,
    2.987_435_000_850_315e72,
    3.639_798_544_059_521e73,
    2.838_063_250_807_799e74,
    3.494_206_602_297_140_4e75,
    2.752_921_353_283_565_2e76,
    3.424_322_470_251_197_4e77,
    2.725_392_139_750_729_5e78,
    3.424_322_470_251_197_3e79,
    2.752_646_061_148_236_6e80,
    3.492_808_919_656_221_4e81,
    2.835_225_442_982_684e82,
    3.632_521_276_442_470_4e83,
    2.976_986_715_131_818e84,
    3.850_472_553_029_018_6e85,
    3.185_375_785_191_045e86,
    4.158_510_357_271_34e87,
    3.472_059_605_858_239_4e88,
    4.574_361_392_998_474_4e89,
    3.853_986_162_502_645_7e90,
    5.123_284_760_158_291e91,
    4.355_004_363_627_989_5e92,
    5.840_544_626_580_451e93,
    5.008_255_018_172_188e94,
    6.775_031_766_833_324e95,
    5.859_658_371_261_46e96,
    7.994_537_484_863_323e97,
    6.972_993_461_801_137e98,
    9.593_444_981_835_987e99,
    8.437_322_088_779_376e100,
    1.170_400_287_783_990_5e102,
    1.037_790_616_919_863_4e103,
    1.451_296_356_852_148_2e104,
    1.297_238_271_149_829e105,
    1.828_633_409_633_706_6e106,
    1.647_492_604_360_283e107,
    2.340_650_764_331_144_5e108,
    2.125_265_459_624_765_3e109,
    3.042_845_993_630_488e110,
    2.784_097_752_108_442e111,
    4.016_556_711_592_244e112,
    3.702_850_010_304_228_4e113,
    5.382_185_993_533_607e114,
    4.998_847_513_910_708e115,
    7.319_772_951_205_705e116,
    6.848_421_094_057_67e117,
    1.010_128_667_266_387_2e119,
    9.519_305_320_740_162e119,
    1.414_180_134_172_942_3e121,
    1.342_222_050_224_362_8e122,
    2.008_135_790_525_578e123,
    1.919_377_531_820_838_8e124,
    2.891_715_538_356_832_2e125,
    2.783_097_421_140_216e126,
    4.221_904_686_000_975_3e127,
    4.091_153_209_076_118e128,
    6.248_418_935_281_443e129,
    6.095_818_281_523_415e130,
    9.372_628_402_922_166e131,
    9.204_685_605_100_357e132,
    1.424_639_517_244_169_2e134,
    1.408_316_897_580_354_7e135,
    2.193_944_856_556_020_4e136,
    2.182_891_191_249_549_7e137,
    3.422_553_976_227_391_5e138,
    3.427_139_170_261_793_3e139,
    5.407_635_282_439_279e140,
    5.449_151_280_716_251_3e141,
    8.652_216_451_902_847e142,
    8.773_133_561_953_163e143,
    1.401_659_065_208_261e145,
    1.430_020_770_598_365_8e146,
    2.298_720_866_941_548_4e147,
    2.359_534_271_487_303_5e148,
    3.815_876_639_122_97e149,
    3.940_422_233_383_797e150,
    6.410_672_753_726_59e151,
    6.659_313_574_418_617e152,
    1.089_814_368_133_520_2e154,
    1.138_742_621_225_583_4e155,
    1.874_480_713_189_655e156,
    1.970_024_734_720_259_3e157,
    3.261_596_440_949_999_6e158,
    3.447_543_285_760_453_7e159,
    5.740_409_736_071_999e160,
    6.102_151_615_796_004e161,
    1.021_792_933_020_815_8e163,
    1.092_285_139_227_484_6e164,
    1.839_227_279_437_468_5e165,
    1.977_036_102_001_747_2e166,
    3.347_393_648_576_192_5e167,
    3.617_976_066_663_197e168,
    6.159_204_313_380_195e169,
    6.693_255_723_326_915e170,
    1.145_612_002_288_716_2e172,
    1.251_638_820_262_133_2e173,
    2.153_750_564_302_786_4e174,
    2.365_597_370_295_431_4e175,
    4.092_126_072_175_294e176,
    4.518_290_977_264_274e177,
    7.856_882_058_576_564e178,
    8.720_301_586_120_05e179,
    1.524_235_119_363_853_6e181,
    1.700_458_809_293_409_6e182,
    2.987_500_833_953_153e183,
    3.349_903_854_308_017e184,
    5.915_251_651_227_243e185,
    6.666_308_670_072_953e186,
    1.183_050_330_245_448_6e188,
    1.339_928_042_684_663_6e189,
    2.389_761_667_095_806_2e190,
    2.720_053_926_649_867_3e191,
    4.875_113_800_875_445e192,
    5.576_110_549_632_228e193,
    1.004_273_442_980_341_6e195,
    1.154_254_883_773_871_3e196,
    2.088_888_761_399_110_6e197,
    2.412_392_707_087_391e198,
    4.386_666_398_938_132e199,
    5.090_148_611_954_395e200,
    9.299_732_765_748_84e201,
    1.084_201_654_346_286_1e203,
    1.990_142_811_870_251_8e204,
    2.331_033_556_844_515e205,
    4.298_708_473_639_744e206,
    5.058_342_818_352_597_5e207,
    9.371_184_472_534_642e208,
    1.107_777_077_219_219e210,
    2.061_660_583_957_621_2e211,
    2.448_187_340_654_474e212,
    4.576_886_496_385_918_6e213,
    5.459_457_769_659_476e214,
    1.025_222_575_190_445_8e216,
    1.228_377_998_173_382e217,
    2.317_003_019_930_407_7e218,
    2.788_418_055_853_577_4e219,
    5.282_766_885_441_329e220,
    6.385_477_347_904_692_4e221,
    1.215_036_383_651_505_8e223,
    1.475_045_267_365_984e224,
    2.818_884_410_071_493e225,
    3.436_855_472_962_743e226,
    6.596_189_519_567_294_5e227,
    8.076_610_361_462_445e228,
    1.556_700_726_617_881_5e230,
    1.914_156_655_666_599_6e231,
    3.704_947_729_350_558e232,
    4.574_834_407_043_173e233,
    8.891_874_550_441_339e234,
    1.102_535_092_097_404_6e236,
    2.151_833_641_206_804e237,
    2.679_160_273_796_693_4e238,
    5.250_474_084_544_601_5e239,
    6.563_942_670_801_899e240,
    1.291_616_624_797_972e242,
    1.621_293_839_688_068_9e243,
    3.203_209_229_498_971e244,
    4.037_021_660_823_292e245,
    8.008_023_073_747_427e246,
    1.013_292_436_866_646_2e248,
    2.018_021_814_584_351_4e249,
    2.563_629_865_272_615e250,
    5.125_775_409_044_253e251,
    6.537_256_156_445_168_4e252,
    1.312_198_504_715_328_7e254,
    1.680_074_832_206_408_2e255,
    3.385_472_142_165_548e256,
    4.351_393_815_414_597_3e257,
    8.802_227_569_630_426e258,
    1.135_713_785_823_209_9e260,
    2.306_183_623_243_171_4e261,
    2.986_927_256_715_042e262,
    6.088_324_765_361_972e263,
    7.915_357_230_294_861_5e264,
    1.619_494_387_586_284_6e266,
    2.113_400_380_488_728e267,
    4.340_244_958_731_243e268,
    5.685_047_023_514_678e269,
    1.171_866_138_857_435_5e271,
    1.540_647_743_372_477_7e272,
    3.187_475_897_692_225e273,
    4.205_968_339_406_864e274,
    8.733_683_959_676_696e275,
    1.156_641_293_336_887_8e277,
    2.410_496_772_870_768e278,
    3.203_896_382_543_179e279,
    6.701_181_028_580_735e280,
    8.938_870_907_295_469e281,
    1.876_330_688_002_606e283,
    2.511_822_724_950_026_8e284,
    5.291_252_540_167_348e285,
    7.108_458_311_608_576e286,
    1.502_715_721_407_527e288,
    2.025_910_618_808_444e289,
    4.297_766_963_225_527_5e290,
    5.814_363_475_980_234e291,
    1.237_756_885_408_951_7e293,
    1.680_351_044_558_287_8e294,
    3.589_494_967_685_96e295,
    4.889_821_539_664_618e296,
    1.048_132_530_564_300_3e298,
    1.432_717_711_121_733e299,
    3.081_509_639_859_043e300,
    4.226_517_247_809_112e301,
    9.121_268_533_982_767e302,
    1.255_275_622_599_306_4e304,
];

#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref LANCZOS_7_C: [f64; 9] = [
    0.999_999_999_999_809_9,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref GSTAR_A_DATA: [f64; 30] = [
    2.167_864_478_664_630_4,
    -0.055_332_490_187_455_84,
    0.018_003_924_314_607_2,
    -0.005_809_192_694_689_377_6,
    0.001_865_236_894_884_003_4,
    -0.000_597_465_241_139_555_3,
    0.000_191_251_699_077_833_55,
    -0.000_061_249_965_469_446_86,
    0.000_019_638_896_331_308_425,
    -6.306_774_125_463_718e-6,
    2.028_869_840_586_139_2e-6,
    -6.538_489_666_083_846e-7,
    2.110_869_805_890_886_5e-7,
    -6.826_071_491_227_495e-8,
    2.210_856_087_588_056_2e-8,
    -7.171_033_193_025_546e-9,
    2.329_089_298_398_540_8e-9,
    -7.574_037_159_850_559e-10,
    2.465_826_722_259_433_3e-10,
    -8.036_224_317_165_988e-11,
    2.621_561_682_634_159_3e-11,
    -8.559_615_502_594_875e-12,
    2.797_083_149_948_796_2e-12,
    -9.147_177_121_188_62e-13,
    2.993_472_019_806_34e-13,
    -9.802_657_590_975_345e-14,
    3.211_677_366_776_715e-14,
    -1.051_803_533_387_814_7e-14,
    3.414_440_572_018_525_3e-15,
    -1.011_515_394_308_118_7e-15,
];

#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref GSTAR_B_DATA: [f64; 30] = [
    0.005_750_227_727_311_434,
    0.000_449_668_953_496_568_5,
    -0.000_167_276_315_318_871_74,
    0.000_061_513_701_491_315_48,
    -0.000_022_372_655_171_152_5,
    8.050_740_535_664_795e-6,
    -2.867_107_710_758_339_6e-6,
    1.010_672_705_374_274_7e-6,
    -3.526_555_847_759_506_4e-7,
    1.217_921_604_641_940_2e-7,
    -4.161_964_018_079_537e-8,
    1.406_628_350_079_520_6e-8,
    -4.698_257_038_053_71e-9,
    1.549_124_866_462_061_4e-9,
    -5.034_093_631_939_488e-10,
    1.608_444_867_373_603_3e-10,
    -5.034_973_319_683_546e-11,
    1.535_715_493_976_213_7e-11,
    -4.523_380_965_577_565e-12,
    1.266_442_917_925_444_8e-12,
    -3.264_828_793_744_932_6e-13,
    7.152_827_272_608_614e-14,
    -9.483_173_525_256_604e-15,
    -2.312_400_199_141_320_8e-15,
    2.840_661_327_717_039e-15,
    -1.724_537_032_161_881_6e-15,
    8.650_792_312_867_111e-16,
    -3.950_656_366_542_755_6e-16,
    1.677_934_213_207_476_2e-16,
    -6.048_315_303_441_477e-17,
];

// Chebyshev expansion for log(gamma(x)/gamma(8)) for 5<x<10 over -1 < t < 1
#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref GAMMA_5_10_DATA: [f64; 24] = [
    -1.528_559_409_666_157_9,
    4.825_915_230_059_59,
    0.227_771_232_097_761_5,
    -0.013_886_766_568_561_788,
    0.001_270_487_649_520_108_3,
    -0.000_139_384_124_025_499_36,
    0.000_016_970_924_299_232_27,
    -2.210_852_882_021_058e-6,
    3.019_660_285_420_231e-7,
    -4.270_567_500_007_911_6e-8,
    6.202_642_381_805_14e-9,
    -9.199_397_320_888_091e-10,
    1.387_555_125_802_814_6e-10,
    -2.121_886_149_190_678_8e-11,
    3.282_173_604_038_144e-12,
    -5.126_000_100_995_379e-13,
    8.071_353_255_487_464e-14,
    -1.279_852_237_656_920_8e-14,
    2.041_771_160_085_25e-15,
    -3.274_523_950_299_235_5e-16,
    5.275_941_842_203_658e-17,
    -8.535_414_715_169_523e-18,
    1.385_863_970_388_807_8e-18,
    -2.257_439_880_773_862_6e-19,
];

/// Chebyshev fit for f(y) = Re(Psi(1 + i*y)) + Euler - y^2/(1+y^2)-y^2/(2(4+y^2))
/// 1 < y < 10 ==>
/// y(x)  = (9x + 11)/2,  -1 < x < 1
/// x(y) = (2y - 11)/9
///
/// g(x) := f(y(x))
#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref R1PY_DATA: [f64; 30] = [
    1.598_883_282_449_769_6,
    0.679_056_253_532_134_6,
    -0.068_485_802_980_122_52,
    -0.005_788_184_183_095_867,
    0.008_511_258_167_108_616,
    -0.004_042_656_134_699_694,
    0.001_352_328_406_159_402_6,
    -0.000_311_646_563_930_660_6,
    0.000_018_507_563_785_249_135,
    0.000_028_348_705_427_529_85,
    -0.000_019_487_536_014_574_535,
    8.070_978_871_083_448e-6,
    -2.298_356_432_134_051_7e-6,
    3.050_662_959_960_475e-7,
    1.304_223_863_241_836_5e-7,
    -1.230_865_718_104_895e-7,
    5.771_085_571_068_243e-8,
    -1.827_555_934_245_096_3e-8,
    3.102_047_130_062_659e-9,
    6.898_932_748_059_381e-10,
    -8.718_229_025_892_306e-10,
    4.406_914_771_024_361_3e-10,
    -1.472_731_109_919_853_6e-10,
    2.758_968_252_326_264_5e-11,
    4.187_182_675_697_586e-12,
    -6.567_346_048_726_008_6e-12,
    3.448_790_088_672_321_3e-12,
    -1.180_725_141_744_869_1e-12,
    2.379_831_434_396_959e-13,
    2.166_363_041_081_883e-15,
];

/// Chebyshev fits from SLATEC code for psi(x)
///
/// Series for PSI        on the interval  0.         to  1.00000D+00
///                                       with weighted error   2.03E-17
///                                        log weighted error  16.69
///                              significant figures required  16.39
///                                   decimal places required  17.37
#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref PSICS_DATA: [f64; 23] = [
    -0.038_057_080_835_217_92,
    0.491_415_393_029_387_14,
    -0.056_815_747_821_244_73,
    0.008_357_821_225_914_313,
    -0.001_333_232_857_994_342,
    0.000_220_313_287_069_308,
    -0.000_037_040_238_178_456,
    0.000_006_283_793_654_854,
    -0.000_001_071_263_908_506,
    0.000_000_183_128_394_654,
    -0.000_000_031_353_509_361,
    0.000_000_005_372_808_776,
    -0.000_000_000_921_168_141,
    0.000_000_000_157_981_265,
    -0.000_000_000_027_098_646,
    0.000_000_000_004_648_722,
    -0.000_000_000_000_797_527,
    0.000_000_000_000_136_827,
    -0.000_000_000_000_023_475,
    0.000_000_000_000_004_027,
    -0.000_000_000_000_000_691,
    0.000_000_000_000_000_118,
    -0.000_000_000_000_000_020,
];

/// Chebyshev fits from SLATEC code for psi(x)
/// Series for APSI       on the interval  0.         to  2.50000D-01
///                                       with weighted error   5.54E-17
///                                        log weighted error  16.26
///                              significant figures required  14.42
///                                   decimal places required  16.86
#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref APSICS_DATA: [f64; 16] = [
    -0.020_474_904_467_818_5,
    -0.010_180_127_153_485_9,
    0.000_055_971_872_538_7,
    -0.000_001_291_717_657_0,
    0.000_000_057_285_860_6,
    -0.000_000_003_821_353_9,
    0.000_000_000_339_743_4,
    -0.000_000_000_037_483_8,
    0.000_000_000_004_899_0,
    -0.000_000_000_000_734_4,
    0.000_000_000_000_123_3,
    -0.000_000_000_000_022_8,
    0.000_000_000_000_004_5,
    -0.000_000_000_000_000_9,
    0.000_000_000_000_000_2,
    -0.000_000_000_000_000_0,
];

#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref PSI_TABLE: [f64; 101] = [
    0.0, /* Infinity */
    /* psi(0) */
    -EUL_GAMMA, /* psi(1) */
    0.422_784_335_098_467_13, /* ...    */
    0.922_784_335_098_467_1,
    1.256_117_668_431_800_5,
    1.506_117_668_431_800_5,
    1.706_117_668_431_800_5,
    1.872_784_335_098_467_2,
    2.015_641_477_955_61,
    2.140_641_477_955_61,
    2.251_752_589_066_721,
    2.351_752_589_066_721,
    2.442_661_679_975_812,
    2.525_995_013_309_145_3,
    2.602_918_090_232_222_4,
    2.674_346_661_660_793_6,
    2.741_013_328_327_460_5,
    2.803_513_328_327_460_5,
    2.862_336_857_739_225_f64,
    2.917_892_413_294_780_8,
    2.970_523_992_242_149,
    3.020_523_992_242_149,
    3.068_143_039_861_196_6,
    3.113_597_585_315_742,
    3.157_075_846_185_307_5,
    3.198_742_512_851_974,
    3.238_742_512_851_974,
    3.277_204_051_313_512_3,
    3.314_241_088_350_549_5,
    3.349_955_374_064_835,
    3.384_438_132_685_525,
    3.417_771_466_018_858_3,
    3.450_029_530_534_987_3,
    3.481_279_530_534_987_3,
    3.511_582_560_838_017_6,
    3.540_994_325_543_9,
    3.569_565_754_115_328_3,
    3.597_343_531_893_106_4,
    3.624_370_558_920_133,
    3.650_686_348_393_817_7,
    3.676_327_374_034_843,
    3.701_327_374_034_843,
    3.725_717_617_937_282,
    3.749_527_141_746_806,
    3.772_782_955_700_294_3,
    3.795_510_228_427_567_2,
    3.817_732_450_649_789_4,
    3.839_471_581_084_572,
    3.860_748_176_829_252_6,
    3.881_581_510_162_586,
    3.901_989_673_427_892,
    3.921_989_673_427_892,
    3.941_597_516_565_147,
    3.960_828_285_795_916_5,
    3.979_696_210_324_218,
    3.998_214_728_842_736_8,
    4.016_396_547_024_555,
    4.034_253_689_881_698,
    4.051_797_549_530_820_5,
    4.069_038_928_841_166,
    4.085_988_081_383_538_5,
    4.102_654_748_050_205,
    4.119_048_190_673_156,
    4.135_177_222_931_221,
    4.151_050_238_804_236_5,
    4.166_675_238_804_236_5,
    4.182_059_854_188_852,
    4.197_211_369_340_366_5,
    4.212_136_742_474_695,
    4.226_842_624_827_636,
    4.241_335_378_450_825,
    4.255_621_092_736_539,
    4.269_705_599_778_792,
    4.283_594_488_667_681,
    4.297_293_118_804_667,
    4.310_806_632_318_181,
    4.324_139_965_651_515,
    4.337_297_860_388_356_5,
    4.350_284_873_375_37,
    4.363_105_386_195_882,
    4.375_763_614_043_984,
    4.388_263_614_043_984,
    4.400_609_293_056_33,
    4.412_804_415_007_549,
    4.424_852_607_778_633,
    4.436_757_369_683_395,
    4.448_522_075_565_748,
    4.460_149_982_542_492,
    4.471_644_235_416_055,
    4.483_007_871_779_692,
    4.494_243_826_835_872,
    4.505_354_937_946_983,
    4.516_343_948_935_994,
    4.527_213_514_153_385,
    4.537_966_202_325_428,
    4.548_604_500_197_769,
    4.559_130_815_987_242,
    4.569_547_482_653_909,
    4.579_856_761_004_424,
    4.590_060_842_637_078,
    4.600_161_852_738_087,
];

#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref PSI_1_TABLE: [f64; 101] = [
    0.0, /* psi(1,0) Infinity */
    std::f64::consts::PI * std::f64::consts::PI / 6.0, /* psi(1,1) */
    0.644_934_066_848_226_4, /* ...      */
    0.394_934_066_848_226_46,
    0.283_822_955_737_115_3,
    0.221_322_955_737_115_33,
    0.181_322_955_737_115_32f64,
    0.153_545_177_959_337_56,
    0.133_137_014_694_031_4,
    0.117_512_014_694_031_43,
    0.105_166_335_681_685_75,
    0.095_166_335_681_685_74,
    0.086_901_872_871_768_38,
    0.079_957_428_427_323_95,
    0.074_040_268_664_010_34,
    0.068_938_227_847_683_8,
    0.064_493_783_403_239_36,
    0.060_587_533_403_239_364,
    0.057_127_325_790_782_61,
    0.054_040_906_037_696_19,
    0.051_270_822_935_203_12,
    0.048_770_822_935_203_12,
    0.046_503_249_239_057_99,
    0.044_437_133_536_578_65,
    0.042_546_774_368_336_69,
    0.040_810_663_257_225_58,
    0.039_210_663_257_225_58,
    0.037_731_373_316_397_18,
    0.036_359_631_203_914_326,
    0.035_084_120_999_832_69,
    0.033_895_060_357_739_946,
    0.032_783_949_246_628_835,
    0.031_743_366_520_302_09,
    0.030_766_804_020_302_09,
    0.029_848_530_374_755_718,
    0.028_983_478_471_641_53,
    0.028_167_151_941_029_287,
    0.027_395_547_002_757_682,
    0.026_665_086_812_838_03,
    0.025_972_566_037_214_76,
    0.025_315_103_841_291_03,
    0.024_690_103_841_291_028,
    0.024_095_219_843_670_565,
    0.023_528_326_419_634_284,
    0.022_987_493_536_995_02,
    0.022_470_964_611_375_183,
    0.021_977_137_450_881_357,
    0.021_504_547_658_820_865,
    0.021_051_854_132_338_295,
    0.020_617_826_354_560_515,
    0.020_201_333_226_697_125,
    0.019_801_333_226_697_127,
    0.019_416_865_714_201_932,
    0.019_047_043_228_994_83,
    0.018_691_044_652_989_135,
    0.018_348_109_124_868_422,
    0.018_017_530_612_471_726,
    0.017_698_653_061_451_318,
    0.017_390_866_050_063_2,
    0.017_093_600_889_540_015,
    0.016_806_327_117_635_387,
    0.016_528_549_339_857_61,
    0.016_259_804_378_825_63,
    0.015_999_658_697_243_943,
    0.015_747_706_064_338_93,
    0.015_503_565_439_338_93,
    0.015_266_879_048_806_387,
    0.015_037_310_637_419_792,
    0.014_814_543_874_220_862,
    0.014_598_280_898_442_315,
    0.014_388_240_990_859_875,
    0.014_184_159_358_206_813,
    0.013_985_786_019_583_524,
    0.013_792_884_785_015_624,
    0.013_605_232_317_385_673,
    0.013_422_617_269_905_76,
    0.013_244_839_492_127_984,
    0.013_071_709_298_222_167,
    0.012_903_046_791_897_322,
    0.012_738_681_242_916_388,
    0.012_578_450_510_661_943,
    0.012_422_200_510_661_943,
    0.012_269_784_720_386_069,
    0.012_121_063_720_980_953,
    0.011_975_904_771_931_745,
    0.011_834_181_415_922_674,
    0.011_695_773_111_424_404,
    0.011_560_564_890_764_59,
    0.011_428_447_041_643_173,
    0.011_299_314_810_238_213,
    0.011_173_068_124_213_722,
    0.011_049_611_334_090_265,
    0.010_928_852_971_573_661,
    0.010_810_705_523_558_537,
    0.010_695_085_220_633_345,
    0.010_581_911_839_012_7,
    0.010_471_108_514_912_978,
    0.010_362_601_570_468_534,
    0.010_256_320_350_360_127, /* ...        */
    0.010_152_197_068_394_28, /* psi(1,99)  */
    0.010_050_166_663_333_571, /* psi(1,100) */
];

#[allow(
clippy::excessive_precision,
clippy::unreadable_literal,
clippy::unseparated_literal_suffix
)]
pub(crate) static ref BERN: [f64; 21] = [
    0.0, /* no element 0 */
    8.333_333_333_333_333e-2,
    -1.388_888_888_888_889e-3,
    3.306_878_306_878_307e-5,
    -8.267_195_767_195_768e-7,
    2.087_675_698_786_81e-8,
    -5.284_190_138_687_493e-10,
    1.338_253_653_068_467_9e-11,
    -3.389_680_296_322_582_7e-13,
    8.586_062_056_277_845e-15,
    -2.174_868_698_558_062e-16,
    5.509_002_828_360_229_5e-18,
    -1.395_446_468_581_252_2e-19,
    3.534_707_039_629_467e-21,
    -8.953_517_427_037_546e-23,
    2.267_952_452_337_683e-24,
    -5.744_724_395_202_645e-25,
    1.455_172_475_614_865e-27,
    -3.685_994_940_665_310_3e-29,
    9.336_734_257_095_045e-31,
    -2.365_022_415_700_63e-32,
];
}