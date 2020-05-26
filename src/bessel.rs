pub(crate) mod cyl_bessel_j;
pub(crate) mod cyl_bessel_y;
pub(crate) mod bessel_data;
mod bessel_helpers;
mod olver;

#[cfg(test)]
mod test {
    use crate::bessel::cyl_bessel_j::*;
    use crate::consts::SQRT_DLB_EPS;
    use crate::test_utils::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
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
        test_check_result_and_code!(cyl_bessel_jv_e, (0.0001, 1.0),         0.7652115411876708497,  TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (0.0001, 10.0),       -0.2459270166445205,     TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (0.0009765625, 10.0), -0.2458500798634692,     TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (0.75, 1.0),           0.5586524932048917478,  TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (0.75, 10.0),         -0.04968928974751508135, TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, ( 1.0, 0.001), 0.0004999999375000026,     TOL0, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, ( 1.0,   1.0), 0.4400505857449335160,     TOL0, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, ( 1.75,  1.0), 0.168593922545763170103,     TOL1, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (30.0,   1.0), 3.482869794251482902e-42,  TOL0, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (30.0, 100.0), 0.08146012958117222297,    TOL1, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (10.0,   1.0), 2.6306151236874532070e-10, TOL0, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (10.0, 100.0), -0.05473217693547201474,   TOL2, SpecFunCode::Success);
        test_check_result_and_code!(cyl_bessel_jv_e, (10.2, 100.0), -0.03548919161046526864,   TOL2, SpecFunCode::Success);
    }
}
