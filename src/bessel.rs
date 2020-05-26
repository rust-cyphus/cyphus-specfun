pub(crate) mod besselj;
pub(crate) mod bessely;
pub(crate) mod data;
mod helpers;
mod olver;

#[cfg(test)]
mod test {
    use crate::bessel::besselj::*;
    use crate::consts::SQRT_DLB_EPS;
    use crate::test_utils::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;

    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * SQRT_DLB_EPS;
    const TOL4: f64 = 16384.0 * f64::EPSILON;

    #[test]
    fn test_besselj0_e() {
        test_check_result_and_code!(
            besselj0_e,
            (0.1),
            0.997_501_562_066_04,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj0_e,
            (2.0),
            0.223_890_779_141_235_67,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj0_e,
            (100.0),
            0.019_985_850_304_223_122,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj0_e,
            (1.0e+10),
            2.175_591_750_246_892e-6,
            SQRT_TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_besselj1_e() {
        test_check_result_and_code!(
            besselj1_e,
            (0.1),
            0.049_937_526_036_242,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj1_e,
            (2.0),
            0.576_724_807_756_873_4,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj1_e,
            (100.0),
            -0.077_145_352_014_112_16,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besselj1_e,
            (1.0e+10),
            -7.676_508_175_684_158e-6,
            TOL4,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_besseljn_e() {
        test_check_result_and_code!(
            besseljn_e,
            (4, 0.1),
            2.602_864_854_568_403e-7,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besseljn_e,
            (5, 2.0),
            0.007_039_629_755_871_685,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besseljn_e,
            (10, 20.0),
            0.186_482_558_023_945_1,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            besseljn_e,
            (100, 100.0),
            0.096_366_673_295_861_56,
            TOL0,
            SpecFunCode::Success
        );
    }
}
