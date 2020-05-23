mod airy_core;
mod airy_data;
mod airy_deriv;
mod airy_zeros;
use crate::result::SpecFunResult;

use airy_core::*;
use airy_deriv::*;

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
}

impl Airy for f64 {
    fn airy_ai_e(&self) -> SpecFunResult<f64> {
        airy_ai_e(*self)
    }
    fn airy_ai_scaled_e(&self) -> SpecFunResult<f64> {
        airy_ai_scaled_e(*self)
    }
    fn airy_ai(&self) -> f64 {
        airy_ai_e(*self).val
    }
    fn airy_ai_scaled(&self) -> f64 {
        airy_ai_scaled_e(*self).val
    }
    fn airy_bi_e(&self) -> SpecFunResult<f64> {
        airy_bi_e(*self)
    }
    fn airy_bi_scaled_e(&self) -> SpecFunResult<f64> {
        airy_bi_scaled_e(*self)
    }
    fn airy_bi(&self) -> f64 {
        airy_bi_e(*self).val
    }
    fn airy_bi_scaled(&self) -> f64 {
        airy_bi_scaled_e(*self).val
    }
}

#[cfg(test)]
mod test {
    use super::airy_deriv::*;
    use super::airy_zeros::*;
    use super::*;
    use crate::result::SpecFunCode;
    use crate::test_check_result_and_code;
    use crate::test_utils::*;
    const TOL0: f64 = 2.0 * f64::EPSILON;
    const SQRT_TOL0: f64 = 2.0 * crate::consts::SQRT_DLB_EPS;
    const TOL1: f64 = 16.0 * f64::EPSILON;
    const TOL2: f64 = 256.0 * f64::EPSILON;
    const TOL3: f64 = 2048.0 * f64::EPSILON;
    const TOL4: f64 = 16384.0 * f64::EPSILON;
    const TOL5: f64 = 131072.0 * f64::EPSILON;

    #[test]
    fn test_airy_ai_e() {
        test_check_result_and_code!(
            airy_ai_e,
            (-500.0),
            0.0725901201040411396,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (-5.0),
            0.3507610090241142,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (-0.3000000000000094),
            0.4309030952855831,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (0.6999999999999907),
            0.1891624003981519,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            airy_ai_e,
            (1.649999999999991),
            0.0583105861872088521,
            TOL0,
            SpecFunCode::Success
        );

        test_check_result_and_code!(
            airy_ai_e,
            (2.54999999999999),
            0.01446149513295428,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (3.499999999999987),
            0.002584098786989702,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_e,
            (5.39999999999998),
            4.272986169411866e-05,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_ai_scaled_e() {
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (-5.0),
            0.3507610090241142,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (0.6999999999999907),
            0.2795125667681217,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (1.649999999999991),
            0.2395493001442741,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (2.54999999999999),
            0.2183658595899388,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (3.499999999999987),
            0.2032920808163519,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_scaled_e,
            (5.39999999999998),
            0.1836050093282229,
            TOL0,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_airy_bi_e() {
        test_check_result_and_code!(
            airy_bi_e,
            (-500.0),
            -0.094688570132991028,
            TOL4,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (-5.0),
            -0.1383691349016005,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (0.6999999999999907),
            0.9733286558781599,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (1.649999999999991),
            2.196407956850028,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (2.54999999999999),
            6.973628612493443,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (3.499999999999987),
            33.05550675461069,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_e,
            (5.39999999999998),
            1604.476078241272,
            TOL1,
            SpecFunCode::Success
        );
    }
    #[test]
    fn test_airy_bi_scaled_e() {
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (-5.0),
            -0.1383691349016005,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (0.6999999999999907),
            0.6587080754582302,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (1.649999999999991),
            0.5346449995597539,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (2.54999999999999),
            0.461835455542297,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (3.499999999999987),
            0.4201771882353061,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_scaled_e,
            (5.39999999999998),
            0.3734050675720473,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_ai_deriv_e() {
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (-5.0),
            0.3271928185544435,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (-0.5500000000000094),
            -0.1914604987143629,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (0.4999999999999906),
            -0.2249105326646850,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (1.899999999999992),
            -0.06043678178575718,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (3.249999999999988),
            -0.007792687926790889,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_e,
            (5.199999999999981),
            -0.0001589434526459543,
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_ai_deriv_scaled_e() {
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (-5.0),
            0.3271928185544435,
            TOL1,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (0.5499999999999906),
            -0.2874057279170166,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (1.499999999999991),
            -0.3314199796863637,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (2.49999999999999),
            -0.3661089384751620,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (3.649999999999986),
            -0.3974033831453963,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_ai_deriv_scaled_e,
            (6.299999999999977),
            -0.4508799189585947,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_bi_deriv_e() {
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (-5.0),
            0.778411773001899,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (-0.5500000000000094),
            0.5155785358765014,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (0.4999999999999906),
            0.5445725641405883,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (1.899999999999992),
            3.495165862891568,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (3.249999999999988),
            36.55485149250338,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_e,
            (5.199999999999981),
            2279.748293583233,
            TOL1,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_bi_deriv_scaled_e() {
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (-5.0),
            0.778411773001899,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (0.5499999999999906),
            0.4322811281817566,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (1.499999999999991),
            0.5542307563918037,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (2.49999999999999),
            0.6755384441644985,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (3.649999999999986),
            0.7613959373000228,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_bi_deriv_scaled_e,
            (6.299999999999977),
            0.8852064139737571,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_zero_ai_e() {
        test_check_result_and_code!(
            airy_zero_ai_e,
            (2),
            -4.087949444130970617,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_e,
            (50),
            -38.02100867725525443,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_e,
            (100),
            -60.45555727411669871,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_e,
            (110),
            -64.43135670991324811,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_zero_bi_e() {
        test_check_result_and_code!(
            airy_zero_bi_e,
            (2),
            -3.271093302836352716,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (50),
            -37.76583438165180116,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (100),
            -60.25336482580837088,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (110),
            -64.2355167606561537,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (111),
            -64.6268994819519378,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_e,
            (200),
            -95.88699147356682665,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_zero_ai_deriv_e() {
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (2),
            -3.248197582179836561,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (50),
            -37.76565910053887108,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (100),
            -60.25329596442479317,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (110),
            -64.23545617243546956,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_ai_deriv_e,
            (1000),
            -280.9378080358935071,
            TOL0,
            SpecFunCode::Success
        );
    }

    #[test]
    fn test_airy_zero_bi_deriv_e() {
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (2),
            -4.073155089071828216,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (50),
            -38.02083574095788210,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (100),
            -60.45548887257140819,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (110),
            -64.43129648944845060,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (111),
            -64.82208737584206093,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (200),
            -96.04731050310324450,
            TOL0,
            SpecFunCode::Success
        );
        test_check_result_and_code!(
            airy_zero_bi_deriv_e,
            (1000),
            -281.0315164471118527,
            TOL0,
            SpecFunCode::Success
        );
    }
}
