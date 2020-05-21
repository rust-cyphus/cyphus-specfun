mod core;
mod data;
use crate::result::SpecFunResult;

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
        core::airy_ai_e(*self)
    }
    fn airy_ai_scaled_e(&self) -> SpecFunResult<f64> {
        core::airy_ai_scaled_e(*self)
    }
    fn airy_ai(&self) -> f64 {
        core::airy_ai_e(*self).val
    }
    fn airy_ai_scaled(&self) -> f64 {
        core::airy_ai_scaled_e(*self).val
    }
    fn airy_bi_e(&self) -> SpecFunResult<f64> {
        core::airy_bi_e(*self)
    }
    fn airy_bi_scaled_e(&self) -> SpecFunResult<f64> {
        core::airy_bi_scaled_e(*self)
    }
    fn airy_bi(&self) -> f64 {
        core::airy_bi_e(*self).val
    }
    fn airy_bi_scaled(&self) -> f64 {
        core::airy_bi_scaled_e(*self).val
    }
}
