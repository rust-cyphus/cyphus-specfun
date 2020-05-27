pub mod poly;
pub mod airy;
pub mod bessel;
pub mod cheb;
pub(crate) mod consts;
pub mod debye;
pub(crate) mod elementary;
pub mod exp;
pub mod gamma;
pub mod logarithm;
pub mod result;
pub(crate) mod test_utils;
pub mod trig;
pub mod zeta;

#[cfg(test)]
mod test {
    #[test]
    fn mytest() {}
}
