/// Euler–Mascheroni constant
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const EUL_GAMMA: f64 = 0.57721566490153286061;

// 3.0_f64.sqrt()
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const SQRT_3: f64 = 1.7320508075688772935;

// Log constants

/// Natural log of pi
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const LNPI: f64 = 1.1447298858494001741;
/// (2.0*f64::const::PI).sqrt().ln()
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const LN_SQRT_PI: f64 = 0.91893853320467274178;
/// Natural log of f64::MAX
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const LN_DBL_MAX: f64 = 7.0978271289338397e+02;
/// Natural log of f64::MIN
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const LN_DBL_MIN: f64 = -7.0839641853226408e+02;
/// Natural log of f64::EPSILON
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const LN_DBL_EPS: f64 = -3.6043653389117154e+01;

// Roots of f64::EPSILON

/// f64::EPSILON^(1/2)
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const SQRT_DLB_EPS: f64 = 1.4901161193847656e-08;
/// f64::EPSILON^(1/3)
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT3_DBL_EPS: f64 = 6.0554544523933429e-06;
/// f64::EPSILON^(1/4)
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT4_DBL_EPS: f64 = 1.2207031250000000e-04;
/// f64::EPSILON^(1/5)
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT5_DBL_EPS: f64 = 7.4009597974140505e-04;
/// f64::EPSILON^(1/6)
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT6_DBL_EPS: f64 = 2.4607833005759251e-03;

// Roots of f64::MAX

/// f64::MAX.sqrt()
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const SQRT_DBL_MAX: f64 = 1.3407807929942596e+154;
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT3_DBL_MAX: f64 = 5.6438030941222897e+102;
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT4_DBL_MAX: f64 = 1.1579208923731620e+77;
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT5_DBL_MAX: f64 = 4.4765466227572707e+61;
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT6_DBL_MAX: f64 = 2.3756689782295612e+51;

// Roots of f64::MIN_POSITIVE

/// f64::MIN.sqrt()
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const SQRT_DBL_MIN: f64 = 1.4916681462400413e-154;
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT3_DBL_MIN: f64 = 2.8126442852362996e-103;
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT4_DBL_MIN: f64 = 1.2213386697554620e-77;
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT5_DBL_MIN: f64 = 2.9476022969691763e-62;
#[allow(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    clippy::unseparated_literal_suffix
)]
pub(crate) const ROOT6_DBL_MIN: f64 = 5.3034368905798218e-52;
