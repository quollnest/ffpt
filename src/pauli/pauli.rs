
///
/// Original code by Jannis Ruh
///
///

macro_rules! const_pauli {
    ($($name:ident,)*) => {$(
        /// Pauli
        #[doc = stringify!($name)]
        /// .
        const $name: Self;
    )*};
}

macro_rules! new_pauli {
    ($(($name:ident, $gate:ident),)*) => {$(
        /// Create a new
        #[doc = stringify!($gate)]
        /// Pauli.
        fn $name() -> Self where Self: Sized {
            Self::$gate
        }
    )*};
}

macro_rules! plus {
    ($(($name:ident, $left:ident, $right:ident),)*) => {$(
        /// Add `other`'s
        #[doc = stringify!($right)]
        /// component onto `self`'s
        #[doc = stringify!($left)]
        /// component in place.
        fn $name(&mut self, other: &Self);
    )*};
}

/// The interface we need for the Pauli tracking.
///
/// Note that we only implement some of the gate conjugations, since many are redundant;
/// also you may want to implement some of the default gate conjugations directly for
/// performance reasons; compare the documentation of [Tracker].
///
/// [Tracker]: crate::tracker::Tracker
pub trait Pauli {
    
    const_pauli!(I, X, Y, Z,);
    new_pauli!((new_i, I), (new_x, X), (new_y, Y), (new_z, Z),);

    /// Create a the new Pauli (X if x) * (Z if z), neglecting phases.
    ///
    /// # Examples
    /// ```
    /// # fn main() { #![cfg_attr(coverage_nightly, coverage(off))]
    /// # use pauli_tracker::pauli::{Pauli, PauliDense};
    /// assert_eq!(PauliDense::new_product(false, false), PauliDense::new_i());
    /// assert_eq!(PauliDense::new_product(false, true), PauliDense::new_x());
    /// assert_eq!(PauliDense::new_product(true, false), PauliDense::new_z());
    /// assert_eq!(PauliDense::new_product(true, true), PauliDense::new_y());
    /// # }
    /// ```
    fn new_product(z: bool, x: bool) -> Self;

    /// Multiply `self` with `other` in place (i.e., adding on the tableau
    /// representation).
    fn multiply(&mut self, other: Self)
    where
        // this is a potential breaking change, but much better than before where we had
        // no default implementation
        Self: Sized,
    {
        #[allow(deprecated)]
        self.add(other);
    }

    /// Add the `other` Pauli to `self` in place.
    #[deprecated(since = "0.4.2", note = "use `multiply` instead")]
    fn add(&mut self, other: Self);

    /// Conjugate the Pauli with the I (identity gate). This does nothing!
    #[inline(always)]
    fn id(&mut self) {}

    /// Conjugate the Pauli with the S gate ignoring phases.
    fn s(&mut self);
    /// Conjugate the Pauli with the H gate ignoring phases.
    fn h(&mut self);
    /// Conjugate the Pauli with the SH gate ignoring phases.
    fn sh(&mut self) {
        self.h();
        self.s();
    }
    /// Conjugate the Pauli with the HS gate ignoring phases.
    fn hs(&mut self) {
        self.s();
        self.h();
    }
    /// Conjugate the Pauli with the SHS gate ignoring phases.
    fn shs(&mut self) {
        self.s();
        self.h();
        self.s();
    }

    plus!((xpx, X, X), (xpz, X, Z), (zpx, Z, X), (zpz, Z, Z),);

    /// Get the Pauli's X component.
    ///
    /// # Examples
    /// ```
    /// # fn main() { #![cfg_attr(coverage_nightly, coverage(off))]
    /// # use pauli_tracker::pauli::{Pauli, PauliDense};
    /// let pauli = PauliDense::new_y();
    /// assert_eq!(pauli.get_x(), true);
    /// # }
    /// ```
    fn get_x(&self) -> bool;

    /// Get the Pauli's Z component.
    ///
    /// # Examples
    /// ```
    /// # fn main() { #![cfg_attr(coverage_nightly, coverage(off))]
    /// # use pauli_tracker::pauli::{Pauli, PauliDense};
    /// let pauli = PauliDense::new_y();
    /// assert_eq!(pauli.get_z(), true);
    /// # }
    /// ```
    fn get_z(&self) -> bool;

    /// Set whether the Pauli products contains X.
    ///
    /// # Examples
    /// ```
    /// # fn main() { #![cfg_attr(coverage_nightly, coverage(off))]
    /// # use pauli_tracker::pauli::{Pauli, PauliDense};
    /// let mut pauli = PauliDense::new_y();
    /// pauli.set_x(false);
    /// assert_eq!(pauli, Pauli::new_z());
    /// # }
    /// ```
    fn set_x(&mut self, x: bool);

    /// Set whether the Pauli products contains Z.
    ///
    /// # Examples
    /// ```
    /// # fn main() { #![cfg_attr(coverage_nightly, coverage(off))]
    /// # use pauli_tracker::pauli::{Pauli, PauliDense};
    /// let mut pauli = PauliDense::new_y();
    /// pauli.set_z(false);
    /// assert_eq!(pauli, Pauli::new_x());
    /// # }
    /// ```
    fn set_z(&mut self, z: bool);

    /// Translate into the tableau encoding
    fn tableau_encoding(&self) -> u8;
}

