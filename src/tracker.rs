use crate::pauli::{Pauli, PauliString, PauliTuple};
//use crate::clifford;
//track_pauli!((track_x, X), (track_y, Y), (track_z, Z),);

macro_rules! remove {
    ($((
        $name:ident,
        $correction:literal
    ),)*) => {$(
        /// "Remove" the
        #[doc=$correction]
        /// Pauli stack from the qu`bit`.
        #[allow(unused_variables)]
        fn $name(&mut self, bit: usize) {
            panic!(
                "the default implementation exists only to not make a major breaking
                change in this trait"
            );
        }
    )*}
}

macro_rules! movements {
    ($((
        $name:ident,
        $from_doc:literal,
        $to_doc:literal
    ),)*) => {$(
        /// "Move" the
        #[doc=$from_doc]
        /// Pauli stack from the `origin` qubit to the `destination` qubit, transforming
        /// it to an
        #[doc=$to_doc]
        /// stack.
        fn $name(&mut self, source: usize, destination: usize);
    )*}
}


macro_rules! track_pauli {
    ($(($name:ident, $gate:ident),)*) => {$(
        /// Track a new frame consisting of the Pauli
        #[doc = stringify!($gate)]
        /// at qu`bit`.
        fn $name(&mut self, bit: usize) {
            //self.track_pauli(bit, Self::Pauli::$gate );
            self.track_pauli(bit, Pauli::$gate)
        }
    )*};
}

pub struct MissingBit(pub usize);

pub trait Tracker {
    /// The storage type used to store the tracked Paulis for each qubit, e.g.,
    /// [PauliStack](crate::pauli::PauliStack) for the [Frames](frames::Frames) tracker or
    /// just a simple [Pauli] for the [Live](live::Live) tracker (in this case it's a
    /// stack with one element ...).
    type Stack;

    /// The type of Pauli representation used for operations like
    /// [track_pauli](Self::track_pauli). It is usally the type that is the most
    /// compatible with [Self::Stack].

    /// Insert a new qu`bit` into the tracker. If the qu`bit` is already present, the old
    /// value is overwritten and returned.
    fn new_qubit(&mut self, bit: usize) -> Option<Self::Stack>;

    /// Track a new frame consisting of the Pauli gate `pauli` at qu`bit`.
    ///
    /// If qu`bit` is not tracked, the method does not error, but simply tracks an empty
    /// frame.
    fn track_pauli(&mut self, bit: usize, pauli: PauliTuple);

    /// Track a new frame including multiple Pauli gates, i.e., i.e., do
    /// [Tracker::track_pauli] for multiple Paulis but all within the same frame.
    fn track_pauli_string(&mut self, string: PauliString);

    track_pauli!((track_x, X), (track_y, Y), (track_z, Z),);

    fn h(&mut self, bit: usize);
    fn s(&mut self, bit: usize);
    fn cz(&mut self, bit_a: usize, bit_b: usize);
    //clifford::trait_gates!();

    movements!(
        (move_x_to_x, "X", "X"),
        (move_x_to_z, "X", "Z"),
        (move_z_to_x, "Z", "X"),
        (move_z_to_z, "Z", "Z"),
    );

    remove!((remove_x, "X"), (remove_z, "Z"),);

    /// Remove the Pauli stack on qu`bit`, if it is present.
    fn measure(&mut self, bit: usize) -> Result<Self::Stack, MissingBit>;
}


#[cfg(test)]
pub mod tests {
    use super::*;
    pub mod utils {
        use super::*;
        use crate::{ff::fastframe::FastFrames, pauli::{PauliDense, PauliString}};

        // when we update the results here and use this module in the test of the tracker
        // implementors, the type system ensures that we test all gates/actions

        //                 name for debugging, expected results
        pub type SingleResults = (&'static str, [u8; 4]);
        pub type DoubleResults = (&'static str, [(u8, u8); 16]);
        pub type SingleAction<T> = fn(&mut T, usize);
        pub type DoubleAction<T> = fn(&mut T, usize, usize);

        // the following expected results are proven in ./docs/conjugation_rules.pdf
        //
        // instead of writing out all the SingleResults and DoubleResults, we make use
        // of homomorphy and just define the results on a basis
        //
        // the encoding is according to crate::pauli::tableau_encoding, i.e., 0=I, 2=X,
        // 3=Y, 1=Z

        pub const N_SINGLES: usize = 20;
        #[rustfmt::skip]
        const SINGLE_GENERATORS: [(&str, [u8; 2]); N_SINGLES] =
            // (name, result: [conjugate Z, conjugate X])
            [
                ("I",    [1, 2]),
                ("X",    [1, 2]),
                ("Y",    [1, 2]),
                ("Z",    [1, 2]),
                ("S",    [1, 3]),
                ("SDG",  [1, 3]),
                ("SZ",   [1, 3]),
                ("SZDG", [1, 3]),
                ("H_xy", [1, 3]),
                ("H",    [2, 1]),
                ("SY",   [2, 1]),
                ("SYDG", [2, 1]),
                ("SH",   [3, 1]),
                ("HS",   [2, 3]),
                ("SHS",  [3, 2]),
                ("SX",   [3, 2]),
                ("SXDG", [3, 2]),
                ("H_yz", [3, 2]),
                // these here are not conjugations with unitary operators, however it
                // still works, because the operation is a homomorphism
                ("remove_z", [0, 2]),
                ("remove_x", [1, 0]),
            ];

        /*macro_rules! single_actions {
            ($tracker:ty) => {
                [
                    <$tracker>::id,
                    <$tracker>::x,
                    <$tracker>::y,
                    <$tracker>::z,
                    <$tracker>::s,
                    <$tracker>::sdg,
                    <$tracker>::sz,
                    <$tracker>::szdg,
                    <$tracker>::hxy,
                    <$tracker>::h,
                    <$tracker>::sy,
                    <$tracker>::sydg,
                    <$tracker>::sh,
                    <$tracker>::hs,
                    <$tracker>::shs,
                    <$tracker>::sx,
                    <$tracker>::sxdg,
                    <$tracker>::hyz,
                    <$tracker>::remove_z,
                    <$tracker>::remove_x,
                ]
            };
        }        pub(crate) use single_actions;*/

        pub const N_DOUBLES: usize = 13;
        #[rustfmt::skip]
        const DOUBLE_GENERATORS: [(&str, [(u8, u8); 4]); N_DOUBLES] = [
            //+ (name, result: [conjugate Z1, conjugate Z2, conjugate X1, conjugate X2])
            // the left tuple entry of the results belongs to the second qubit (q0) in the
            // function call and the right entry to the first one (q1), i.e., q1 controls
            // q0
            ("cz",          [(1, 0), (0, 1), (2, 1), (1, 2)]),
            ("cx",          [(1, 1), (0, 1), (2, 0), (2, 2)]),
            ("cy",          [(1, 1), (0, 1), (2, 1), (3, 2)]),
            ("swap",        [(0, 1), (1, 0), (0, 2), (2, 0)]),
            ("zcz",         [(1, 0), (1, 1), (2, 2), (0, 2)]),
            ("zcx",         [(1, 2), (2, 1), (2, 0), (0, 2)]),
            ("zcy",         [(1, 2), (3, 1), (2, 2), (0, 2)]),
            ("iswap",       [(0, 1), (1, 0), (1, 3), (3, 1)]),
            ("iswapdg",     [(0, 1), (1, 0), (1, 3), (3, 1)]),
            // cf comment above for remove_*
            ("move_x_to_x", [(1, 0), (0, 1), (2, 0), (2, 0)]),
            ("move_x_to_z", [(1, 0), (0, 1), (2, 0), (1, 0)]),
            ("move_z_to_x", [(1, 0), (2, 0), (2, 0), (0, 2)]),
            ("move_z_to_z", [(1, 0), (1, 0), (2, 0), (0, 2)]),
        ];

        /*macro_rules! double_actions {
            ($tracker:ty) => {
                [
                    <$tracker>::cz,
                    <$tracker>::cx,
                    <$tracker>::cy,
                    <$tracker>::swap,
                    <$tracker>::zcz,
                    <$tracker>::zcx,
                    <$tracker>::zcy,
                    <$tracker>::iswap,
                    <$tracker>::iswapdg,
                    <$tracker>::move_x_to_x,
                    <$tracker>::move_x_to_z,
                    <$tracker>::move_z_to_x,
                    <$tracker>::move_z_to_z,
                ]
            };
        }
        pub(crate) use double_actions;*/
        // Old version:
        //
        //pub fn single_check<T, R>(runner: R, actions: [SingleAction<T>; N_SINGLES])
        //
        #[cfg_attr(coverage_nightly, coverage(off))]
        pub fn single_check<R>(runner: R, actions: [SingleAction<FastFrames>; N_SINGLES])
        where
            R: Fn(SingleAction<FastFrames>, SingleResults),
        {
            for (action, result_generator) 
                in actions.into_iter().zip(SINGLE_GENERATORS) {

                let mut results = [0; 4];
                for (i, r) in results.iter_mut().enumerate() {
                    *r = (if (i & 1) > 0 {
                        result_generator.1[0]
                    } else {
                        0
                    }) ^ (if (i & 2) > 0 {
                        result_generator.1[1]
                    } else {
                        0
                    })
                }
                (runner)(action, (result_generator.0, results))
            }
        }

        //Old version
        //pub fn double_check<T, R>(runner: R, actions: [DoubleAction<T>; N_DOUBLES])
        
        #[cfg_attr(coverage_nightly, coverage(off))]
        pub fn double_check<R>(runner: R, actions: [DoubleAction<FastFrames>; N_DOUBLES])
        where
            R: Fn(DoubleAction<FastFrames>, DoubleResults),
        {
            for (action, result_generator) 
                in actions.into_iter().zip(DOUBLE_GENERATORS) {

                let mut results = [(0, 0); 16];
                for (i, r) in (0..).zip(results.iter_mut()) {
                    // cf. the masks below in double_init
                    let a = if (i & 1) > 0 {
                        result_generator.1[0]
                    } else {
                        (0, 0)
                    };
                    let b = if (i & 2) > 0 {
                        result_generator.1[2]
                    } else {
                        (0, 0)
                    };
                    let c = if (i & 4) > 0 {
                        result_generator.1[1]
                    } else {
                        (0, 0)
                    };
                    let d = if (i & 8) > 0 {
                        result_generator.1[3]
                    } else {
                        (0, 0)
                    };
                    *r = (a.0 ^ b.0 ^ c.0 ^ d.0, a.1 ^ b.1 ^ c.1 ^ d.1)
                }
                (runner)(action, (result_generator.0, results))
            }
        }

        #[cfg_attr(coverage_nightly, coverage(off))]
        pub fn single_init<T: From<PauliTuple>>(input: u8) -> PauliString {
            let pd: PauliDense = PauliDense::try_from(input).unwrap().into();
            let pt = PauliTuple::new_product(pd.get_z(), 
                pd.get_x());
            vec![(0, pt)]
        }

        #[cfg_attr(coverage_nightly, coverage(off))]
        pub fn double_init<T: From<PauliTuple>>(input: u8) -> PauliString {
            // masks to decode p in 0..16 into two paulis and vice versa
            const SECOND: u8 = 12; // = 1100
            const FIRST: u8 = 3; // = 0011
            const SECOND_SHIFT: u8 = 2;
            let pd1: PauliDense = 
                 PauliDense::try_from((input & SECOND) >> SECOND_SHIFT)
                 .unwrap()
                 .into();
            let pd2: PauliDense = 
                 PauliDense::try_from(input & FIRST)
                 .unwrap()
                 .into();
            let pt1 = PauliTuple::new_product(
                pd1.get_z(), 
                pd1.get_x()
            );
            let pt2 = PauliTuple::new_product(
                pd2.get_z(), 
                pd2.get_x()
            );
            vec![
                (1, pt1),
                (0, pt2),
            ]
        }

        #[cfg_attr(coverage_nightly, coverage(off))]
        pub fn double_output<T: Into<PauliDense>>(
            frame: impl IntoIterator<Item = (usize, T)>,
        ) -> (u8, u8) {
            let mut output = [0, 0];
            for (i, p) in frame {
                output[i] = p.into().storage()
            }
            (output[0], output[1])
        }
    }

    mod defaults {
        use std::collections::HashMap;

        use coverage_helper::test;

        use super::{
            super::*,
            utils::{DoubleAction, DoubleResults, N_DOUBLES},
        };
        use crate::{
            cols::Base, ff::fastframe::FastFrames, pauli::{PauliDense, PauliString}
        };

        #[derive(Debug)]
        struct DefaultTester {
            paulis: HashMap<usize, PauliTuple>,
            skip_it: bool,
        }

        impl DefaultTester {
            fn init(n: usize) -> Self {
                Self {
                    paulis: HashMap::from_iter((0..n).map(|i| (i, PauliTuple::I))),
                    skip_it: false,
                }
            }
        }

        impl Tracker for DefaultTester {
            type Stack = PauliTuple;

            fn new_qubit(&mut self, bit: usize) -> Option<Self::Stack> {
                self.paulis.insert(bit, PauliTuple::I)
            }
            fn track_pauli(&mut self, _: usize, _: PauliTuple) {
                todo!()
            }
            fn track_pauli_string(&mut self, string: PauliString) {
                for (bit, pauli) in string {
                    if let Some(p) = self.paulis.get_mut(&bit) {
                        p.multiply(pauli)
                    }
                }
            }

            fn h(&mut self, bit: usize) {
                self.paulis.get_mut(&bit).unwrap().h()
            }
            fn s(&mut self, bit: usize) {
                self.paulis.get_mut(&bit).unwrap().s()
            }
            fn cz(&mut self, bit_a: usize, bit_b: usize) {
                unimplemented!()
                /*let (a, b) = 
                    self.paulis.get_two_mut(bit_a, bit_b).unwrap();
                a.zpx(b);
                b.zpx(a);*/
            }

            fn move_x_to_x(&mut self, _: usize, _: usize) {
                self.skip_it = true
            }
            fn move_x_to_z(&mut self, _: usize, _: usize) {
                self.skip_it = true
            }
            fn move_z_to_x(&mut self, _: usize, _: usize) {
                self.skip_it = true
            }
            fn move_z_to_z(&mut self, _: usize, _: usize) {
                self.skip_it = true
            }
            fn remove_x(&mut self, _: usize) {
                self.skip_it = true
            }
            fn remove_z(&mut self, _: usize) {
                self.skip_it = true
            }

            fn measure(&mut self, _: usize) -> Result<Self::Stack, MissingBit> {
                todo!()
            }
        }

        use super::*;
        use crate::tracker::tests::utils::{N_SINGLES, SingleAction, SingleResults};

        type ActionS = SingleAction<FastFrames>;
        type ActionD = DoubleAction<FastFrames>;

        #[cfg_attr(coverage_nightly, coverage(off))]
        fn single_runner(action: ActionS, result: SingleResults) {
            for (input, check) in (0u8..).zip(result.1) {
                let mut tracker = FastFrames::init(2);
                tracker.track_pauli_string(utils::single_init(input));
                (action)(&mut tracker, 0);
                /*if tracker.skip_it {
                    tracker.skip_it = false;
                    return;
                }*/
                let computed = tracker.paulis.get(&0).unwrap()
                    .storage();
                assert_eq!(
                    computed, check,
                    "gate: {}, input: {}, expected: {}, computed: {}",
                    result.0, input, check, computed
                );
            }
        }

        #[test]
        fn single_actions() {
            let actions: [ActionS; N_SINGLES] = utils::single_actions!(DefaultTester);
            utils::single_check(single_runner, actions);
        }

        #[cfg_attr(coverage_nightly, coverage(off))]
        fn double_runner(action: ActionD, result: DoubleResults) {
            for (input, check) in (0u8..).zip(result.1) {
                let mut tracker = FastFrames::init(2);
                tracker.track_pauli_string(utils::double_init(input));
                (action)(&mut tracker, 1, 0);
                /*if tracker.skip_it {
                    tracker.skip_it = false;
                    return;
                }*/
                let computed = utils::double_output(tracker.paulis);
                assert_eq!(
                    computed, check,
                    "{}, {}, {:?}, {:?}",
                    result.0, input, check, computed
                );
            }
        }

        #[test]
        fn double_actions() {
            let actions: [ActionD; N_DOUBLES] = utils::double_actions!(DefaultTester);
            utils::double_check(double_runner, actions);
        }
    }
}
