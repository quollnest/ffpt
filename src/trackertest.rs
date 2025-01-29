
#[cfg(test)]
mod tests {
    use coverage_helper::test;

    use super::*;

    // we only check the basic functionality here, more complicated circuits are tested
    // in [super::circuit] to test the tracker and the circuit at once

    mod action_definition_check {
        //use super::{super::*, test, *};
        use crate::{
            //collection::BufferedVector,
            ff::fastframe::FastFrames, pauli::PauliDense, tracker::tests::utils::{
                self, DoubleAction, DoubleResults, SingleAction, SingleResults, N_DOUBLES, N_SINGLES
            }
        };

        // maybe todo: in the following functions there's a pattern behind how we encode
        // one-qubit and two-qubit actions, it's like a "TwoBitVec"; one could probably
        // implement that in connection with [Pauli]

        type ThisTracker = FastFrames;
        // type ThisTracker =
        //     Frames<Vector<crate::boolean_vector::bitvec_simd::SimdBitVec>>;

        #[test]
        fn single() {
            type Action = SingleAction<ThisTracker>;

            const ACTIONS: [Action; N_SINGLES] = utils::single_actions!(ThisTracker);

            #[cfg_attr(coverage_nightly, coverage(off))]
            fn runner(action: Action, result: SingleResults) {
                let mut tracker: ThisTracker = FastFrames::init(1);
                for input in (0..4).rev() {
                    tracker.track_pauli_string(vec![(
                        0,
                        PauliDense::try_from(input).unwrap().into(),
                    )]);
                }
                (action)(&mut tracker, 0);
                for (input, check) in (0u8..).zip(result.1) {
                    let computed = tracker
                        //.pop_frame::<PauliDense>()
                        //.unwrap()
                        .first()
                        .unwrap()
                        .1
                        .storage();
                    assert_eq!(
                        computed, check,
                        "gate: {}, input: {}, expected: {}, computed: {}",
                        result.0, input, check, computed
                    );
                }
            }

            utils::single_check(runner, ACTIONS)
        }

        #[test]
        fn double() {
            type Action = DoubleAction<ThisTracker>;

            const ACTIONS: [Action; N_DOUBLES] = utils::double_actions!(ThisTracker);

            #[cfg_attr(coverage_nightly, coverage(off))]
            fn runner(action: Action, result: DoubleResults) {
                let mut tracker: ThisTracker = FastFrames::init(2);
                for pauli in (0..16).rev() {
                    tracker.track_pauli_string(utils::double_init(pauli));
                }
                (action)(&mut tracker, 1, 0);
                for (input, check) in (0u8..).zip(result.1) {
                    let computed =
                        utils::double_output(tracker.pop_frame::<PauliDense>().unwrap());
                    assert_eq!(
                        computed, check,
                        "gate: {}, input: {}, expected: {:?}, computed: {:?}",
                        result.0, input, check, computed
                    );
                }
            }

            utils::double_check(runner, ACTIONS);
        }
    }
}
