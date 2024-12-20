
use hashbrown::hash_map::DefaultHashBuilder;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
//use thiserror::Error;

//use super::{MissingBit, PauliString, Tracker};
/*use crate::{
    boolean_vector::BooleanVector,
    collection::{Base, Full, Init, IterableBase, Map},
    pauli::{Pauli, PauliStack, PauliTuple},
};*/

//use bitvec::prelude::{BitVec, BitSlice};

use crate::pauli::{Pauli, PauliTuple};

use bitvec::vec::BitVec;
use bitvec::ptr::BitRef;
use bitvec::slice::BitSlice;
use bitvec::field::BitField;

use core::ops::{BitXorAssign, BitOrAssign, BitXor, BitOr};


///
///Retrieves bit slices from the ffstacks
///for each respective qbit
///
#[derive(Debug)]
pub struct FFBitSlice<'a> {
    x: &'a BitSlice,
    z: &'a BitSlice
}

impl<'a> FFBitSlice<'a> {
    pub fn new(x: &'a BitSlice, z: &'a BitSlice) -> Self {
        FFBitSlice { x, z }
    }
}

#[derive(Debug)]
struct FFCopySpace {
    pub copy_a: BitVec,
    pub copy_b: BitVec,
    pub copy_c: BitVec,
    pub copy_d: BitVec,
}

impl FFCopySpace {
    pub fn new(width: usize) -> Self {
        let mut ffc = FFCopySpace {
            copy_a: BitVec::repeat(false, width),
            copy_b: BitVec::repeat(false, width),
            copy_c: BitVec::repeat(false, width),
            copy_d: BitVec::repeat(false, width),
        };
        ffc.copy_a.force_align();
        ffc.copy_b.force_align();
        ffc.copy_c.force_align(); 
        ffc.copy_d.force_align();
        return ffc;
    }
}


#[derive(Debug)]
pub struct FFStacks {
    x: BitVec,
    z: BitVec,
    width: usize,
    copy_space: FFCopySpace,
}

///
/// Each qbit has its own vector slice
/// that can be operated on
///
/// For track_pauli, the jump between locations will
/// likely trigger a ram hit while cz's, cy's and all that will be
/// fairly cheap
///
/// TODO: Allow frame_idx to be the segment instead of width
/// 
///
impl FFStacks {

    pub fn with_capacity(num_keys: usize, width: usize) -> Self {
        let mut ff = FFStacks {
            x: BitVec::repeat(false, num_keys),
            z: BitVec::repeat(false, num_keys),
            width,
            copy_space: FFCopySpace::new(width)
        };
        ff.x.force_align();
        ff.z.force_align();
        return ff;
    }
        
    //TODO: Sice the bitvec is bit index, we want to likely buffer operations
    //into a vector operation
    //TODO: Align the bitvectors on WORD lines, buffer and SIMD it
    pub fn encode(&mut self, qubit: usize, frame_offset: usize, 
        pauli: PauliTuple) {
        self.x.set((qubit * self.width) + frame_offset, pauli.get_x());
        self.z.set((qubit * self.width) + frame_offset, pauli.get_z());
    }
   
    

    ///Performs an xor_inplace on a qubit, z xor x
    ///TODO: Validate this operation
    ///
    pub fn s(&mut self, qubit: usize, frame_idx: usize) {
        //xor 1 bit
        let offset = (qubit * self.width) + frame_idx;
        let x_ref = &mut self.x[offset..offset+1];
        x_ref.bitxor_assign(&self.z[offset..offset+1])    
    }

    /// Conjugate the PauliStack with the Hadamard gate ignoring phases.
    /// TODO: Validate this operation
    pub fn h(&mut self, qubit: usize, frame_idx: usize) {
        //self.swap_space.copy_within(self.x[qubit..], self.width);

        let offset = (qubit * self.width) + frame_idx;
        //eprintln!("off:{}, qbit:{}, w:{}, f:{}", offset, qubit, self.width, frame_idx);
        let z_temp = *self.z.get(offset).unwrap();
        let x_temp = *self.x.get(offset).unwrap();
        self.z.set(offset, x_temp);
        self.x.set(offset, z_temp);
    }
    pub fn sh(&mut self, qubit: usize, frame_idx: usize) {
        self.h(qubit, frame_idx);
        self.s(qubit, frame_idx);
    }

    pub fn hs(&mut self, qubit: usize, frame_idx: usize) {
        self.s(qubit, frame_idx);
        self.h(qubit, frame_idx);
    }

    pub fn shs(&mut self, qubit: usize, frame_idx: usize) {

        let offset = (qubit * self.width) + frame_idx;
        let z_ref = &mut self.z[offset..offset+1];
        z_ref.bitxor_assign(&self.x[offset..offset+1])    

    }

    pub fn get_slice_pair(&mut self, bit1: usize, bit2: usize) 
        -> (FFBitSlice, FFBitSlice) {
        let offset_a = bit1 * self.width;
        let offset_b = bit2 * self.width;
        let bit_1_slice = FFBitSlice::new(
                &self.x[offset_a..offset_a+self.width],
                &self.z[offset_a..offset_a+self.width]
            );
        let bit_2_slice = FFBitSlice::new(
                &self.x[offset_b..offset_a+self.width],
                &self.z[offset_b..offset_a+self.width]
            );

        return (bit_1_slice, bit_2_slice);
    }

    pub fn sum_up(&self, qubit: usize, filter: &[bool]) -> PauliTuple {
        unimplemented!()
    }

    ///
    /// Gets the offset of the bitvec as a slice using the qubit
    /// and internally the width of the vector as given on construction
    ///
    pub fn bitvecs_slice_of_qubit_ref(&self, qubit: usize) {
        unimplemented!()
    }

    pub fn cx(&mut self, control_bit: usize, target_bit: usize) {

        let width = self.width;
        let control_loc = control_bit * width;
        let target_loc = target_bit * width;
        let control_end = control_loc + width;
        let target_end = target_loc + width;

        //eprintln!("c:{} t:{} w:{} cloc:{} tloc:{}", 
        //    control_bit, target_bit, width, control_loc, target_loc);
        //TODO: Consider re-use of memory instead of new bitvec allocations
        //Pref Safe version, TODO: Reduce the cost of an extra copy
        //let tx_mut_copy = &mut self.x[target_loc..target_end].to_bitvec();
        //let cz_mut_copy = &mut self.x[control_loc..control_end].to_bitvec();

        self.copy_space.copy_a.copy_from_bitslice(&self.x[target_loc..target_end]);
        self.copy_space.copy_b.copy_from_bitslice(&self.x[control_loc..control_end]);

        self.copy_space.copy_a.bitxor_assign(&self.x[control_loc..control_end]);
        self.copy_space.copy_b.bitxor_assign(&self.z[target_loc..target_end]);

        self.x[target_loc..target_end].copy_from_bitslice(self.copy_space.copy_a.as_bitslice());
        self.z[control_loc..control_end].copy_from_bitslice(self.copy_space.copy_a.as_bitslice());
            
    }

    pub fn cz(&mut self, bit_a: usize, bit_b: usize) {
               
        let width = self.width;
        let a_loc = bit_a * width;
        let b_loc = bit_b * width;

        let b_x = &mut self.x[b_loc..width].to_bitvec();
        let a_x = &mut self.x[b_loc..width].to_bitvec();


        self.z[a_loc..width].copy_from_bitslice(b_x.as_bitslice());
        self.z[b_loc..width].copy_from_bitslice(a_x.as_bitslice());

    }

    pub fn cy(&mut self, control_bit: usize, target_bit: usize) {
        let width = self.width;
        let control_loc = control_bit * width;
        let target_loc = target_bit * width;

        //TODO: Consider reusing memory
        let tx_mut_copy = &mut self.x[target_loc..width].to_bitvec();
        let tz_mut_copy = &mut self.z[target_loc..width].to_bitvec();
        
        let cz_mut_copy = &mut self.z[control_loc..width].to_bitvec();
        
        cz_mut_copy.bitxor_assign(&self.z[target_loc..width]);
        cz_mut_copy.bitxor_assign(&self.x[target_loc..width]);
        tz_mut_copy.bitxor_assign(&self.x[control_loc..width]);
        tx_mut_copy.bitxor_assign(&self.x[control_loc..width]);

        
        self.x[target_loc..width].copy_from_bitslice(tx_mut_copy.as_bitslice());
        self.z[target_loc..width].copy_from_bitslice(tz_mut_copy.as_bitslice());
        self.z[control_loc..width].copy_from_bitslice(cz_mut_copy.as_bitslice());

    }

    pub fn swap(&mut self, bit_a: usize, bit_b: usize) {
        
        let width = self.width;
        let control_loc = bit_a * width;
        let target_loc = bit_b * width;

        
        self.copy_space.copy_a.copy_from_bitslice(& self.x[control_loc..width]);
        self.copy_space.copy_b.copy_from_bitslice(& self.z[control_loc..width]);
        self.copy_space.copy_c.copy_from_bitslice(& self.x[target_loc..width]);
        self.copy_space.copy_d.copy_from_bitslice(& self.z[target_loc..width]);

        self.x[target_loc..width].copy_from_bitslice(
            self.copy_space.copy_a.as_bitslice());
        self.z[target_loc..width].copy_from_bitslice(
            self.copy_space.copy_b.as_bitslice());
        self.x[control_loc..width].copy_from_bitslice(
            self.copy_space.copy_c.as_bitslice());
        self.z[control_loc..width].copy_from_bitslice(
            self.copy_space.copy_d.as_bitslice());
    
    }

    pub fn iswap(&mut self, bit_a: usize, bit_b: usize) {
        let width = self.width;
        let control_loc = bit_a * width;
        let target_loc = bit_b * width;

        
        self.copy_space.copy_a.copy_from_bitslice(& self.x[control_loc..width]);
        self.copy_space.copy_b.copy_from_bitslice(& self.z[control_loc..width]);
        self.copy_space.copy_c.copy_from_bitslice(& self.x[target_loc..width]);
        self.copy_space.copy_d.copy_from_bitslice(& self.z[target_loc..width]);

        //TODO: Validate the results
        self.copy_space.copy_d.bitxor_assign(self.copy_space.copy_c.as_bitslice());
        self.copy_space.copy_d.bitxor_assign(self.copy_space.copy_a.as_bitslice());
        self.copy_space.copy_b.bitxor_assign(self.copy_space.copy_c.as_bitslice());
        self.copy_space.copy_b.bitxor_assign(self.copy_space.copy_a.as_bitslice());

        self.x[target_loc..width].copy_from_bitslice(
            self.copy_space.copy_a.as_bitslice());
        self.z[target_loc..width].copy_from_bitslice(
            self.copy_space.copy_b.as_bitslice());
        self.x[control_loc..width].copy_from_bitslice(
            self.copy_space.copy_c.as_bitslice());
        self.z[control_loc..width].copy_from_bitslice(
            self.copy_space.copy_d.as_bitslice());

    }



    pub fn zcx(&mut self, control_bit: usize, target_bit: usize) {
        let width = self.width;
        let control_loc = control_bit * width;
        let target_loc = target_bit * width;

        
        self.copy_space.copy_a.copy_from_bitslice(& self.x[target_loc..width]);
        self.copy_space.copy_b.copy_from_bitslice(& self.x[control_loc..width]);

        //TODO: Validate the results
        self.copy_space.copy_d.bitxor_assign(self.copy_space.copy_c.as_bitslice());
        self.copy_space.copy_d.bitxor_assign(self.copy_space.copy_a.as_bitslice());

        self.x[target_loc..width].copy_from_bitslice(
            self.copy_space.copy_a.as_bitslice());
        self.x[control_loc..width].copy_from_bitslice(
            self.copy_space.copy_b.as_bitslice());
    }
    
    pub fn zcy(&mut self, control_bit: usize, target_bit: usize) {
        let width = self.width;
        let control_loc = control_bit * width;
        let target_loc = target_bit * width;

        
        self.copy_space.copy_a.copy_from_bitslice(& self.x[control_loc..width]);
        self.copy_space.copy_b.copy_from_bitslice(& self.z[target_loc..width]);
        self.copy_space.copy_c.copy_from_bitslice(& self.x[target_loc..width]);

        //TODO: Validate the results
        self.copy_space.copy_a.bitxor_assign(&self.z[target_loc..width]);
        self.copy_space.copy_a.bitxor_assign(&self.x[target_loc..width]);
        self.copy_space.copy_b.bitxor_assign(&self.z[control_loc..width]);
        self.copy_space.copy_c.bitxor_assign(&self.z[control_loc..width]);

        self.x[target_loc..width].copy_from_bitslice(
            self.copy_space.copy_a.as_bitslice());
        self.z[target_loc..width].copy_from_bitslice(
            self.copy_space.copy_b.as_bitslice());
        self.x[control_loc..width].copy_from_bitslice(
            self.copy_space.copy_c.as_bitslice());
    }

}

#[derive(Debug)]
pub struct FastFrames {
    ffs: FFStacks,
    frames_num: usize,
    nqbits: usize,
    //TODO: Supply an inverse mapping


}

impl FastFrames {

    pub fn new(num_keys: usize) -> Self {
       
        FastFrames {
            //stacks: vec![PauliTuple::new_i(); num_keys*num_keys],
            ffs: FFStacks::with_capacity(
                        num_keys*num_keys,
                        num_keys),
            frames_num: 0,
            nqbits: num_keys,
        }

    }

    pub fn init(nqbits: usize) -> Self {
        FastFrames::new(nqbits)
    }

    pub fn as_storage(&self) -> &Self {
        &self
    }

    pub fn len(&self) -> usize {
        return self.nqbits
    }

    pub fn frames_num(&self) -> usize {
        self.frames_num
    }

    // TODO: This should be known on construction
    // Transform this later on
    //pub fn bulk_copy_pauli_i(&self) -> Vec<PauliTuple> {
    //    self.cloneable_stack.clone()
    //}  
    //TODO: Validate the tracking of pauli states
    pub fn track_pauli_ff(&mut self, qubit: usize, pauli: PauliTuple) {
        if self.nqbits == 0 {
            return;
        } 
        self.ffs.encode(qubit, self.frames_num, pauli);
        self.frames_num += 1;
    }

    pub fn track_x(&mut self, qubit: usize, pauli: PauliTuple) {
        self.track_pauli_ff(qubit, pauli);
    }

    pub fn track_y(&mut self, qubit: usize, pauli: PauliTuple) {
        self.track_pauli_ff(qubit, pauli);
    }

    pub fn track_z(&mut self, qubit: usize, pauli: PauliTuple) {
        self.track_pauli_ff(qubit, pauli);
    }

    pub fn transpose<T>(len: usize) {
        unimplemented!()
    }

    // Fixing the PauliStack operations that will operate only on sets of PauliTuples
    // Make sure the mapping is correct
    //TODO: Fix the call that PauliStack does
    //Get the PauliStack?
    pub fn s(&mut self, bit: usize) {     
        self.ffs.s(bit, self.frames_num);
    }
    
    pub fn h(&mut self, bit: usize) {     
        self.ffs.h(bit, self.frames_num);
    }

    pub fn sh(&mut self, bit: usize) {     
        self.ffs.sh(bit, self.frames_num);
    }

    pub fn hs(&mut self, bit: usize) {     
        self.ffs.hs(bit, self.frames_num);
    }
    
    pub fn shs(&mut self, bit: usize) {
        self.ffs.shs(bit, self.frames_num);
    } 
 
    pub fn cz(&mut self, bit_a: usize, bit_b: usize) {
        self.ffs.cz(bit_a, bit_b);
    }

    pub fn cx(&mut self, control: usize, target: usize) {
        self.ffs.cx(control, target);
    }

    pub fn cy(&mut self, control_bit: usize, target_bit: usize) {
        self.ffs.cy(control_bit, target_bit);
    }

}

