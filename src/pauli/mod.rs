
pub mod pauli;
pub mod pauli_tuple;
pub mod pauli_dense;
//mod pauli_test;
//pub mod pauli_string;

pub use pauli::Pauli;
pub use pauli_tuple::PauliTuple;
pub use pauli_dense::PauliDense;


pub type PauliString = Vec<(usize, PauliTuple)>;


pub mod tableau_encoding {
    /// Code for the identity.
    pub const I: u8 = 0;
    /// Code for the Pauli X gate.
    pub const X: u8 = 2;
    /// Code for the Pauli Y gate.
    pub const Y: u8 = 3;
    /// Code for the Pauli Z gate.
    pub const Z: u8 = 1;
}

