"""
Module for implementing black-box quantum state preparation based on Sanders, Low, Scherer, Berry.

Implementation using 2's complement bit adder for comparisons.
"""

import tequila as tq
# from typing import Union


def prepare(target_coeff: list[int], ancilla: list) -> tq.QCircuit:
    """Prepare the quantum state given by the coefficient vector target_coeff in the output register ancilla[0],
    using additional ancilla registers ancilla[1] (data), ancilla[2] (ref), ancilla[3] (flag).

    Preconditions:
        - len(ancilla) == 4
        - isinstance(ancilla[0], list)
        - isinstance(ancilla[1], list)
        - isinstance(ancilla[2], list)
        - not isinstance(ancilla[3], list)
    """


def prep_uniform_superposition(reg) -> tq.QCircuit:
    """Prepare the given register in a uniform superposition of computational basis states, assuming all
    qubits in register start in the zero state."""
    circ = tq.QCircuit()
    for anc in reg:
        circ += tq.gates.H(target=anc)
    return circ


def target_proc() -> tq.QCircuit:
    """Process the given coefficient target vector and prepare the amplitudes to be written to data register."""


def write_amp(target_coeff, out_reg, data_reg):
    """Write target """


def compare_op(data_reg, ref_reg, flag):
    """Compare the registers data_reg, ref_reg and store result in flag"""


if __name__ == "__main__":
    ...
