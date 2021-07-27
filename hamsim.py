"""
CSC299 ROP
Application: Hamiltonian simulation

While there are better methods for Hamiltonian simulation, this file is simply intended to show
one possible use of the LCU algorithm.
"""

from typing import Union
import tequila as tq
from lcu_v1 import LCU
from math import log10 as log
from math import factorial


def ham_sim(ancilla: Union[list[Union[str, int]], str, int],
            unitaries: list[tuple[float, tq.QCircuit]],
            error: float, segments: int, time: float) -> tq.QCircuit:
    """..."""
    segment_evol = segmented_ham_sim(ancilla=ancilla, unitaries=unitaries,
                                     error=error, time=time / segments)

    circ = tq.QCircuit()

    for _ in range(segments):
        circ += segment_evol

    return circ


def segmented_ham_sim(ancilla: Union[list[Union[str, int]], str, int],
                      unitaries: list[tuple[float, tq.QCircuit]],
                      error: float, segments: int, time: float) -> tq.QCircuit:
    """..."""
    lcu = LCU(ancilla=ancilla, unitaries=unitaries)
    lcu_ham = lcu.full_circuit

    num_terms = _num_powers(error=error, segments=segments)

    lin_comb_powers = []
    # Compute coefficients and powers of H
    for k in range(num_terms):
        power = tq.QCircuit()
        k_factorial = factorial(k)

        # How to make positive coefficients, instead of complex?
        # coeff = ((-1j * time) ** k) / k_factorial
        coeff = (time ** k) / k_factorial

        for _ in range(k):
            power += lcu_ham

        lin_comb_powers.append((coeff, power))

    segment_lcu = LCU(ancilla=..., unitaries=lin_comb_powers)

    return segment_lcu.full_circuit


def _num_powers(error: float, segments: int) -> int:
    """Return K = (log (segments / error)) / (log (log (segments / error)))"""
    return int((log(segments / error)) / (log(log(segments / error))))
