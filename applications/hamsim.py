"""
Automatized non-unitary embeddings: Time-independent Hamiltonian simulation

While there are better methods for Hamiltonian simulation, this file is simply intended to show
one possible use of the LCU algorithm.
"""

from typing import Union
import tequila as tq
from lcu_v1 import LCU
from math import log10 as log
from math import factorial


def ham_power(ancilla: list[Union[list[Union[str, int]]]],
              unitaries: list[tuple[float, tq.QCircuit]], power: int) -> tq.QCircuit:
    """Return the circuit corresponding to H^k, where H is the Hamiltonian corresponding to
    the linear combination expressed in unitaries, and k = power.

    Preconditions:
        - len(ancilla) == power
        - all(2 ** (len(anc) - 1) < len(unitaries) <= 2 ** len(anc) for anc in ancilla)
        - ancilla has no repeated qubits
    """
    circ = tq.QCircuit()
    anc0 = ancilla[0]
    lcu = LCU(anc0, unitaries).full_circuit

    for i in range(power):
        anc = ancilla[i]
        qubit_map = {anc0[k]: anc[k] for k in range(len(anc0))}
        lcu = lcu.map_qubits(qubit_map)
        circ += lcu

    return circ


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
