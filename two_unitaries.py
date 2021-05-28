"""
CSC299 ROP: Simple implementation
"""

import tequila as tq
import copy
from typing import Iterable
import numpy as np


def prepare_2unitary(ancillary, sum_of_unitaries: list[tuple[float, tq.QCircuit]]) -> tq.QCircuit:
    """
    Prepare operator, when the Hamiltonian can be expressed as the linear combination of two
    unitary operators.
    Requires only one ancillary qubit.

    Preconditions:
        - ...
    """
    alpha_0, alpha_1 = sum_of_unitaries[0][0], sum_of_unitaries[1][0]

    theta = -2 * np.arcsin(np.sqrt(alpha_1 / (alpha_0 + alpha_1)))

    return tq.gates.Ry(target=ancillary, angle=theta)


def select_2unitary(ancillary, unitary_0: tq.QCircuit, unitary_1: tq.QCircuit) -> tq.QCircuit:
    """
    Select operator, when the Hamiltonian can be expressed as the linear combination of two
    unitary operators.
    Requires only one ancillary qubit.

    Returns ...
    """
    impl_1 = control_unitary(ancilla=ancillary, unitary=unitary_1)

    x_gate = tq.gates.X(target=ancillary)
    control = control_unitary(ancilla=ancillary, unitary=unitary_0)

    impl_0 = x_gate + control + x_gate

    return impl_1 + impl_0


def control_unitary(ancilla, unitary: tq.QCircuit) -> tq.QCircuit:
    """Return controlled version of unitary

    SHOULD NOT mutate unitary

    Preconditions:
        - ancilla and unitary cannot have any common qubits
    """
    gates = unitary.gates
    cgates = []
    for gate in gates:
        cgate = copy.deepcopy(gate)
        if isinstance(ancilla, Iterable):
            control_lst = list(cgate.control) + list(ancilla)
        else:
            control_lst = list(cgate.control) + [ancilla]
        cgate._control = tuple(control_lst)
        cgate.finalize()
        cgates.append(cgate)

    return tq.QCircuit(gates=cgates)


def example_function() -> tq.QCircuit:
    """Test example for LCU with two unitary matrices"""
    identity = tq.gates.X(1) + tq.gates.X(1)
    unitaries = [(0.5, identity), (0.5, tq.gates.Z(1))]
    return algorithm_2unitary(unitaries=unitaries)


def lcu_2unitary(nonunitary: tq.numpy.array, sum_of_unitaries: list[tuple[float, tq.QCircuit]]) \
        -> tq.QCircuit:
    """Linear combinations of unitaries for the case with two unitaries.

    Preconditions:
        - len(unitaries) == 2
    """
    if check_weighted_sum(nonunitary, sum_of_unitaries):
        return algorithm_2unitary(sum_of_unitaries)
    else:
        raise ValueError


def algorithm_2unitary(unitaries: list[tuple[float, tq.QCircuit]]) -> tq.QCircuit:
    """The main part.

    TODO: write proper docstring
    """
    prepare = prepare_2unitary(0, unitaries)
    circ = prepare + select_2unitary(0, unitaries[0][1], unitaries[1][1]) + prepare.dagger()
    return circ


def check_weighted_sum(nonunitary: tq.numpy.array, unitaries: list[tuple[complex, tq.QCircuit]]) \
        -> bool:
    """
    Return whether the matrix given is indeed the weighted sum of the unitaries.
    """
    ...
