"""
CSC299 ROP: Simple implementation
"""

import tequila as tq
import copy
from typing import Iterable
from numpy import arcsin, sqrt, floor, pi


def prepare_2unitary(ancillary, sum_of_unitaries: list[tuple[float, tq.QCircuit]]) -> tq.QCircuit:
    """
    Prepare operator, when the Hamiltonian can be expressed as the linear combination of two
    unitary operators.
    Requires only one ancillary qubit.

    Preconditions:
        - ...
    """
    alpha_0, alpha_1 = sum_of_unitaries[0][0], sum_of_unitaries[1][0]

    theta = -2 * arcsin(sqrt(alpha_1 / (alpha_0 + alpha_1)))

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


def check_weighted_sum(nonunitary: tq.numpy.array, unitaries: list[tuple[float, tq.QCircuit]]) \
        -> bool:
    """
    Return whether the matrix given is indeed the weighted sum of the unitaries.

    TODO
    """
    ...


def amp_amp(unitaries: list[tuple[float, tq.QCircuit]], walk_op: tq.QCircuit, ancilla) \
        -> tq.QCircuit:
    """Amplitude amplification procedure obtained by repeating the amplitude amplification
    step for a total of s times where s is the result of function _num_iter()
    """
    amplification_operator = amp_amp_op(walk_op, ancilla)
    s = _num_iter(unitaries)

    sum_of_steps = tq.QCircuit()
    for _ in range(s):
        sum_of_steps += amplification_operator

    return walk_op + sum_of_steps


def amp_amp_op(walk_op: tq.QCircuit, ancilla) -> tq.QCircuit:
    """Return WRW.dagger()R,
     where R is the reflect operator returned by the func reflect_operator"""
    anc_qubits = ancilla if isinstance(ancilla, list) else [ancilla]
    state_qubits = [qubit for qubit in walk_op.qubits if qubit not in anc_qubits]

    reflect = reflect_operator(state_qubits=state_qubits, ancilla=ancilla)

    return reflect + walk_op.dagger() + reflect + walk_op


def reflect_operator(state_qubits, ancilla) -> tq.QCircuit:
    """
    Return the reflection operator R = (I - 2P) \\otimes I_N,
    where:
        - I is the identity operator over the ancilla,
        - P is the projector onto the 0 state for the ancilla,
        - I_N is the identity operator over the state register

    """
    return tq.gates.X(target=ancilla) + tq.gates.X(control=ancilla, target=state_qubits) \
           + tq.gates.X(target=ancilla)


def _num_iter(unitaries: list[tuple[float, tq.QCircuit]]) -> int:
    """Return the number of times to apply the amplitude amplificiation to maximize
    success probability"""
    s = sum(pair[0] for pair in unitaries)
    denom = 4 * arcsin(1 / s)
    return floor(pi / denom)

# Testing amplitude amplification

# TODO


def test_algorithm(op: tq.QCircuit,
                   unitaries: list[tuple[float, tq.QCircuit]], ancilla) -> bool:
    """Test whether this works at all. If it doesn't, welp, I'm in trouble."""

    ...

# Main tests


if __name__ == '__main__':
    ...
