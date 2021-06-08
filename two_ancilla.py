"""
CSC299 ROP
LCU implementation with two qubit ancilla
"""

import tequila as tq
from typing import Iterable, Optional
import numpy as np
from math import sqrt
import random
import time
import copy


# def prepare_operator_2ancilla(ancilla: list, unitaries: list[tuple[float, tq.QCircuit]]) \
#         -> tq.QCircuit:
#     """Return the circuit corresponding to the prepare operator, for the case of two ancilla."""
#     # TODO
#     # Define required state
#     coefficients = [unit[0] for unit in unitaries]
#     normalize = sqrt(sum(coefficients))
#
#     coefficients = [coeff / normalize for coeff in coefficients]
#
#     if len(coefficients) < 4:
#         coefficients.append(0)
#
#     # use trig identities for cosAcosB, sinAsinB, cosAsinB, sinAcosB.
#
#     ...


def prepare_operator_optimize_2anc(ancilla: list, unitaries: list[tuple[float, tq.QCircuit]]) \
        -> tq.QCircuit:
    """Return the circuit corresponding to the prepare operator.

    Preconditions:
        - ...
        - TODO
    """
    # TODO
    m = len(ancilla)

    # Define required state
    coefficients = [unit[0] for unit in unitaries]
    normalize = sqrt(sum(coefficients))

    coefficients = [coeff / normalize for coeff in coefficients]

    if len(coefficients) < 2 ** m:
        coefficients.append(0)

    wfn_target = tq.QubitWaveFunction.from_array(np.asarray(coefficients))
    wfn_target = wfn_target.normalize()

    # Create general parametric circuit
    th, phi, lam, pqc = param_circ(ancilla)
    n_th, n_phi, n_lam = len(th), len(phi), len(lam)

    # Initialize random wfn
    th0 = {key: random.uniform(0, np.pi) for key in th}
    phi0 = {key: random.uniform(0, np.pi) for key in phi}
    lam0 = {key: random.uniform(0, np.pi) for key in lam}
    all_values = {**th0, **phi0, **lam0}

    # Define (in)fidelity
    wfn_pqc = tq.simulate(pqc, variables=all_values)
    inf = fidelity(wfn_target, pqc, is_inner_pdt=False)

    # Minimize objective
    max_angles = 2 * np.pi
    min_angles = 0
    bnds_list = [[min_angles, max_angles]]
    for _ in range(len(all_values)):
        bnds_list.append([min_angles, max_angles])
    bnds = dict(zip([str(th[i]) for i in range(0, n_th)], bnds_list))
    bnds = {**bnds, **dict(zip([str(phi[i]) for i in range(0, n_phi)], bnds_list))}
    bnds = {**bnds, **dict(zip([str(lam[i]) for i in range(0, n_lam)], bnds_list))}

    # TODO


def fidelity(wfn_target: tq.wavefunction.qubit_wavefunction.QubitWaveFunction,
             qc: tq.QCircuit, is_inner_pdt: Optional[bool] = True) -> float:
    """Return the fidelity ..."""
    if is_inner_pdt:
        wfn_qc = tq.simulate(qc)
        return abs(wfn_target.inner(wfn_qc)) ** 2

    else:
        rho_targ = tq.paulis.Projector(wfn=wfn_target)
        objective = tq.Objective.ExpectationValue(U=qc, H=rho_targ)
        return float(tq.simulate(objective))


def prepare_operator_optimization(ancilla: list, unitaries: list[tuple[float, tq.QCircuit]]) \
        -> tq.QCircuit:
    """Return the circuit corresponding to the prepare operator.

    Preconditions:
        - ...
        - TODO
    """
    # TODO
    m = len(ancilla)

    # Define required state
    coefficients = [unit[0] for unit in unitaries]
    normalize = sqrt(sum(coefficients))

    coefficients = [coeff / normalize for coeff in coefficients]

    if len(coefficients) < 2 ** m:
        coefficients.append(0)

    wfn_target = tq.QubitWaveFunction.from_array(np.asarray(coefficients))
    wfn_target = wfn_target.normalize()

    # Create parametrized circuit

    # Objective function: infidelity

    # Minimize infidelity
    ...


def param_circ(anc: list) \
        -> tuple[list[tq.Variable], list[tq.Variable], list[tq.Variable], tq.QCircuit]:
    """Return a parameterized circuit with ...

    Preconditions:
        - TODO
    """

    # define variables
    th = [tq.Variable(name='theta_{}'.format(i)) for i in range(0, 4)]
    phi = [tq.Variable(name='phi_{}'.format(i)) for i in range(0, 4)]
    lam = [tq.Variable(name='lam_{}'.format(i)) for i in range(0, 4)]

    # PQC
    pqc = unitary_gate(th[0], phi[0], lam[0], anc[0]) + unitary_gate(th[1], phi[1], lam[1], anc[1])
    pqc += tq.gates.CNOT(control=anc[0], target=anc[1])
    pqc += unitary_gate(th[2], phi[2], lam[2], anc[0]) + unitary_gate(th[3], phi[3], lam[3], anc[1])

    return (th, phi, lam, pqc)


def unitary_gate(th, phi, lam, q0) -> tq.QCircuit:
    """Return a particular quantum gate that is not included in the basic gate set"""
    ugate = tq.gates.Rz(target=q0, angle=phi) + tq.gates.Ry(target=q0, angle=th) + tq.gates.Rz(
        target=q0, angle=lam)
    return ugate


# draw the circuit for some (th,phi,lam) for qubit "0"
# tq.draw(unitary_gate(np.pi, np.pi / 4.0, np.pi / 2.0, 0))


def select_operator(ancilla: list, sum_of_unitaries: list[tuple[float, tq.QCircuit]]) \
        -> tq.QCircuit:
    """Return the circuit corresponding to the select operator

    Preconditions:
        - 2 ** (len(ancilla) - 1) < len(sum_of_unitaries) <= 2 ** len(ancilla)
    """
    unitaries = [pair[1] for pair in sum_of_unitaries]
    circuit = tq.QCircuit()

    for i in range(len(unitaries)):
        circuit += controlled_unitary(ancilla, unitaries[i], i)

    return circuit


def controlled_unitary(ancilla: list, unitary: tq.QCircuit, n: int) -> tq.QCircuit:
    """Return controlled version of unitary

    SHOULD NOT mutate unitary

    Preconditions:
        - ancilla and unitary cannot have any common qubits
        - 0 <= n < 2 ** len(ancilla)
    """
    m = len(ancilla)

    binary_list = _num_to_binary_list(m, n)
    # print(binary_list)

    # assert all([digit == 0 or digit == 1 for digit in binary_list])

    circuit = tq.QCircuit()
    for i in range(len(binary_list)):
        if binary_list[i] == 0:
            circuit += tq.gates.X(target=ancilla[m-i-1])
    reverse_gates = circuit.dagger()

    circuit += _control_unitary(ancilla, unitary)
    circuit += reverse_gates

    return circuit


def _control_unitary(ancilla, unitary: tq.QCircuit) -> tq.QCircuit:
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


def _num_to_binary_list(m: int, n: int) -> list[int]:
    """Return the binary representation of n in the form of a list of length m.

    Note: bin(int) exists but returns str

    Preconditions:
        - 2 ** m > n
    """
    # binary = bin(n)[2:]
    # binary_list = [int(digit) for digit in binary]

    binary = tq.BitString.from_int(integer=n, nbits=m)
    binary_list = binary.array

    # if len(binary_list) < m:
    #     k = len(binary_list)
    #     extend = [0 for _ in range(m - k + 1)]
    #     binary_list = extend + binary_list

    return binary_list


def example_func_select() -> tq.QCircuit:
    """Test select"""
    ancilla = [0, 1, 2]
    # identity = tq.QCircuit()
    unitaries = [(0.5, tq.gates.H(3)), (0.5, tq.gates.Z(3)), (0.5, tq.gates.X(4)),
                 (0.5, tq.gates.Y(4)), (0.5, tq.gates.H(3))]
    return select_operator(ancilla=ancilla, sum_of_unitaries=unitaries)


# if __name__ == '__main__':
#     print(example_func_select())
