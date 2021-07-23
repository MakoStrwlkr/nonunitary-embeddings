"""
CSC299 ROP
LCU implementation with two qubit ancilla
"""

import tequila as tq
# from tequila import TequilaWarning
from typing import Iterable
# import numpy as np
from numpy import arcsin, floor, asarray, pi
from math import sqrt
import random
# import time
import copy

TARGET_FIDELITY = 0.995


class LCUFidelityWarning(Warning):
    """Warning raised when resulting state has low fidelity with target state."""

    def __str__(self) -> str:
        """Return a string representation of this warning."""
        return 'Resulting state may have a lower fidelity with target state than ideally expected'


def lcu(ancilla: list, unitaries: list[tuple[float, tq.QCircuit]], prepare: tq.QCircuit = None) \
        -> tq.QCircuit:
    """Return the circuit corresponding to the prepare operator.

    Preconditions:
        - TODO
    """
    if prepare is None:
        prep, infid = _prepare_operator_optimize_2anc(ancilla, unitaries)
        if 1 - infid.energy < TARGET_FIDELITY:
            raise LCUFidelityWarning
    else:
        prep = prepare

        m = len(ancilla)

        # Define required state
        coefficients = [unit[0] for unit in unitaries]
        normalize = sqrt(sum(coefficients))

        coefficients = [coeff / normalize for coeff in coefficients]

        if len(coefficients) < 2 ** m:
            coefficients.append(0)

        wfn_target = tq.QubitWaveFunction.from_array(asarray(coefficients))
        wfn_target = wfn_target.normalize()

        if tq.simulate(fidelity(wfn_target, prepare)) < TARGET_FIDELITY:
            raise LCUFidelityWarning

    return prep + select_operator(ancilla, unitaries) + prep.dagger()


def _prepare_operator_optimize_2anc(ancilla: list, unitaries: list[tuple[float, tq.QCircuit]]) \
        -> tuple[tq.QCircuit, tq.optimizers]:
    """Return the circuit corresponding to the prepare operator.

    Preconditions:
        - ...
        - TODO
    """
    m = len(ancilla)

    # Define required state
    coefficients = [unit[0] for unit in unitaries]
    normalize = sqrt(sum(coefficients))

    coefficients = [sqrt(coeff) / normalize for coeff in coefficients]

    if len(coefficients) < 2 ** m:
        coefficients.append(0)

    wfn_target = tq.QubitWaveFunction.from_array(asarray(coefficients))
    wfn_target = wfn_target.normalize()

    # Create general parametric circuit
    th, phi, lam, pqc = param_circ(ancilla)
    n_th, n_phi, n_lam = len(th), len(phi), len(lam)

    # Initialize random wfn
    th0 = {key: random.uniform(0, pi) for key in th}
    phi0 = {key: random.uniform(0, pi) for key in phi}
    lam0 = {key: random.uniform(0, pi) for key in lam}
    initial_values = {**th0, **phi0, **lam0}

    # Define (in)fidelity
    # wfn_pqc = tq.simulate(pqc, variables=initial_values)
    inf = 1.0 - fidelity(wfn_target, pqc)

    # Define bounds (if supported)
    min_angles, max_angles = 0, 4 * pi
    bnds_list = [[min_angles, max_angles]]
    for _ in range(len(initial_values)):
        bnds_list.append([min_angles, max_angles])
    th_dict = dict(zip([str(th[i]) for i in range(0, n_th)], bnds_list))
    phi_dict = dict(zip([str(phi[i]) for i in range(0, n_phi)], bnds_list))
    lam_dict = dict(zip([str(lam[i]) for i in range(0, n_lam)], bnds_list))

    bnds = {**th_dict, **phi_dict, **lam_dict}

    # Minimize objective
    # t0 = time.time()
    infid = tq.minimize(objective=inf, initial_values=initial_values, method='TNC',
                        method_bounds=bnds, silent=True)
    # t1 = time.time()

    return pqc, infid


def fidelity(wfn_target: tq.wavefunction.qubit_wavefunction.QubitWaveFunction,
             qc: tq.QCircuit) -> tq.Objective:
    """Return the fidelity between wfn_target and the result of applying """
    rho_targ = tq.paulis.Projector(wfn=wfn_target)
    objective = tq.Objective.ExpectationValue(U=qc, H=rho_targ)
    return objective


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
            circuit += tq.gates.X(target=ancilla[m - i - 1])
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


# Amp amp

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
    if isinstance(state_qubits, list) and isinstance(ancilla, list):
        qubits = list(set(state_qubits + ancilla))
    elif isinstance(state_qubits, list) and not isinstance(ancilla, list):
        qubits = list(set(state_qubits + [ancilla]))
    elif not isinstance(state_qubits, list) and isinstance(ancilla, list):
        qubits = list(set([state_qubits] + ancilla))
    else:
        qubits = list(set([state_qubits] + [ancilla]))

    z_gate, cz_gate = tq.gates.Z(target=qubits), tq.gates.Z(control=ancilla, target=state_qubits)

    return z_gate + cz_gate


def _num_iter(unitaries: list[tuple[float, tq.QCircuit]]) -> int:
    """Return the number of times to apply the amplitude amplificiation to maximize
    success probability"""
    s = sum(pair[0] for pair in unitaries)
    alpha = arcsin(1 / s)
    frac = (pi / 2) / alpha
    return floor(0.5 * (frac - 1))

# if __name__ == '__main__':
#     print(example_func_select())
