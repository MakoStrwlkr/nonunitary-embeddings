"""
CSC299 ROP: General implementation of LCU circuit
"""

# Imports (must add to requirements)
# Requirements: Python 3.9+

import tequila as tq
# import numpy as np
from numpy import pi, sqrt, asarray
import time
from typing import Optional, Any, Iterable, Union
import copy


# Implement complete algorithm (with amp_amp)


# Implement LCU circuit

def lcu(ancilla: Union[str, int, list], unitaries) -> tq.QCircuit:
    """Return the circuit ...

    TODO: Complete docstring

    Preconditions:
        - TODO
    """
    # Must also check if length of ancilla is appropriate. If not, add qubits to ancilla

    # Check if ancilla has only 1 qubit
    if (not isinstance(ancilla, list)) or len(ancilla) == 1:
        # Special case of 1 qubit ancilla
        ...

    # Check if ancilla has only 2 qubits
    elif len(ancilla) == 2:
        # Special case of 2 qubit ancilla
        ...

    # General case
    else:
        ...


# Implement Prepare operator

def prepare_opertor(ancilla: list, unitaries: list[tuple[float, tq.QCircuit]],
                    steps: Optional[int] = 4, debug: Optional[bool] = False) -> tq.QCircuit:
    """Return the circuit corresponding to the select operator

    Preconditions:
        - 2 ** (len(ancilla) - 1) < len(sum_of_unitaries) <= 2 ** len(ancilla)
        - int > 0
    """
    # Define m for convenience
    m = len(ancilla)

    # Define required state
    coefficients = [unit[0] for unit in unitaries]
    normalize = sqrt(sum(coefficients))

    coefficients = [coeff / normalize for coeff in coefficients]

    if len(coefficients) < 2 ** m:
        coefficients.append(0)

    wfn_target = tq.QubitWaveFunction.from_array(asarray(coefficients)).normalize()

    # Define zero state
    zero_state_coeff = [1.0] + [0 for _ in range(len(coefficients) - 1)]
    zero_state = tq.QubitWaveFunction.from_array(asarray(zero_state_coeff))

    # Define generators
    generator_1 = tq.paulis.KetBra(bra=wfn_target, ket=zero_state.normalize())
    generator_2 = tq.paulis.KetBra(ket=wfn_target, bra=zero_state.normalize())

    g = 1.0j * (generator_1 - generator_2)

    circ = tq.gates.Trotterized(generator=g, angle=pi, steps=steps)
    wfn = tq.simulate(circ)

    # print("steps = {}, F = {}".format(steps, abs(wfn.inner(wfn_target)) ** 2))

    circ = tq.QCircuit()
    projector = tq.paulis.Projector(wfn=wfn_target)
    for step in range(steps):
        for ps in g.paulistrings:
            t = tq.Variable((str(ps), step))
            circ += tq.gates.ExpPauli(paulistring=ps, angle=pi / steps * t)
    expect = tq.ExpectationValue(H=projector, U=circ)

    if debug:
        t0 = time.time()
        result = tq.minimize(1 - expect, initial_values=1.0, silent=True)
        t1 = time.time()

        print("steps = {}, F = {}".format(steps, 1.0 - result.energy))
        print("time taken = {}".format(t1 - t0))
    else:

        result = tq.minimize(1 - expect, initial_values=1.0, silent=True)
        print("steps = {}, F = {}".format(steps, 1.0 - result.energy))

    return circ


# Implement Select operator

def select_operator(ancilla: list, sum_of_unitaries: list[tuple[float, tq.QCircuit]]) \
        -> tq.QCircuit:
    """Return the circuit corresponding to the select operator

    Preconditions:
        - 2 ** (len(ancilla) - 1) < len(sum_of_unitaries) <= 2 ** len(ancilla)
    """
    unitaries = [pair[1] for pair in sum_of_unitaries]
    circuit = tq.QCircuit()

    for i in range(len(unitaries)):
        circuit += _controlled_unitary(ancilla, unitaries[i], i)

    return circuit


# Implement amplitude amplification


# Implement special case for exactly 1 qubit in ancilla


# Implement special case of Prepare for exactly 2 qubits in ancilla


# Implement required functions

def _controlled_unitary(ancilla: list, unitary: tq.QCircuit, n: int) -> tq.QCircuit:
    """Return controlled version of unitary with all the qubits in ancilla as the control qubits.

    This function is not an in-place function and it DOES NOT MUTATE unitary;
    Instead, it returns a new tequila circuit.

    Preconditions:
        - ancilla and unitary cannot have any common qubits
        - 0 <= n < 2 ** len(ancilla)
    """
    m = len(ancilla)

    binary_list = _num_to_binary_list(m, n)

    circuit = tq.QCircuit()
    for i in range(len(binary_list)):
        if binary_list[i] == 0:
            circuit += tq.gates.X(target=ancilla[m - i - 1])
    reverse_gates = circuit.dagger()

    circuit += _control_unitary(ancilla, unitary) + reverse_gates

    return circuit


def _control_unitary(ancilla: Union[str, int, list], unitary: tq.QCircuit) -> tq.QCircuit:
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
    binary = tq.BitString.from_int(integer=n, nbits=m)
    binary_list = binary.array

    return binary_list


def fidelity(wfn_target: tq.wavefunction.qubit_wavefunction.QubitWaveFunction,
             qc: tq.QCircuit) -> tq.Objective:
    """Return the fidelity between wfn_target and the result of applying """
    rho_targ = tq.paulis.Projector(wfn=wfn_target)
    objective = tq.Objective.ExpectationValue(U=qc, H=rho_targ)
    return objective


def param_circ(anc: list) \
        -> tuple[list[tq.Variable], list[tq.Variable], list[tq.Variable], tq.QCircuit]:
    """Return a parameterized quantum circuit which acts over two qubits
    with variables theta, phi, lambda for each qubit.

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


# Implement testing functions

# Test 1 qubit ancilla

# Test 2 qubit ancilla

# Test general case

# Test amplitude amplification


# Main testing

if __name__ == '__main__':
    # circuit = tq.gates.X(0) + tq.gates.Y(1) + tq.gates.Z(2) + tq.gates.H(3) + tq.gates.X(4)
    prep_0 = prepare_opertor(ancilla=[5, 6, 7],
                             unitaries=[(1, tq.gates.X(0)),
                                        (1, tq.gates.Y(1)),
                                        (1, tq.gates.Z(2)),
                                        (1, tq.gates.H(3))],
                             steps=1, debug=True)
    prep_1 = prepare_opertor(ancilla=[5, 6, 7],
                             unitaries=[(1, tq.gates.X(0)),
                                        (1, tq.gates.Y(1)),
                                        (1, tq.gates.Z(2)),
                                        (1, tq.gates.H(3))],
                             steps=2, debug=True)
    prep_2 = prepare_opertor(ancilla=[5, 6, 7],
                             unitaries=[(1, tq.gates.X(0)),
                                        (1, tq.gates.Y(1)),
                                        (1, tq.gates.Z(2)),
                                        (1, tq.gates.H(3))],
                             steps=4, debug=True)
    prep_3 = prepare_opertor(ancilla=[5, 6, 7],
                             unitaries=[(1, tq.gates.X(0)),
                                        (1, tq.gates.Y(1)),
                                        (1, tq.gates.Z(2)),
                                        (1, tq.gates.H(3))],
                             steps=8, debug=True)
