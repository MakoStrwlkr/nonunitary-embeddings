"""
CSC299 ROP: General implementation of LCU circuit
"""

# Imports (must add to requirements)
# Requirements: Python 3.9+

import tequila as tq
from warnings import warn
from numpy import pi, sqrt, arcsin, asarray, floor
import time
from typing import Optional, Iterable, Union
import copy
from random import uniform

# Global variables
TARGET_FIDELITY = 0.99


# LCU Class

class LCU:
    """
    A custom data type for the LCU circuit

    Public (immutable) properties:
        - ancilla: list of qubits in the ancilla register
        - prep_circ: circuit corresponding to prepare operator for LCU circuit
        - select_circ: circuit corresponding to select operator in LCU circuit
        - lcu_circ: LCU circuit, excluding amplitude amplification procedure
        - amp_amp: amplitude amplification operator, signifying one step of amplitude amplification
        - full_circ: full LCU circuit, including amplitude amplification procedure

    """
    # Private Instance Attributes:
    #   - ancilla: list of qubits in the ancilla register
    #   - prep_circ: circuit corresponding to prepare operator for LCU circuit
    #   - select_circ: circuit corresponding to select operator in LCU circuit
    #   - lcu_circ: LCU circuit, excluding amplitude amplification procedure
    #   - amp_amp: amplitude amplification operator, signifying one step of amplitude amplification
    #   - full_circ: full LCU circuit, including amplitude amplification procedure

    _ancilla: list[Union[int, str]]
    _prep_circ: tq.QCircuit
    _select_circ: tq.QCircuit
    _lcu_circ: tq.QCircuit
    _amp_amp: tq.QCircuit
    _full_circ: tq.QCircuit

    def __init__(self, ancilla: Union[list[Union[str, int]], str, int],
                 unitaries: list[tuple[float, tq.QCircuit]],
                 prepare: Optional[tq.QCircuit] = None) -> None:
        """Initialize the LCU object such that it stores the various parts of the LCU algorithm.

        Preconditions:
            - TODO
        """
        self._ancilla = list(ancilla)

        if prepare:
            self._prep_circ = prepare
        else:
            self._prep_circ = prepare_operator(ancilla, unitaries)

        self._select_circ = select_operator(ancilla, unitaries)

        self._lcu_circ = lcu(ancilla, unitaries, prepare)

        self._amp_amp = amp_amp_op(self._lcu_circ, ancilla)

        self._full_circ = ...

    @property
    def ancilla(self) -> list:
        """Return a list of the qubits in the ancilla register"""
        return self._ancilla

    @property
    def prepare(self) -> tq.QCircuit:
        """Return the prepare operator for the LCU algorithm"""
        return self._prep_circ

    @property
    def select(self) -> tq.QCircuit:
        """Return the select operator for the LCU algorithm"""
        return self._select_circ

    @property
    def lcu(self) -> tq.QCircuit:
        """Return the circuit for the LCU algorithm,
        excluding the amplitude amplification procedure"""
        return self._lcu_circ

    @property
    def amp_amp(self) -> tq.QCircuit:
        """Return the amplitude amplification operator used in the LCU algorithm,
        which symbolizes one step in the generalized amplitude amplification procedure"""
        return self._amp_amp

    @property
    def full_circuit(self) -> tq.QCircuit:
        """Return the full circuit for the LCU algorithm,
        including the amplitude amplification procedure"""
        return self._full_circ


# Warning(s) class(es) - Add more as required

class LCUFidelityWarning(UserWarning):
    """Warning raised when resulting state has low fidelity with target state."""

    def __str__(self) -> str:
        """Return a string representation of this warning."""
        return 'Resulting state may have a lower fidelity with target state than ideally expected'


class LCUMismatchWarning(UserWarning):
    """Warning raised when length of ancilla does not match with the number of unitaries."""

    def __str__(self) -> str:
        """Return a string representation of this warning."""
        return 'Length of given ancilla does not match with the number of unitaries'


class UnnecessaryLCUWarning(UserWarning):
    """Warning raised when implementing the LCU algorithm is not necessary."""

    def __str__(self) -> str:
        """Return a string representation of this warning."""
        return 'Implementing the LCU algorithm is not necessary in this case'


# Implement complete algorithm (with amp_amp)

def nonunitary_embedding(ancilla: Union[str, int, list[Union[str, int]]],
                         unitaries: list[tuple[float, tq.QCircuit]],
                         prepare: Optional[tq.QCircuit] = None,
                         tolerance: Optional[float] = TARGET_FIDELITY,
                         easter_egg: Optional[bool] = False) -> tq.QCircuit:
    """Return the complete circuit of the algorithm for nonunitary embeddings, which includes
    the lcu circuit followed by the amplitude amplification procedure.

    Preconditions:
        - TODO

    """

    lcu_circ = lcu(ancilla=ancilla, unitaries=unitaries,
                   prepare=prepare, tolerance=tolerance, easter_egg=easter_egg)

    amp_amp_circ = amp_amp(unitaries=unitaries, walk_op=lcu_circ, ancilla=ancilla)

    return lcu_circ + amp_amp_circ


# Implement LCU circuit

def lcu(ancilla: Union[str, int, list[Union[str, int]]],
        unitaries: list[tuple[float, tq.QCircuit]],
        prepare: Optional[tq.QCircuit] = None,
        tolerance: Optional[float] = TARGET_FIDELITY,
        easter_egg: Optional[bool] = True) -> tq.QCircuit:
    """Return the circuit for the linear combinations os unitaries algorithm, excluding the
    procedure for amplitude amplification.

    Preconditions:
        - prepare is not a parametrized circuit
        - ...
    """
    anc = ancilla if isinstance(ancilla, list) else [ancilla]

    if len(unitaries) > 2 ** len(anc):
        warn('Size of ancilla is too small. Adding extra qubits...', LCUMismatchWarning)

        # TODO: Add extra qubits

    if prepare is not None:
        wfn_target = _target_wfn(ancilla=ancilla, unitaries=unitaries)

        # If satisfactory fidelity not achieved, raise a warning
        fid = tq.simulate(fidelity(wfn_target, prepare))
        if fid < tolerance:
            warn(f'Target fidelity of {tolerance} not achieved. Fidelity was {fid}',
                 LCUFidelityWarning)

        # if someone mistakenly tries to do this
        if len(unitaries) == 1:

            # Small Easter Egg... If you spotted this, congratulations! Have a bit of humor!
            if easter_egg:
                warn("Uhh... You really don't need to use the LCU algorithm for this...",
                     UnnecessaryLCUWarning)

            return unitaries[0][1]

        # Check if ancilla has only 1 qubit
        elif (not isinstance(ancilla, list)) or len(ancilla) == 1:
            # Special case of 1 qubit ancilla
            unitary_0, unitary_1 = unitaries[0][1], unitaries[1][1]
            return prepare + _select_1ancilla(ancilla, unitary_0, unitary_1) + prepare.dagger()

        return prepare + select_operator(ancilla, unitaries) + prepare.dagger()
    else:
        return _lcu_no_prepare(ancilla, unitaries)


def _lcu_no_prepare(ancilla: Union[str, int, list[Union[str, int]]],
                    unitaries: list[tuple[float, tq.QCircuit]],
                    tolerance: Optional[float] = TARGET_FIDELITY,
                    easter_egg: Optional[bool] = False) -> tq.QCircuit:
    """Return the circuit ...

    TODO: Complete docstring

    Preconditions:
        - ...
    """
    # Must also check if length of ancilla is appropriate. If not, add qubits to ancilla.
    # Raise warning if not.

    # if someone mistakenly tries to do this (stupidly)
    if len(unitaries) == 1:

        # Small Easter Egg... If you spotted this, congratulations! Have a bit of humor!
        if easter_egg:
            print("Uhh... You really don't need to use the LCU algorithm for this...")
            return unitaries[0][1]

        raise UnnecessaryLCUWarning

    # Check if ancilla has only 1 qubit
    elif (not isinstance(ancilla, list)) or len(ancilla) == 1:
        # Special case of 1 qubit ancilla
        return lcu_1ancilla(ancilla, unitaries)

    # Check if ancilla has only 2 qubits
    elif len(ancilla) == 2:
        # Special case for 2 qubit ancilla
        prep = _prepare_2ancilla(ancilla=ancilla, unitaries=unitaries, tolerance=tolerance)

    # General case
    else:
        prep = prepare_operator(ancilla, unitaries)

    return prep + select_operator(ancilla, unitaries) + prep.dagger()


# Implement Prepare operator

def prepare_operator(ancilla: list, unitaries: list[tuple[float, tq.QCircuit]],
                     steps: Optional[int] = 1, tolerance: Optional[float] = TARGET_FIDELITY,
                     debug: Optional[bool] = False) -> tq.QCircuit:
    """Return the circuit corresponding to the prepare operator

    Preconditions:
        - 2 ** (len(ancilla) - 1) < len(sum_of_unitaries) <= 2 ** len(ancilla)
        - int > 0
    """
    # wfn_target = _target_wfn(ancilla=ancilla, unitaries=unitaries)

    m = len(ancilla)

    # Define required state
    coefficients = [unit[0] for unit in unitaries]
    normalize = sqrt(sum(coefficients))

    coefficients = [coeff / normalize for coeff in coefficients]

    if len(coefficients) < 2 ** m:
        extension = [0 for _ in range(2 ** m - len(coefficients) + 1)]
        coefficients.extend(extension)

    wfn_target = tq.QubitWaveFunction.from_array(asarray(coefficients)).normalize()

    # Define zero state
    m = len(ancilla)

    # Define required state
    coefficients = [unit[0] for unit in unitaries]
    normalize = sqrt(sum(coefficients))

    coefficients = [coeff / normalize for coeff in coefficients]

    if len(coefficients) < 2 ** m:
        extension = [0 for _ in range(2 ** m - len(coefficients) + 1)]
        coefficients.extend(extension)

    # Define zero state

    zero_state_coeff = [1.0] + [0 for _ in range(len(coefficients) - 1)]
    zero_state = tq.QubitWaveFunction.from_array(asarray(zero_state_coeff))

    # Define generators
    generator_1 = tq.paulis.KetBra(bra=wfn_target, ket=zero_state.normalize())
    generator_2 = tq.paulis.KetBra(ket=wfn_target, bra=zero_state.normalize())

    g = 1.0j * (generator_1 - generator_2)

    # Don't remove this line! It's required!!
    # assert g.is_hermitian()

    tq.simulate(tq.gates.Trotterized(generator=g, angle=pi / 2, steps=1))

    circ = tq.QCircuit()
    projector = tq.paulis.Projector(wfn=wfn_target)

    for step in range(steps):
        for ps in g.paulistrings:
            t = tq.Variable((str(ps), step))
            circ += tq.gates.ExpPauli(paulistring=ps, angle=pi / steps * t)

    # Define objective function to be fidelity, based on the expectation value
    expect = tq.ExpectationValue(H=projector, U=circ)

    if debug:
        t0 = time.time()

        # Minimize the infidelity
        result = tq.minimize(1 - expect, initial_values=1.0, silent=True)
        t1 = time.time()

        print("steps = {}, F = {}".format(steps, 1.0 - result.energy))
        print("time taken = {}".format(t1 - t0))
    else:

        result = tq.minimize(1 - expect, initial_values=1.0, silent=True)

    # If satisfactory fidelity not achieved, raise a warning
    fid = 1 - result.energy
    if fid < tolerance:
        warn(f'Target fidelity of {tolerance} not achieved. Fidelity was {fid}',
             LCUFidelityWarning)

    return circ.map_variables(variables=result.variables)


# Implement Select operator

def select_operator(ancilla: list, unitaries: list[tuple[float, tq.QCircuit]]) \
        -> tq.QCircuit:
    """Return the circuit corresponding to the select operator

    Preconditions:
        - 2 ** (len(ancilla) - 1) < len(sum_of_unitaries) <= 2 ** len(ancilla)
    """
    unitaries = [pair[1] for pair in unitaries]
    circuit = tq.QCircuit()

    for i in range(len(unitaries)):
        circuit += _controlled_unitary(ancilla, unitaries[i], i)

    return circuit


# Implement amplitude amplification

# TODO

def amp_amp(unitaries: list[tuple[float, tq.QCircuit]], walk_op: tq.QCircuit, ancilla) \
        -> tq.QCircuit:
    """Return the amplitude amplification procedure obtained by repeating
    the amplitude amplification step for a total of s times where s is the
    result of function _num_iter()
    """
    amplification_operator = amp_amp_op(walk_op, ancilla)
    s = _num_iter(unitaries)

    sum_of_steps = tq.QCircuit()
    for _ in range(s):
        sum_of_steps += amplification_operator

    return sum_of_steps


def amp_amp_op(walk_op: tq.QCircuit, ancilla: Union[list[Union[str, int]], str, int]) \
        -> tq.QCircuit:
    """Return WRW.dagger()R,
     where R is the reflect operator returned by the func reflect_operator"""
    anc_qubits = ancilla if isinstance(ancilla, list) else [ancilla]
    state_qubits = [qubit for qubit in walk_op.qubits if qubit not in anc_qubits]

    reflect = reflect_operator(state_qubits=state_qubits, ancilla=ancilla)

    return reflect + walk_op.dagger() + reflect + walk_op


def reflect_operator(state_qubits, ancilla: Union[list[Union[str, int]], str, int]) \
        -> tq.QCircuit:
    """
    Return the reflection operator R = (I - 2P) \\otimes I_N,
    where:
        - I is the identity operator over the ancilla,
        - P is the projector onto the 0 state for the ancilla,
        - I_N is the identity operator over the state register

    """
    x_gate, cnot_gate = tq.gates.X(target=ancilla), tq.gates.X(control=ancilla, target=state_qubits)

    return x_gate + cnot_gate + x_gate


# Implement special case for exactly 1 qubit in ancilla

def lcu_1ancilla(ancilla: Union[str, int],
                 unitaries: list[tuple[float, tq.QCircuit]],
                 tolerance: Optional[float] = TARGET_FIDELITY) -> tq.QCircuit:
    """The main part.

    TODO: write proper docstring
    """
    prepare = _prepare_1ancilla(ancilla, unitaries)

    wfn_target = _target_wfn(ancilla=list(ancilla), unitaries=unitaries)

    fid = tq.simulate(fidelity(wfn_target, prepare))

    if fid < tolerance:
        warn(f'Target fidelity of {tolerance} not achieved. Fidelity was {fid}.',
             LCUFidelityWarning)

    circ = prepare + _select_1ancilla(ancilla, unitaries[0][1], unitaries[1][1]) + prepare.dagger()
    return circ


def _prepare_1ancilla(ancilla: Union[list[Union[str, int]], str, int],
                      unitaries: list[tuple[float, tq.QCircuit]]) -> tq.QCircuit:
    """
    Prepare operator, when the Hamiltonian can be expressed as the linear combination of two
    unitary operators.
    Requires only one ancillary qubit.

    Preconditions:
        - ...

    TODO: Complete docstring
    """
    alpha_0, alpha_1 = unitaries[0][0], unitaries[1][0]

    theta = -2 * arcsin(sqrt(alpha_1 / (alpha_0 + alpha_1)))

    return tq.gates.Ry(target=ancilla, angle=theta)


def _select_1ancilla(ancillary, unitary_0: tq.QCircuit, unitary_1: tq.QCircuit) -> tq.QCircuit:
    """
    Select operator, when the Hamiltonian can be expressed as the linear combination of two
    unitary operators.
    Requires only one ancillary qubit.

    Returns ...

    TODO: Complete docstring
    """
    impl_1 = _control_unitary(ancilla=ancillary, unitary=unitary_1)

    x_gate = tq.gates.X(target=ancillary)
    control = _control_unitary(ancilla=ancillary, unitary=unitary_0)

    impl_0 = x_gate + control + x_gate

    return impl_1 + impl_0


# Implement special case of Prepare for exactly 2 qubits in ancilla


def lcu_2ancilla(ancilla: list[Union[str, int]],
                 unitaries: list[tuple[float, tq.QCircuit]],
                 prepare: Optional[tq.QCircuit] = None,
                 tolerance: Optional[float] = TARGET_FIDELITY) -> tq.QCircuit:
    """Return the circuit corresponding to the prepare operator.

    Preconditions:
        - 0 < tolerance < 1
        - prepare is not a parametrized circuit
        - TODO
    """
    if prepare is None:
        prep, infid = _prepare_2ancilla(ancilla, unitaries)
        if 1 - infid.energy < tolerance:
            raise LCUFidelityWarning

    else:
        prep = prepare

        wfn_target = _target_wfn(ancilla=ancilla, unitaries=unitaries)

        # If satisfactory fidelity not achieved, raise a warning
        fid = tq.simulate(fidelity(wfn_target, prepare))

        if fid < tolerance:
            warn('...', LCUFidelityWarning)

    return prep + select_operator(ancilla, unitaries) + prep.dagger()


def _prepare_2ancilla(ancilla: list[Union[str, int]],
                      unitaries: list[tuple[float, tq.QCircuit]],
                      debug: Optional[bool] = False,
                      tolerance: Optional[float] = TARGET_FIDELITY) \
        -> tuple[tq.QCircuit, tq.optimizers]:
    """Return the circuit corresponding to the prepare operator.

    Preconditions:
        - len(ancilla) == 2
    """
    # Define required state
    wfn_target = _target_wfn(ancilla=ancilla, unitaries=unitaries)

    # Create general parametric circuit
    th, phi, lam, pqc = param_circ(ancilla)
    n_th, n_phi, n_lam = len(th), len(phi), len(lam)

    # Initialize parameters randomly
    th0 = {key: uniform(0, pi) for key in th}
    phi0 = {key: uniform(0, pi) for key in phi}
    lam0 = {key: uniform(0, pi) for key in lam}
    initial_values = {**th0, **phi0, **lam0}

    # Define (in)fidelity
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

    if debug:
        t0 = time.time()
        result = tq.minimize(objective=inf, initial_values=initial_values, method='TNC',
                             method_bounds=bnds, silent=True)
        t1 = time.time()

        print('time taken: {}'.format(t1 - t0))
    else:
        result = tq.minimize(objective=inf, initial_values=initial_values, method='TNC',
                             method_bounds=bnds, silent=True)

    # If satisfactory fidelity not achieved, raise a warning
    fid = 1 - result.energy

    if fid < tolerance:
        warn(f'Target fidelity of {tolerance} not achieved. Fidelity was {fid}',
             LCUFidelityWarning)

    return pqc, result


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
        - len(anc) = 2
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
    ugate = tq.QCircuit()
    ugate += tq.gates.Rz(target=q0, angle=phi) + tq.gates.Ry(target=q0, angle=th)
    ugate += tq.gates.Rz(target=q0, angle=lam)
    return ugate


def _target_wfn(ancilla: list[Union[str, int]],
                unitaries: list[tuple[float, tq.QCircuit]]) -> tq.QubitWaveFunction:
    """Return the ideal target wavefunction expected after applying the prepare operator to the
    zero state of the ancilla register"""
    m = len(ancilla)

    # Define required state
    coefficients = [unit[0] for unit in unitaries]
    normalize = sqrt(sum(coefficients))

    coefficients = [coeff / normalize for coeff in coefficients]

    if len(coefficients) < 2 ** m:
        extension = [0 for _ in range(2 ** m - len(coefficients) + 1)]
        coefficients.extend(extension)

    wfn_target = tq.QubitWaveFunction.from_array(asarray(coefficients)).normalize()

    return wfn_target


def _num_iter(unitaries: list[tuple[float, tq.QCircuit]]) -> int:
    """Return the number of times to apply the amplitude amplificiation to maximize
    success probability"""
    s = sum(pair[0] for pair in unitaries)
    denom = 4 * arcsin(1 / s)
    return floor(pi / denom)


def _add_qubits(anc: list[Union[str, int]], length: int) -> list[Union[str, int]]:
    """Add more items to """


# Implement testing functions

# TODO

# Test 1 qubit ancilla

# Test 2 qubit ancilla

# Test general case

# Test amplitude amplification


# Main testing

if __name__ == '__main__':
    # circuit = tq.gates.X(0) + tq.gates.Y(1) + tq.gates.Z(2) + tq.gates.H(3) + tq.gates.X(4)
    prep_0 = prepare_operator(ancilla=[5, 6, 7],
                              unitaries=[(0.1, tq.gates.X(0)),
                                         (0.2, tq.gates.Y(1)),
                                         (0.3, tq.gates.Z(2)),
                                         (0.4, tq.gates.H(3))],
                              steps=1, debug=True)
    prep_1 = prepare_operator(ancilla=[5, 6, 7],
                              unitaries=[(1, tq.gates.X(0)),
                                         (1, tq.gates.Y(1)),
                                         (1, tq.gates.Z(2)),
                                         (1, tq.gates.H(3))],
                              steps=2, debug=True)
    prep_2 = prepare_operator(ancilla=[5, 6, 7],
                              unitaries=[(1, tq.gates.X(0)),
                                         (1, tq.gates.Y(1)),
                                         (1, tq.gates.Z(2)),
                                         (1, tq.gates.H(3))],
                              steps=4, debug=True)
    prep_3 = prepare_operator(ancilla=[5, 6, 7],
                              unitaries=[(1, tq.gates.X(0)),
                                         (1, tq.gates.Y(1)),
                                         (1, tq.gates.Z(2)),
                                         (1, tq.gates.H(3))],
                              steps=8, debug=True)
