"""
CSC299 ROP
Testing file
"""

# Import statements
import tequila as tq
from math import sqrt, isclose
from numpy import arcsin, asarray
from lcu_v1 import LCU, lcu_1ancilla, prepare_operator, reflect_operator
from typing import Union
from random import random


# Implement testing functions

# TODO

############################################################################################
# Test 1 qubit ancilla
############################################################################################

#####################
# Test prepare
#####################

# For reference, prepare operator is replicated here
def _prepare_1ancilla(ancilla: Union[list[Union[str, int]], str, int],
                      unitaries: list[tuple[float, tq.QCircuit]]) -> tq.QCircuit:
    """
    Prepare operator, when the Hamiltonian can be expressed as the linear combination of two
    unitary operators.
    Requires only one ancillary qubit.

    Preconditions:
        - if isinstance(ancilla, list): len(ancilla) == 1
        - all(coeff > 0 for coeff in [pair[0] for pair in unitaries])
    """
    alpha_0, alpha_1 = unitaries[0][0], unitaries[1][0]

    theta = 2 * arcsin(sqrt(alpha_1 / (alpha_0 + alpha_1)))

    return tq.gates.Ry(target=ancilla, angle=theta)


def test_prepare_1anc() -> None:
    """Test whether the prepare operator for 1 qubit ancilla works as intended

    Preconditions:
        - 0 not in state_prep.qubits
    """
    num0 = random()
    num1 = 1 - num0
    unitaries = [(num0, tq.gates.X(1)), (num1, tq.gates.Z(1))]
    prepare = _prepare_1ancilla(0, unitaries)

    wfn_target = tq.QubitWaveFunction.from_array(asarray([sqrt(num0), sqrt(num1)]))

    projector = tq.paulis.Projector(wfn=wfn_target)

    expval = tq.ExpectationValue(H=projector, U=prepare)

    p = tq.simulate(expval)

    # print(p)
    assert isclose(p, 1, abs_tol=0.001)


#####################
# Test select
#####################

# For reference, select operator is replicated below
def _select_1ancilla(ancillary: Union[str, int], unitary_0: tq.QCircuit, unitary_1: tq.QCircuit) \
        -> tq.QCircuit:
    """
    Select operator, when the Hamiltonian can be expressed as the linear combination of two
    unitary operators.
    Requires only one ancillary qubit.

    Returns the select operator.

    Preconditions:
        - ancillary not in unitary_0.qubits
        - ancillary not in unitary_1.qubits
    """
    impl_1 = unitary_1.add_controls(ancillary, inpl=False)

    x_gate = tq.gates.X(target=ancillary)
    control = unitary_0.add_controls(ancillary, inpl=False)

    impl_0 = x_gate + control + x_gate

    return impl_1 + impl_0


def test_select_1anc() -> None:
    """Test whether the select operator for the 1 qubit ancilla works as expected."""
    ...

# Test LCU


############################################################################################
# Test 2 qubit ancilla
############################################################################################

# Test prepare


############################################################################################
# Test general case
############################################################################################

# Test prepare

def test_prepare_general() -> None:
    """Test whether the general case for the prepare operator works as expected."""
    n = 0
    while n == 0:
        n = int(4 * random())

    length = 2 ** n
    ancilla = list(range(length))
    coeffs = [0.5 + 0.49 * random() for _ in range(length)]

    unitaries = [(coeff, tq.gates.X(16)) for coeff in coeffs]

    norm = sqrt(sum(coeffs))
    coefficients = [sqrt(coeff) / norm for coeff in coeffs]
    wfn_target = tq.QubitWaveFunction.from_array(asarray(coefficients))

    prepare = prepare_operator(ancilla=ancilla, unitaries=unitaries)
    print(tq.simulate(prepare))
    expval = tq.ExpectationValue(H=projector, U=amp_amp)

    p = tq.simulate(expval)

    print(p)
    print(wfn_target)
    fid = abs(wfn_target.inner(tq.simulate(prepare))) ** 2

    assert isclose(fid, 1, abs_tol=0.001)


# Test select


############################################################################################
# Amplitude amplification
############################################################################################

def amp_amp_op(walk_op: tq.QCircuit, ancilla) -> tq.QCircuit:
    """Return W R W.dagger() R,
     where R is the reflect operator returned by the function reflect_operator"""
    anc_qubits = ancilla if isinstance(ancilla, list) else [ancilla]
    state_qubits = [qubit for qubit in walk_op.qubits if qubit not in anc_qubits]

    reflect = reflect_operator(state_qubits=state_qubits, ancilla=ancilla)

    return reflect + walk_op.dagger() + reflect + walk_op


#####################
# Test reflection
#####################

def test_reflection() -> None:
    """Test whether the reflection operator works as intended."""
    circ = tq.gates.H(0) + reflect_operator(1, 0)
    wfn_target = tq.QubitWaveFunction.from_array(asarray([-1 / sqrt(2), 0, 1 / sqrt(2), 0]))

    projector = tq.paulis.Projector(wfn=wfn_target)

    expval = tq.ExpectationValue(H=projector, U=circ)

    p = tq.simulate(expval)

    assert isclose(p, 1, rel_tol=0.001)


#####################
# Test amp_amp_op
#####################

def test_amp_amp_stationary() -> None:
    """Test whether the non-trivial stationary angles for amp_amp_op are as expected"""
    # Create uniform superposition of ancilla: qubit 0
    state_prep = tq.gates.H(0)

    # Define target state
    wfn_target = tq.QubitWaveFunction.from_array(asarray([- 1 / sqrt(2), 0, 1 / sqrt(2), 0]))

    # Define amplitude amplification circuit
    amp_amp = amp_amp_op(state_prep, 0)
    circ = amp_amp + state_prep

    fid = abs(wfn_target.inner(tq.simulate(circ))) ** 2

    assert isclose(fid, 1.0, abs_tol=0.001)


def test_amp_amp_change() -> None:
    """Test whether the amp_amp_op works as expected for non-stationary angles"""
    # TODO
    # Create uniform superposition of ancilla: qubit 0
    state_prep = tq.gates.H(0)

    alpha_0 = random()

    while alpha_0 not in {0, 1, 1 / sqrt(2)}:
        alpha_0 = random()

    alpha_1 = sqrt(1 - (alpha_0 ** 2))

    # Define target state
    wfn_target = tq.QubitWaveFunction.from_array(asarray([1 / sqrt(2), 0, - 1 / sqrt(2), 0]))

    # Define amplitude amplification circuit
    amp_amp = amp_amp_op(state_prep, 0)
    projector = tq.paulis.Projector(wfn=wfn_target)

    expval = tq.ExpectationValue(H=projector, U=amp_amp)

    p = tq.simulate(expval)

    print(p)


if __name__ == '__main__':
    # import pytest
    # pytest.main()
    test_prepare_general()
