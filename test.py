"""
CSC299 ROP
Testing file
"""

# Import statements
import tequila as tq
from math import sqrt, isclose
from numpy import arcsin, asarray
from lcu_v1 import amp_amp_op, LCU, lcu_1ancilla, reflect_operator
from typing import Union, Optional
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

    theta = -2 * arcsin(sqrt(alpha_1 / (alpha_0 + alpha_1)))

    return tq.gates.Ry(target=ancilla, angle=theta)


def test_prepare_1anc() -> None:
    """Test whether the prepare operator for 1 qubit ancilla works as intended

    Preconditions:
        - 0 not in state_prep.qubits
    """
    num0 = random()
    num1 = sqrt(1 - (num0 ** 2))
    unitaries = [(num0, tq.gates.X(1)), (num1, tq.gates.Z(1))]
    prepare = _prepare_1ancilla(0, unitaries)

    wfn_target = tq.QubitWaveFunction.from_array(asarray([]))

    projector = tq.paulis.Projector(wfn=wfn_target)

    expval = tq.ExpectationValue(H=projector, U=...)

    assert ...

# Test select

# Test LCU

############################################################################################
# Test 2 qubit ancilla
############################################################################################

# Test prepare

# Test select

############################################################################################
# Test general case
############################################################################################

# Test prepare

# Test select


############################################################################################
# Amplitude amplification
############################################################################################

# Test reflection

def test_reflection() -> None:
    """Test whether the reflection operator works as intended."""
    circ = tq.gates.H(0) + reflect_operator(1, 0)
    wfn_target = tq.QubitWaveFunction.from_array(asarray([-1 / sqrt(2), 0, 1 / sqrt(2), 0]))

    projector = tq.paulis.Projector(wfn=wfn_target)

    expval = tq.ExpectationValue(H=projector, U=circ)

    p = tq.simulate(expval)

    assert isclose(p, 1, rel_tol=0.001)


# Test amp_amp_op

def test_amp_amp_stationary() -> Union[float, tq.QubitWaveFunction]:
    """Test whether the non-trivial stationary angles for amp_amp_op are as expected"""
    # Create uniform superposition of ancilla: qubit 0
    state_prep = tq.gates.H(0)

    # Define target state
    wfn_target = tq.QubitWaveFunction.from_array(asarray([- 1 / sqrt(2), 0, 1 / sqrt(2), 0]))

    # Define amplitude amplification circuit
    amp_amp = amp_amp_op(state_prep, 0)
    projector = tq.paulis.Projector(wfn=wfn_target)

    expval = tq.ExpectationValue(H=projector, U=amp_amp)

    p = tq.simulate(expval)

    # assert isclose(p, 1.0, abs_tol=0.001)
    return p


def test_amp_amp_change() -> None:
    """Test whether the amp_amp_op works as expected for non-stationary angles"""
    # TODO
    # Create uniform superposition of ancilla: qubit 0
    state_prep = tq.gates.H(0)

    # Define target state
    wfn_target = tq.QubitWaveFunction.from_array(asarray([1 / sqrt(2), 0, - 1 / sqrt(2), 0]))

    # lcu_object = LCU(ancilla=0, unitaries=[(1 / sqrt(2), tq.gates.X(1)),
    #                                        (1 / sqrt(2), tq.gates.Z(1))])
    # amp_amp, lcu_op = lcu_object.amp_amp, lcu_object.lcu

    # Define amplitude amplification circuit
    amp_amp = amp_amp_op(state_prep, 0)
    projector = tq.paulis.Projector(wfn=wfn_target)

    expval = tq.ExpectationValue(H=projector, U=amp_amp)

    p = tq.simulate(expval)

    return p


# if __name__ == '__main__':
#     import pytest
#     pytest.main()
