"""
CSC299 ROP
Testing file
"""

# Import statements
import tequila as tq
from math import sqrt
from lcu_v1 import LCU, amp_amp_op


# Implement testing functions

# TODO

############################################################################################
# Test 1 qubit ancilla
############################################################################################

# Test prepare

# Test select

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

def test_amp_amp_stationary() -> None:
    """Test whether the non-trivial stationary angles for amp_amp_op are as expected"""
    # TODO

    # Create uniform superposition of ancilla: qubit 0
    state_prep = tq.gates.H(0)

    lcu_object = LCU(ancilla=0, unitaries=[(1 / sqrt(2), tq.gates.X(1)),
                                           (1 / sqrt(2), tq.gates.Z(1))])
    amp_amp_op, lcu_op = lcu_object.amp_amp, lcu_object.lcu


def test_amp_amp_change() -> None:
    """Test whether the amp_amp_op works as expected for non-stationary angles"""
    # TODO
    raise NotImplementedError


if __name__ == '__main__':
    import pytest
    pytest.main()
