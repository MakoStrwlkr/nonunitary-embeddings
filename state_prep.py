"""
CSC299 ROP
State preparation module
"""

import tequila as tq
import numpy as np
import random
import time
from typing import Optional


def param_circ(anc: list) -> tq.QCircuit:
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

    return pqc


def unitary_gate(th, phi, lam, q0) -> tq.QCircuit:
    """Return a particular quantum gate that is not included in the basic gate set"""
    ugate = tq.gates.Rz(target=q0, angle=phi) + tq.gates.Ry(target=q0, angle=th) + tq.gates.Rz(
        target=q0, angle=lam)
    return ugate


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


def controlled_phase(phi, q0, q1) -> tq.QCircuit:
    """Return a quantum circuit that performs the controlled-phase shift operation, with
    q0 as control and q1 as target."""
    cph = tq.gates.Rz(target=q0, angle=phi/2) + tq.gates.CRz(control=q0, target=q1, angle=phi/2)
    return cph


# if __name__ == '__main__':
#     import doctest
#     doctest.testmod(verbose=True)
