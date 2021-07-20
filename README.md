# nonunitary-embeddings
Automatized implementations of non-unitary embeddings in quantum computers

## General information
The field of quantum computation evolves rapidly on the hardware as well as on the software side. Elevating a known classical problem to be executable on a quantum computer requires to reformulate it in a physically interpretative way, using unitary operations that determine the evolution of the quantum state, and hermitian operators that determine the measurements on the quantum system. Many operations are easy to construct for classical algorithms but their embedding into a unitary operation are often not intuitive. There are however techniques to embed those operations into unitary operations on potentially larger qubit systems, using so-called ancillary qubits. Here, the non-unitary operation is embedded into a unitary operation on a larger system.  The goal of this project is to implement generalized techniques that allow the embedding of non-unitary operations in an automatized way providing a trustworthy access to them for several potential fields of applications.

## Table of contents

### Complete PDF

* LaTeX file
* Complete PDF (Planned)

### Notebook reports:

* Report 1: Simple implementation with a 1 qubit ancilla
* Report 2: Implementation with a 2 qubit ancilla
* Report 3: Generalized implementation and testing (In progress)
* Report 4: Alternative implementations and extensions (Planned)
* Report 5: Applications of algorithms (Planned)

### Python files
* two_unitaries.py: 1 qubit ancilla
* two_ancilla.py: 2 qubit ancilla
* lcu_v1.py: general case and data class
* test.py: collection of all testing functions

## Requirements
tequila
	
## Setup
TODO

## Plan
TODO

## TODO list
TODO (lol, sorry)
