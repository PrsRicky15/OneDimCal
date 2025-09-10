# oneDimcal
# 1D Hamiltonian Solver on a Grid

A numerical solver for the **one-dimensional Schrödinger equation** using a grid-based Hamiltonian (finite-difference discretization).  
Useful for studying bound states, eigenvalues, and eigenfunctions of arbitrary 1D potentials.

---

## Theory

We solve the stationary Schrödinger equation:
$$
\hat{H}\psi(x) = E\psi(x), \quad 
\hat{H} = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x)
$$

---

## Features

- Build 1D Hamiltonians on a uniform grid
- Predefined potentials
- Dense diagonalization (small systems) and sparse eigensolvers (large systems)
- Expectation value calculations
- Outputs eigenvalues and eigenfunctions for visualization

---
