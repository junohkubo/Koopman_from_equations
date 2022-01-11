# Numerical methods to evaluate Koopman matrix from system equations

Jun Ohkubo

[arXiv: https://arxiv.org/abs/2111.07213](https://arxiv.org/abs/2111.07213)

## Brief introduction

Recently, methods to obtain an approximation of the Koopman operator, i.e., a Koopman matrix, have been developed. The extended dynamical mode decomposition (EDMD) is one of the famous methods, in which a data set of snapshot pairs is used.

Recent developments in the duality for stochastic processes enable us to evaluate the Koopman matrix directly from nonlinear system equations. Without any sampling, a method based on combinatorics, an approximation of the resolvent, and extrapolations, is proposed. 

This approach is complementary to a data-driven approach provided a prior knowledge of the system equations is available. Here, some codes for the method are provided.

1. `dual_computation.c` (CODE 1): A code to evaluate an element of the Koopman matrix (C language)
2. `Makefile`: Makefile for CODE 1
3. `dual_derivation_van_der_Pol.ipynb` (CODE 2): A sample code to derive an event file for CODE 1 (Jupyter Notebook)
4. `van_der_Pol_parameters.inp`: A sample parameter file
5. `EDMD_van_der_Pol.ipynb` (CODE 3): A code for the EDMD, which is avaiable for the comparison.

CODE 3 and file 4 are for the noisy van der Pol system. It is easy to modify CODE 3 to your system.

## Usage (for sample codes)

First, perform CODE 3. Then, a file `van_der_Pol_event.inp` is generated. (`sympy` and some other modules are necessary.)

Second, compile CODE 1 by
```
./make
```

Third, perform CODE 1 as
```
./app_dual_computation van_der_Pol_event.inp van_der_Pol_parameters.inp
```
