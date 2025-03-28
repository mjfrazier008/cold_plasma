# cold_plasma
MATLAB code for spectral calculations of cold plasma interface Hamiltonians. For background see  	
[Topological edge states of continuous Hamiltonians](https://doi.org/10.48550/arXiv.2503.11811) and Appendix E in particular for description of numerical scheme implimented here. 

## spectrum.m:
Main file which generates and plots the spectrum of a user-specified edge Hamiltonain. Set up to utilize MATLAB parallel computing toolbox but can easily be modified to use without (change 'parfor' to 'for' in line 27), though not recommended. 

### system parameters:

$k_z$ (kz): vertical wave number

$\Omega$ (Om): Cyclotron frequency

$\omega_p$ (op): Plasma frequency

### discretization and output parameters:

L: 1/2 Length (centered around 0) of domain

N: Number of discretization points

kmin: Minimum value of $k_x/k_y$

kmax: Maximum value of $k_x/k_y$. 

Nk: Number of linearly spaced points between kmin and kmax to diagonalize $H_I$ at. 

## spectrum.m:
Diagonalizes edge Hamiltonian for only transverse magnetic Hamiltonian ($k_z = 0$). 

## discH_...  :
Functions which generate discretized interface Hamiltonians for diagonalization.

## make_H1, make_Hreg:
Generate constant coefficient Hamiltonians for cold plasma model. 