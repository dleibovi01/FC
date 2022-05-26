# Fourier Continuation method

This repository contains routines and scripts associated to the Fourier Continuation procedure.

The FC_Data folder contains the Fourier Continuation matrices mentioned in equations (26) and (27) in [1],
which are used to generate the smooth periodic continuations.

(The functions in generate_bdry_continuations.m and precomp_fc_data.m were used to generate these matrices.)

The script advection_eqn can be run to study the convergence properties of an FC-based solver, using an
RK4 scheme in time and FC-based differentiation in space, for the case of a 1D Linear advection equation
with a prescribed Dirichlet boundary condition.

The script wave_eqn can be run to study the convergence properties of an FC-based solver, using an
AB4 scheme in time and FC-based differentiation in space, for the case of a 1D second order wave equation
with the option to prescribe either Dirichlet or Neumann boundary conditions at each domain boundary.


Ref : 
[1] Amlani, F., & Bruno, O. P. (2016). An FC-based spectral solver for elastodynamic problems in general 
    three-dimensional domains.  Journal of Computational Physics, 307, 333-354.
    
Special thanks to Thomas Anderson and Jagabandhu Paul.
