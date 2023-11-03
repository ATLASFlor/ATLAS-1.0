# <p align="center">ATLAS-1.0</p> #  
## Atmospheric Lagrangian Dispersion Model for tephra transport and deposition ##
### Reckziegel Florencia, Folch Arnau & Viramonte José. ###


## Introduction ##
ATLAS (Atmospheric Lagrangian dispersion) is an atmospheric dispersion model tailored to volcanic tephra/ash particles.

Using the Advection-Diiffusion-Sedimentation equation, ATLAS describes the plume dispersion and sedimentation, and calculates
load and concentrations isopach. 
The model can be run in differnt scales (regional to global), with different off-line numerical weather prediction models. 
ATLAS can be used in forward mode (to forecast the ash dispersal from a volcano or to describe a past eruption)
or backward mode(to constrain unknown source term characteristics). Multiple source terms can be defined, with different
granulometric characteritics, in the same simulation.

## Dependencies ##
ATLAS requires:

* An MPI implementation(OpenMPI, mpich)
* An MPI-compatible F90 compiler
* [The entire NCO toolkit](https://nco.sourceforge.net/)
* NetCDF-Fortran (Not a part of NCO tools)
* Python3

## Compiling ##
From the top level you can run

		make

to compile ATLAS. This will generate an executable called **P_ATLAS.1.0.x** under the Runs/ folder

## Running ##
ATLAS can be run by doing

		mpirun -n X ../P_ATLAS.1.0.x <simulation>

Where **mpirun** is the executable in your chosen MPI-implementation that runs an MPI executable, **simulation** is the name
of your simulation (See Manual for more details on folder structure and input files), and **X** is the number of concurrent instances 
you wish to run.
For testing purposes, an example is provided along with the necessary Data files, and can be run from Runs/pexample by doing

		mpirun -n X ../P_ATLAS_1.0.x pexample
