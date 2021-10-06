*************************************************************
## Parallelized Finite Difference Time Domain (FDTD) Solver
 - Version 3.0, October 04, 2021
 - Author(s):  Michael Strickland,
               Rafael L. Delgado,
               Sebastian Steinbeißer,
               Johanes H. Weber
 - Email:  mstrick6@kent.edu,
           rafael.delgado@upm.es,
           sebastian.steinbeisser@tum.de,
           dr.rer.nat.weber@gmail.com
*************************************************************

 DESCRIPTION
-------------------------------------------------------------

This code can use several techniques (finite differences and
Fast Fourier Transform) to solve the Schrödinger Equation in
imaginary time for an arbitrary 3d potential. It uses the MPI
(Message Passing Interface) standard. The lattice is equally
partitioned into slices along the `x` direction. Code can
extract ground state and first few excited state wavefunction
and energies.

The Relativistic Kinetic Term for the Schrodinger Equation
can be used by means of the Fast Fourier Transform (FFT),
with strictly periodic boundary conditions.

REQUIREMENTS
-------------------------------------------------------------

The MPI (Message Passing Interface) API must be installed
on your system. Currently tested against MPICH and OpenMPI.
Can run on a single computational node or as many as you
like.

The FFTW3-MPI and GSL (GNU Scientific Library) are also
required to be installed.

 COMPILING
-------------------------------------------------------------

To compile, simply type `make` from the command line.

 USAGE
-------------------------------------------------------------

To run:

```bash
mpirun -np <Number of Worker Nodes> mpisolve
```

All parameters are specified in the params.txt file. They
can also be set via the commandline using, e.g.,

```bash
mpirun -np <Number of Worker Nodes> mpisolve -PARAMNAME [value]
```

Parameters set via the commandline override those set in
the `params.txt` file.

 DEBUGGING
-------------------------------------------------------------

```bash
mpirun N -x DISPLAY run_gdb.csh mpisolve
```

 CONTRIBUTORS
-------------------------------------------------------------
Michael Strickland,
Adrian Dumitru,
Yun Guo,
Rafael L. Delgado,
Sebastian Steinbeißer, and
Johannes H. Weber

 LICENSE
-------------------------------------------------------------

GNU General Public License (GPLv3).
See detailed text in the [LICENSE.md](./LICENSE.md) file.

 ATTRIBUTION
-------------------------------------------------------------

We ask that if you use this code for work which results in a
publication that you cite the following papers:

Rafael L. Delgado, Sebastian Steinbeißer, Michael Strickland, and Johannes H. Weber,  
"The Relativistic Schrödinger Equation through FFTW3: An Extension of quantumfdtd",  
To be submitted.  
https://inspirehep.net/literature/1804287

M. Strickland, and D. Yager-Elorriaga,  
"A Parallel Algorithm for Solving the 3d Schrodinger Equation",  
Journal of Computational Physics 229 (2010), 6015.  
https://inspirehep.net/literature/817954
