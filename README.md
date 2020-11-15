# ``vortex`` code
### Implements Helmholtz-Hodge decomposition for an AMR velocity field.
#### Developed by David Vallés-Pérez, Susana Planelles and Vicent Quilis

``vortex`` has been developed at the Departament d'Astronomia i Astrofísica of the Universitat de València, in the Computational Cosmology group. This project has been supported by the Spanish Ministerio de Ciencia e Innovación (MICINN, grants AYA2016-77237-C3-3-P and PID2019-107427GB-C33) and by the Generalitat Valenciana (grant PROMETEO/2019/071).


##### Repository organisation
The folder ``./src`` contains the source code, written in FORTRAN and parallellised according to the OpenMP standard.
The folder ``./test`` contains Python code to generate a set of idealised tests and to assess the performance of the code in such tests.

##### GENERAL CONSIDERATIONS

Besides the source code, the following files are needed:

1) ``vortex_parameters.dat``. This file dimensions the arrays, and therefore the code needs to be compiled when these parameters are changed.
2) ``vortex.dat``. This file contains runtime parameters. They can be changed once the code has been compiled.
3) Simulation data. By default, we read the simulation data in a folder simu_masclet, which contains the "gas" files. In order to use vortex on other code's outputs, the functions in ``reader.f`` (actual reader of the outputs) and ``nomfile.f`` (names of the simulation data files) need to be adapted.

The outputs of the code are written, by default, inside a folder ``output_files``. This folder needs to be created before running vortex. The output file will be saved as ``velocitiesXXXXX`` (XXXXX is the iteration number). This behaviour can be changed in ``nomfile.f``. As an additional safety measure, the code stops if a file with the same name is already in the folder.

The code is parallelised according to the OpenMP standard directives. For the code to run in parallel, it has to be compiled with the flag ``-fopenmp`` (gfortran), and the environment variable ``OMP_NUM_THREADS`` needs to be set to the number of cores we want to run the code with.

##### Code structure

The routines in this programme are distributed in the following files:

- ``vortex.f``: main program, contains the main execution workflow.
- ``poisson.f``: solve Poisson's equations, both for the coarse grid and for the refinement patches
- ``diff.f`` or ``diff_ho.f``: perform the finite differences of the velocity field and of the potentials (first one to 1st order centered derivatives; second file for higher order derivatives)
- ``interp.f``: linearly interpolate from coarser to finer grids
- ``grids.f``: build the base and AMR grids
- ``overlaps.f``: handle the overlaps between different refinement patches
- ``outliers.f``: find cells where the reconstruction is not reliable and interpolate their values from the parent, coarser cells
- ``boundaries.f``: special treatments for patches boundaries (currently deprecated)
- ``filter.f``: filter out turbulent motions using a multi-scale algorithm
- ``writer.f``: write the output files
- ``reader.f``: read the input data (may need to be adapted to read data from other simulation codes
- ``nomfile.f``: generates the filenames for I/O
