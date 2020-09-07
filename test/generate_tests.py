################################################################################
#                          GENERATE TESTS FOR VORTEX                           #
#                                                                              #
#                                                                              #
#                           Author: David Vallés                               #
################################################################################

import os, sys

initialpath = os.getcwd()
os.chdir('../')  # go to the root of the project
cwd = os.getcwd()
# sys.path.append(os.path.join(cwd, 'src/'))

import numpy as np
import masclet_framework as masclet
from cython_fortran_file import FortranFile as FF

# modules in this package
from masclet_mock_files import *
from masclet_mock_grids import *
from masclet_mock_velocityfields import *

mockit = 4000
ncores = 8
verbose = True

#######################################
# test1a: constant divergence field
#######################################
if verbose:
    print('Starting with test 1a')

output_path = os.path.join(initialpath, 'test_files/1a/simu_masclet')
maxl = 10
nmax = 128
size = 1

if verbose:
    print('Grids')

npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, \
    patchry, patchrz, pare = mock_gridsAprime(nlevels=maxl, Nx=nmax, L=size)

write_mock_grids(mockit, 100, maxl, 1e-10, 0.1, 100000, npatch, patchnx, patchny, patchnz, patchx, patchy, patchz,
                 patchrx, patchry, patchrz, pare, digits=5, path=output_path)

if verbose:
    print('Velocities')

cellsrx, cellsry, cellsrz = masclet.tools_xyz.compute_position_fields(patchnx, patchny, patchnz, patchrx, patchry,
                                                                      patchrz, npatch, size, nmax, ncores=ncores)

vx, vy, vz = constant_div_field(cellsrx, cellsry, cellsrz, 0.01)

write_mock_clus(vx, vy, vz, mockit, digits=5, path=output_path)

os.chdir(initialpath)

if verbose:
    print('Done!')

#######################################
# test1b: constant curl field
#######################################
if verbose:
    print('Starting with test 1b')

output_path = os.path.join(initialpath, 'test_files/1b/simu_masclet')
maxl = 10
nmax = 128
size = 1

if verbose:
    print('Grids')

npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, \
    patchry, patchrz, pare = mock_gridsAprime(nlevels=maxl, Nx=nmax, L=size)

write_mock_grids(mockit, 100, maxl, 1e-10, 0.1, 100000, npatch, patchnx, patchny, patchnz, patchx, patchy, patchz,
                 patchrx, patchry, patchrz, pare, digits=5, path=output_path)

if verbose:
    print('Velocities')

cellsrx, cellsry, cellsrz = masclet.tools_xyz.compute_position_fields(patchnx, patchny, patchnz, patchrx, patchry,
                                                                      patchrz, npatch, size, nmax, ncores=ncores)

vx, vy, vz = constant_curl_field(cellsrx, cellsry, cellsrz, 0.01)

write_mock_clus(vx, vy, vz, mockit, digits=5, path=output_path)

os.chdir(initialpath)


#######################################
# test1c: constant divergence field in realistic grid structure
#######################################
if verbose:
    print('Starting with test 1c')

output_path = os.path.join(initialpath, 'test_files/1c/simu_masclet')
maxl = 9
nmax = 128
size=1

if verbose:
    print('Grids')

## read grids from whim simulation
npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, patchry, patchrz, \
    pare = masclet.read_masclet.read_grids(4000, path='../../../TFM/Accretion/data/simulation_outputs',
                      parameters_path='../../../TFM/Accretion/data', digits=5, read_general=False, read_patchnum=True,
                      read_dmpartnum=False, read_patchcellextension=True, read_patchcellposition=True,
                      read_patchposition=True, read_patchparent=True, nparray=True)
patchrx = patchrx / 40
patchry = patchry / 40
patchrz = patchrz / 40
patchx = patchx+1
patchy = patchy+1
patchz = patchz+1


write_mock_grids(mockit, 100, maxl, 1e-10, 0.1, 100000, npatch, patchnx, patchny, patchnz, patchx, patchy, patchz,
                 patchrx, patchry, patchrz, pare, digits=5, path=output_path)

if verbose:
    print('Velocities')

cellsrx, cellsry, cellsrz = masclet.tools_xyz.compute_position_fields(patchnx, patchny, patchnz, patchrx, patchry,
                                                                      patchrz, npatch, size, nmax, ncores=ncores)

vx, vy, vz = constant_div_field(cellsrx, cellsry, cellsrz, 0.01)

write_mock_clus(vx, vy, vz, mockit, digits=5, path=output_path)

os.chdir(initialpath)

#### END MESSAGE
if verbose:
    print('Done!')