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
from masclet_mock_files import write_mock_grids, write_mock_clus

mockit = 4000

output_path = os.path.join(initialpath, 'test_files')

ncores = 8


def mock_gridsA(nlevels, Nx, L=1):
    npatch = np.zeros(nlevels + 1, dtype='i4')
    npatch[1:] = 8
    numpatches = 8 * nlevels + 1
    patchnx = np.zeros(numpatches, dtype='i4')
    patchny = np.zeros(numpatches, dtype='i4')
    patchnz = np.zeros(numpatches, dtype='i4')
    patchx = np.zeros(numpatches, dtype='i4')
    patchy = np.zeros(numpatches, dtype='i4')
    patchz = np.zeros(numpatches, dtype='i4')
    patchrx = np.zeros(numpatches, dtype='f4')
    patchry = np.zeros(numpatches, dtype='f4')
    patchrz = np.zeros(numpatches, dtype='f4')
    pare = np.zeros(numpatches, dtype='i4')

    patchnx[0] = Nx
    patchny[0] = Nx
    patchnz[0] = Nx
    patchrx[0] = - L / 2 + L / Nx
    patchry[0] = - L / 2 + L / Nx
    patchrz[0] = - L / 2 + L / Nx

    ipatch = 0
    for l in range(1, nlevels + 1):
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    ipatch += 1

                    ipare = max(0, ipatch - 8)
                    pare[ipatch] = ipare

                    patchnx[ipatch] = Nx / 2
                    patchny[ipatch] = Nx / 2
                    patchnz[ipatch] = Nx / 2

                    patchx[ipatch] = (patchnx[ipatch] / 2 * (i + 1)) % patchnx[ipare]
                    patchy[ipatch] = (patchny[ipatch] / 2 * (j + 1)) % patchny[ipare]
                    patchz[ipatch] = (patchnz[ipatch] / 2 * (k + 1)) % patchnz[ipare]

                    patchrx[ipatch] = patchrx[ipare] + (patchx[ipatch] - 1 / 2) * L / Nx / 2 ** (l - 1)
                    patchry[ipatch] = patchry[ipare] + (patchy[ipatch] - 1 / 2) * L / Nx / 2 ** (l - 1)
                    patchrz[ipatch] = patchrz[ipare] + (patchz[ipatch] - 1 / 2) * L / Nx / 2 ** (l - 1)

    return npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, patchry, patchrz, pare

def constant_curl_field(cellsrx, cellsry, cellsrz, constant):
    cellsrho = [np.sqrt(cellrx ** 2 + cellry ** 2) for cellrx, cellry, cellrz in zip(cellsrx, cellsry, cellsrz)]
    v = [constant * rho for rho in cellsrho]
    vx = [vi * y / rho for vi, y, rho in zip(v, cellsry, cellsrho)]
    vy = [- vi * x / rho for vi, x, rho in zip(v, cellsrx, cellsrho)]
    vz = [np.zeros(vxi.shape) for vxi in vx]
    return vx, vy, vz

# test1b: constant curl field
maxl = 10
nmax = 128
size=1

npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, patchry, patchrz, pare = mock_gridsA(nlevels=maxl,
                                                                                                         Nx=nmax,
                                                                                                         L=size)

write_mock_grids(mockit, 100, maxl, 1e-10, 0.1, 100000, npatch, patchnx, patchny, patchnz, patchx, patchy, patchz,
                 patchrx, patchry, patchrz, pare, digits=5, path=output_path)

cellsrx, cellsry, cellsrz = masclet.tools_xyz.compute_position_fields(patchnx, patchny, patchnz, patchrx, patchry,
                                                                      patchrz, npatch, size, nmax, ncores=ncores)

vx, vy, vz = constant_curl_field(cellsrx, cellsry, cellsrz, 0.01)

write_mock_clus(vx, vy, vz, mockit, digits=5, path=output_path)

os.chdir(initialpath)
