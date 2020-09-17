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
sys.path.append(os.path.join(cwd, 'test/'))
solapst_path = os.path.join(cwd, '../../TFM/Compute_solapst/solapst')
# sys.path.append(os.path.join(cwd, 'src/'))

import numpy as np
import masclet_framework as masclet
from cython_fortran_file import FortranFile as FF
import pickle
import yt
from scipy.stats import spearmanr

# modules in this package
from masclet_mock_files import *
from masclet_mock_grids import *
from masclet_mock_velocityfields import *

it = 4000
ncores = 8
verbose = True
use_mysolapst = False
percentile = 90

#######################################
# test
#######################################
if verbose:
    print('#####################')
    print('Starting with the test')
    print('#####################')

output_path = initialpath
maxl = 9
nmax = 128
size = 40

npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, patchry, patchrz, \
pare = masclet.read_masclet.read_grids(it, path=os.path.join(output_path, 'simu_masclet'),
                                       parameters_path=output_path, digits=5, read_general=False,
                                       read_patchnum=True, read_dmpartnum=False, read_patchcellextension=True,
                                       read_patchcellposition=True, read_patchposition=True, read_patchparent=True,
                                       nparray=True)

truevx, truevy, truevz, cr0amr = masclet.read_masclet.read_clus(it, path=os.path.join(output_path, 'simu_masclet'),
                                                                parameters_path=output_path, digits=5,
                                                                output_delta=False, output_pres=False, output_pot=False,
                                                                output_temp=False, output_metalicity=False,
                                                                output_cr0amr=True, output_solapst=False)
# for the usual solapst
solapst = masclet.read_masclet.read_npz_field('solapst_bool' + '{:05d}'.format(it), solapst_path)
# solapst[0] = np.ones(truevx[0].shape, dtype='bool')

div, rotx, roty, rotz, scalarpot, vecpotx, vecpoty, vecpotz, vx, vy, vz, velcompx, velcompy, velcompz, velrotx, \
velroty, velrotz, _ = masclet.read_masclet.read_vortex(it, path=os.path.join(output_path, 'output_files'),
                                                       grids_path=os.path.join(output_path, 'simu_masclet'),
                                                       parameters_path=output_path,
                                                       digits=5, are_divrot=True, are_potentials=True,
                                                       are_velocities=True,
                                                       is_solapst=True, verbose=False)

if use_mysolapst:
    solapst = _
else:
    del _
# solapst[0] = np.ones(truevx[0].shape, dtype='bool')

# test velocities
# first assert the written total velocities are the original ones
for ipatch in range(len(vx)):
    if ((vx[ipatch] - truevx[ipatch]) ** 2 + (vy[ipatch] - truevy[ipatch]) ** 2 + (
            vz[ipatch] - truevz[ipatch]) ** 2).sum() > 0:
        print('Problem with the original velocity in patch {}!'.format(ipatch))

velx = [vcx + vrx for vcx, vrx in zip(velcompx, velrotx)]
vely = [vcy + vry for vcy, vry in zip(velcompy, velroty)]
velz = [vcz + vrz for vcz, vrz in zip(velcompz, velrotz)]
v = [np.sqrt(vxi ** 2 + vyi ** 2 + vzi ** 2) for vxi, vyi, vzi in zip(vx, vy, vz)]

errx = [(vxi / vi) ** 2 * np.abs(vxitilde / vxi - 1) for vi, vxi, vxitilde in zip(v, vx, velx)]
erry = [(vyi / vi) ** 2 * np.abs(vyitilde / vyi - 1) for vi, vyi, vyitilde in zip(v, vy, vely)]
errz = [(vzi / vi) ** 2 * np.abs(vzitilde / vzi - 1) for vi, vzi, vzitilde in zip(v, vz, velz)]
err = [np.sqrt(errxi ** 2 + erryi ** 2 + errzi ** 2) for errxi, erryi, errzi in zip(errx, erry, errz)]

os.chdir(initialpath)

# errclean = masclet.tools.clean_field(err, cr0amr, solapst, npatch, up_to_level=up_to_level)

boundariesx = []
boundariesy = []
boundariesz = []
timeserrboundaries = []
corr_err_vmodulos = []
corr_err_vmincomps = []

for ipatch in range(npatch.sum() + 1):
    A = np.where(err[ipatch] > np.percentile(err[ipatch], percentile))
    Ax = A[0]
    Ay = A[1]
    Az = A[2]

    histAx = np.histogram(Ax, patchnx[ipatch] - 1)
    histAy = np.histogram(Ay, patchny[ipatch] - 1)
    histAz = np.histogram(Az, patchnz[ipatch] - 1)

    boundaryx = (histAx[0][0] + histAx[0][-1]) / histAx[0].sum()
    boundaryy = (histAy[0][0] + histAy[0][-1]) / histAy[0].sum()
    boundaryz = (histAz[0][0] + histAz[0][-1]) / histAz[0].sum()

    boundariesx.append(boundaryx)
    boundariesy.append(boundaryy)
    boundariesz.append(boundaryz)

    n1 = patchnx[ipatch]
    n2 = patchny[ipatch]
    n3 = patchnz[ipatch]

    errboundarycells = np.concatenate([err[ipatch][0, :, :].flatten(), err[ipatch][n1 - 1, :, :].flatten(),
                                       err[ipatch][1:n1, 0, :].flatten(), err[ipatch][1:n1, n2 - 1, :].flatten(),
                                       err[ipatch][1:n1, 1:n2, 0].flatten(), err[ipatch][1:n1, 1:n2, n3-1].flatten()])
    errboundarycells = np.percentile(errboundarycells, percentile)
    errinterior = err[ipatch][1:n1,1:n2,1:n3].flatten()
    errinterior = np.percentile(errinterior, percentile)
    timeserrboundary = errboundarycells / errinterior

    timeserrboundaries.append(timeserrboundary)

    corr_err_vmodulo = spearmanr(err[ipatch].flatten(), v[ipatch].flatten())
    corr_err_vmincomp = spearmanr(err[ipatch].flatten(),
                                    np.minimum(np.minimum(np.abs(vx[ipatch]), np.abs(vy[ipatch])),np.abs(vz[ipatch])).flatten())

    corr_err_vmodulos.append(corr_err_vmodulo)
    corr_err_vmincomps.append(corr_err_vmincomp)


for l in range(maxl + 1):
    if l == 0:
        minipatch = 0
        maxipatch = 0
    else:
        minipatch = npatch[0:l].sum() + 1
        maxipatch = npatch[0:l + 1].sum()
    if l == 0:
        print('Level \t Nerrboundx \t Nerrboundy \t Nerrboundz \t err_bdry / err_interior \t max(previous)')
    print(l, np.mean(boundariesx[minipatch:maxipatch + 1]),
          np.mean(boundariesy[minipatch:maxipatch + 1]), np.mean(boundariesz[minipatch:maxipatch + 1]),
          np.mean(timeserrboundaries[minipatch:maxipatch + 1]), np.max(timeserrboundaries[minipatch:maxipatch + 1]))

for l in range(maxl + 1):
    if l == 0:
        minipatch = 0
        maxipatch = 0
    else:
        minipatch = npatch[0:l].sum() + 1
        maxipatch = npatch[0:l + 1].sum()
    if l == 0:
        print('Level \t Corrcoef')
    print(l, np.mean(corr_err_vmodulos[minipatch:maxipatch + 1]), np.mean(corr_err_vmincomps[minipatch:maxipatch + 1]))

if verbose:
    print('#####################')
    print('Error field is computed!')
    print('#####################')
