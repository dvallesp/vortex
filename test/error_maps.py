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

# modules in this package
from masclet_mock_files import *
from masclet_mock_grids import *
from masclet_mock_velocityfields import *

it = 4000
ncores = 8
verbose = True

use_mysolapst = True
ubox = [-5,5,-5,5,-5,5]
up_to_level = 4
plotname = 'sliceploterror_zoom2.png'

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
#solapst[0] = np.ones(truevx[0].shape, dtype='bool')

# test velocities
# first assert the written total velocities are the original ones
for ipatch in range(len(vx)):
    if ((vx[ipatch] - truevx[ipatch]) ** 2 + (vy[ipatch] - truevy[ipatch]) ** 2 + (
            vz[ipatch] - truevz[ipatch]) ** 2).sum() > 0:
        print('Problem with the original velocity in patch {}!'.format(ipatch))

velx = [vcx + vrx for vcx, vrx in zip(velcompx, velrotx)]
vely = [vcy + vry for vcy, vry in zip(velcompy, velroty)]
velz = [vcz + vrz for vcz, vrz in zip(velcompz, velrotz)]
v = [np.sqrt(vxi**2 + vyi**2 + vzi**2) for vxi, vyi, vzi in zip(vx, vy, vz)]

errx = [(vxi / vi) ** 2 * np.abs(vxitilde / vxi - 1) for vi, vxi, vxitilde in zip(v, vx, velx)]
erry = [(vyi / vi) ** 2 * np.abs(vyitilde / vyi - 1) for vi, vyi, vyitilde in zip(v, vy, vely)]
errz = [(vzi / vi) ** 2 * np.abs(vzitilde / vzi - 1) for vi, vzi, vzitilde in zip(v, vz, velz)]
err = [np.sqrt(errxi ** 2 + erryi ** 2 + errzi ** 2) for errxi, erryi, errzi in zip(errx, erry, errz)]

os.chdir(initialpath)

err = masclet.tools.clean_field(err, cr0amr, solapst, npatch, up_to_level=up_to_level)
u = masclet.tools.uniform_grid_zoom_parallel(err, ubox, up_to_level, npatch, patchnx, patchny, patchnz, patchrx,
                                             patchry, patchrz, size, nmax, ncores=ncores, copies=2, verbose=True)
bbox = np.array([[ubox[0], ubox[1]], [ubox[2], ubox[3]], [ubox[4], ubox[5]]])
data = dict(Error=u)
ds = yt.load_uniform_grid(data, u.shape, 3.08e24, bbox=bbox, nprocs=ncores)
slc = yt.SlicePlot(ds, 'z', 'Error')
slc.save(plotname)


if verbose:
    print('#####################')
    print('Error field is computed!')
    print('#####################')
