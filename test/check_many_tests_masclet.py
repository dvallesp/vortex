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
solapst_path = '/scratch/Simulations/new_cluster/solapst_new/'
# sys.path.append(os.path.join(cwd, 'src/'))

import numpy as np
import masclet_framework as masclet
from cython_fortran_file import FortranFile as FF
import pickle
import gc

# modules in this package
from masclet_mock_files import *
from masclet_mock_grids import *
from masclet_mock_velocityfields import *

itini = 50
itfin = 4050
every = 50
ncores = 8
verbose = True
eps_err = 0.01

percentiles = [5, 25, 50, 75, 95]
use_solapst = True
file_dictionary = {it: {} for it in range(itini, itfin+every, every)}

#######################################
# test
#######################################
for it in range(itini, itfin+every, every):
    if verbose:
        print('#####################')
        print('Starting with the test for it {:05d}'.format(it))
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

    truevx, truevy, truevz = masclet.read_masclet.read_clus(it, path=os.path.join(output_path, 'simu_masclet'),
                                                            parameters_path=output_path, digits=5, output_delta=False,
                                                            output_pres=False, output_pot=False,
                                                            output_temp=False, output_metalicity=False, output_cr0amr=False,
                                                            output_solapst=False)
    # for the usual solapst
    solapst = masclet.read_masclet.read_npz_field('solapst_bool' + '{:05d}'.format(it), solapst_path)
    solapst[0] = np.ones(truevx[0].shape, dtype='bool')

    div, rotx, roty, rotz, scalarpot, vecpotx, vecpoty, vecpotz, vx, vy, vz, velcompx, velcompy, velcompz, velrotx, \
    velroty, velrotz = masclet.read_masclet.read_vortex(it, path=os.path.join(output_path, 'output_files'),
                                                        grids_path=os.path.join(output_path, 'simu_masclet'),
                                                        parameters_path=output_path, digits=5, are_divrot=True,
                                                        are_potentials=True, are_velocities=True, verbose=False)

    # test velocities

    ev = {}
    frac_vcomp = {}
    frac_vrot = {}

    maxvx = max([np.max(abs(vxi)) for vxi in vx])
    maxvy = max([np.max(abs(vyi)) for vyi in vy])
    maxvz = max([np.max(abs(vzi)) for vzi in vz])

    for l in range(maxl + 1):
        if l == 0:
            minipatch = 0
            maxipatch = 0
        else:
            minipatch = npatch[0:l].sum() + 1
            maxipatch = npatch[0:l + 1].sum()

        if verbose:
            print('Level {}, patches {} to {}'.format(l, minipatch, maxipatch))

        if maxipatch < minipatch:
            ev[l] = 0
            frac_vcomp[l] = 0
            frac_vrot[l] = 0
            continue

        if use_solapst:
            solapst_l = np.concatenate([a.flatten() for a in solapst[minipatch:maxipatch + 1]])
        else:
            solapst_l = np.concatenate([np.ones(a.shape, dtype='bool').flatten() for a in solapst[minipatch:maxipatch + 1]])

        vx_l = np.concatenate([a.flatten() for a in vx[minipatch:maxipatch + 1]])[solapst_l]
        vy_l = np.concatenate([a.flatten() for a in vy[minipatch:maxipatch + 1]])[solapst_l]
        vz_l = np.concatenate([a.flatten() for a in vz[minipatch:maxipatch + 1]])[solapst_l]
        v_l = np.sqrt(vx_l ** 2 + vy_l ** 2 + vz_l ** 2)

        velcompx_l = np.concatenate([a.flatten() for a in velcompx[minipatch:maxipatch + 1]])[solapst_l]
        velcompy_l = np.concatenate([a.flatten() for a in velcompy[minipatch:maxipatch + 1]])[solapst_l]
        velcompz_l = np.concatenate([a.flatten() for a in velcompz[minipatch:maxipatch + 1]])[solapst_l]
        velcomp_l = np.sqrt(velcompx_l ** 2 + velcompy_l ** 2 + velcompz_l ** 2)

        velrotx_l = np.concatenate([a.flatten() for a in velrotx[minipatch:maxipatch + 1]])[solapst_l]
        velroty_l = np.concatenate([a.flatten() for a in velroty[minipatch:maxipatch + 1]])[solapst_l]
        velrotz_l = np.concatenate([a.flatten() for a in velrotz[minipatch:maxipatch + 1]])[solapst_l]
        velrot_l = np.sqrt(velrotx_l ** 2 + velroty_l ** 2 + velrotz_l ** 2)

        velrecx_l = velcompx_l + velrotx_l
        velrecy_l = velcompy_l + velroty_l
        velrecz_l = velcompz_l + velrotz_l

        ev[l] = [np.percentile(np.sqrt(((vx_l / v_l) ** 2 * np.abs(velrecx_l - vx_l) / (abs(vx_l) + eps_err * maxvx)) ** 2 +
                                       ((vy_l / v_l) ** 2 * np.abs(velrecy_l - vy_l) / (abs(vy_l) + eps_err * maxvy)) ** 2 +
                                       ((vz_l / v_l) ** 2 * np.abs(velrecz_l - vz_l) / (abs(vz_l) + eps_err * maxvz)) ** 2),
                               perc) for perc in percentiles]
        frac_vcomp[l] = [np.percentile(velcomp_l / v_l, perc) for perc in percentiles]
        frac_vrot[l] = [np.percentile(velrot_l / v_l, perc) for perc in percentiles]

        if verbose:
            print('* At level {}, percentiles: '.format(l), percentiles)
            print('Relative errors in velocity reconstruction: ', ['{:.2e}'.format(i) for i in ev[l]])
            print('Fraction of compressive velocity: ', ['{:.2e}'.format(i) for i in frac_vcomp[l]])
            print('Fraction of rotational velocity: ', ['{:.2e}'.format(i) for i in frac_vrot[l]], '\n')

    os.chdir(initialpath)
    gc.collect()

    file_dictionary[it] = {'percentiles': percentiles, 'ev': ev, 'frac_vcomp': frac_vcomp, 'frac_vrot': frac_vrot}

if verbose:
    print('#####################')
    print('Done with the test!')
    print('#####################')

with open('test_results.pickle', 'wb') as f:
    pickle.dump(file_dictionary, f)

#### END MESSAGE
if verbose:
    print('Done!')
