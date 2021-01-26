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
import pickle

# modules in this package
from masclet_mock_files import *
from masclet_mock_grids import *
from masclet_mock_velocityfields import *

mockit = 4000
ncores = 8
verbose = True
write = True

eps_err = 0.01
percentiles = [5, 25, 50, 75, 95]
weighted = False
file_dictionary = {}

def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

#######################################
# test1d: multifrequency field
#######################################
if verbose:
    print('Starting with test 1d')

output_path = os.path.join(initialpath, 'test_files/1d/')
maxl = 4
nmax = 128
size = 1
uni_nmax = nmax * 2 ** maxl

npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, patchry, patchrz, \
pare = masclet.read_masclet.read_grids(mockit, path=os.path.join(output_path, 'simu_masclet'),
                                       parameters_path=output_path, digits=5, read_general=False,
                                       read_patchnum=True, read_dmpartnum=False, read_patchcellextension=True,
                                       read_patchcellposition=True, read_patchposition=True, read_patchparent=True,
                                       nparray=True)

truevx, truevy, truevz = masclet.read_masclet.read_clus(mockit, path=os.path.join(output_path, 'simu_masclet'),
                                                        parameters_path=output_path, digits=5, output_delta=False,
                                                        output_pres=False, output_pot=False,
                                                        output_temp=False, output_metalicity=False, output_cr0amr=False,
                                                        output_solapst=False)

div, rotx, roty, rotz, scalarpot, vecpotx, vecpoty, vecpotz, vx, vy, vz, velcompx, velcompy, velcompz, velrotx, \
velroty, velrotz = masclet.read_masclet.read_vortex(mockit, path=os.path.join(output_path, 'output_files'),
                                                    grids_path=os.path.join(output_path, 'simu_masclet'),
                                                    parameters_path=output_path, digits=5, are_divrot=True,
                                                    are_potentials=True, are_velocities=True, verbose=False)

truevelcompx, truevelcompy, \
truevelcompz = masclet.read_masclet.read_clus(mockit, path=os.path.join(output_path,'simu_masclet/comp'),
                                              parameters_path=output_path, digits=5, output_delta=False,
                                              output_pres=False, output_pot=False, output_temp=False,
                                              output_metalicity=False, output_cr0amr=False, output_solapst=False)

truevelrotx, truevelroty, \
truevelrotz = masclet.read_masclet.read_clus(mockit, path=os.path.join(output_path,'simu_masclet/sol'),
                                              parameters_path=output_path, digits=5, output_delta=False,
                                              output_pres=False, output_pot=False, output_temp=False,
                                              output_metalicity=False, output_cr0amr=False, output_solapst=False)
# test velocities
# first assert the written total velocities are the original ones
for ipatch in range(len(vx)):
    if ((vx[ipatch] - truevx[ipatch]) ** 2 + (vy[ipatch] - truevy[ipatch]) ** 2 + (
            vz[ipatch] - truevz[ipatch]) ** 2).sum() > 0:
        print('Problem with the original velocity in patch {}!'.format(ipatch))

ev = {}
ev_comp = {}
ev_rot = {}
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
        minipatch = 8 * (l - 1) + 1
        maxipatch = 8 * l

    vx_l = np.concatenate([a.flatten() for a in vx[minipatch:maxipatch + 1]])
    vy_l = np.concatenate([a.flatten() for a in vy[minipatch:maxipatch + 1]])
    vz_l = np.concatenate([a.flatten() for a in vz[minipatch:maxipatch + 1]])
    v_l = np.sqrt(vx_l ** 2 + vy_l ** 2 + vz_l ** 2)

    velcompx_l = np.concatenate([a.flatten() for a in velcompx[minipatch:maxipatch + 1]])
    velcompy_l = np.concatenate([a.flatten() for a in velcompy[minipatch:maxipatch + 1]])
    velcompz_l = np.concatenate([a.flatten() for a in velcompz[minipatch:maxipatch + 1]])
    velcomp_l = np.sqrt(velcompx_l ** 2 + velcompy_l ** 2 + velcompz_l ** 2)

    velrotx_l = np.concatenate([a.flatten() for a in velrotx[minipatch:maxipatch + 1]])
    velroty_l = np.concatenate([a.flatten() for a in velroty[minipatch:maxipatch + 1]])
    velrotz_l = np.concatenate([a.flatten() for a in velrotz[minipatch:maxipatch + 1]])
    velrot_l = np.sqrt(velrotx_l ** 2 + velroty_l ** 2 + velrotz_l ** 2)

    truevelcompx_l = np.concatenate([a.flatten() for a in truevelcompx[minipatch:maxipatch + 1]])
    truevelcompy_l = np.concatenate([a.flatten() for a in truevelcompy[minipatch:maxipatch + 1]])
    truevelcompz_l = np.concatenate([a.flatten() for a in truevelcompz[minipatch:maxipatch + 1]])
    truevelcomp_l = np.sqrt(truevelcompx_l ** 2 + truevelcompy_l ** 2 + truevelcompz_l ** 2)

    truevelrotx_l = np.concatenate([a.flatten() for a in truevelrotx[minipatch:maxipatch + 1]])
    truevelroty_l = np.concatenate([a.flatten() for a in truevelroty[minipatch:maxipatch + 1]])
    truevelrotz_l = np.concatenate([a.flatten() for a in truevelrotz[minipatch:maxipatch + 1]])
    truevelrot_l = np.sqrt(truevelrotx_l ** 2 + truevelroty_l ** 2 + truevelrotz_l ** 2)

    velrecx_l = velcompx_l + velrotx_l
    velrecy_l = velcompy_l + velroty_l
    velrecz_l = velcompz_l + velrotz_l

    if not weighted:
        ev[l] = [np.percentile(np.sqrt(((vx_l / v_l) ** 2 * np.abs(velrecx_l - vx_l) / (abs(vx_l) + eps_err * maxvx)) ** 2 +
                                       ((vy_l / v_l) ** 2 * np.abs(velrecy_l - vy_l) / (abs(vy_l) + eps_err * maxvy)) ** 2 +
                                       ((vz_l / v_l) ** 2 * np.abs(velrecz_l - vz_l) / (abs(vz_l) + eps_err * maxvz)) ** 2),
                               perc) for perc in percentiles]
        ev_comp[l] = [np.percentile(np.sqrt(((truevelcompx_l / truevelcomp_l) ** 2 * np.abs(
            velcompx_l - truevelcompx_l) / (abs(truevelcompx_l) + eps_err * maxvx)) ** 2 +
                                            ((truevelcompy_l / truevelcomp_l) ** 2 * np.abs(
                                                velcompy_l - truevelcompy_l) / (
                                                         abs(truevelcompy_l) + eps_err * maxvy)) ** 2 +
                                            ((truevelcompz_l / truevelcomp_l) ** 2 * np.abs(
                                                velcompz_l - truevelcompz_l) / (
                                                         abs(truevelcompz_l) + eps_err * maxvz)) ** 2),
                                    perc) for perc in percentiles]
        ev_rot[l] = [np.percentile(np.sqrt(((truevelrotx_l / truevelrot_l) ** 2 * np.abs(velrotx_l - truevelrotx_l) / (
                    abs(truevelrotx_l) + eps_err * maxvx)) ** 2 +
                                           ((truevelroty_l / truevelrot_l) ** 2 * np.abs(velroty_l - truevelroty_l) / (
                                                       abs(truevelroty_l) + eps_err * maxvy)) ** 2 +
                                           ((truevelrotz_l / truevelrot_l) ** 2 * np.abs(velrotz_l - truevelrotz_l) / (
                                                       abs(truevelrotz_l) + eps_err * maxvz)) ** 2),
                                   perc) for perc in percentiles]
    else:
        ev[l] = [weighted_quantile(np.sqrt(((vx_l / v_l) ** 2 * np.abs(velrecx_l - vx_l) / (abs(vx_l)*(1-eps_err) + eps_err * maxvx)) ** 2 +
                                       ((vy_l / v_l) ** 2 * np.abs(velrecy_l - vy_l) / (abs(vy_l)*(1-eps_err) + eps_err * maxvy)) ** 2 +
                                       ((vz_l / v_l) ** 2 * np.abs(velrecz_l - vz_l) / (abs(vz_l)*(1-eps_err) + eps_err * maxvz)) ** 2),
                               perc/100, sample_weight=v_l**2, old_style=True) for perc in percentiles]

        ev_comp[l] = [weighted_quantile(np.sqrt(((truevelcompx_l / truevelcomp_l) ** 2 * np.abs(velcompx_l - truevelcompx_l) / (
                    abs(truevelcompx_l) + eps_err * maxvx)) ** 2 +
                                            ((truevelcompy_l / truevelcomp_l) ** 2 * np.abs(velcompy_l - truevelcompy_l) / (
                                                        abs(truevelcompy_l) + eps_err * maxvy)) ** 2 +
                                            ((truevelcompz_l / truevelcomp_l) ** 2 * np.abs(velcompz_l - truevelcompz_l) / (
                                                        abs(truevelcompz_l) + eps_err * maxvz)) ** 2),
                                    perc/100, sample_weight=truevelcomp_l**2, old_style=True) for perc in percentiles]

        ev_rot[l] = [weighted_quantile(np.sqrt(((truevelrotx_l / truevelrot_l) ** 2 * np.abs(velrotx_l - truevelrotx_l) / (
                    abs(truevelrotx_l) + eps_err * maxvx)) ** 2 +
                                           ((truevelroty_l / truevelrot_l) ** 2 * np.abs(velroty_l - truevelroty_l) / (
                                                       abs(truevelroty_l) + eps_err * maxvy)) ** 2 +
                                           ((truevelrotz_l / truevelrot_l) ** 2 * np.abs(velrotz_l - truevelrotz_l) / (
                                                       abs(truevelrotz_l) + eps_err * maxvz)) ** 2),
                                   perc/100, sample_weight=truevelrot_l**2, old_style=True) for perc in percentiles]
    frac_vcomp[l] = [np.percentile(velcomp_l / v_l, perc) for perc in percentiles]
    frac_vrot[l] = [np.percentile(velrot_l / v_l, perc) for perc in percentiles]

    print('* At level {}, percentiles: '.format(l), percentiles)
    print('Relative errors in velocity reconstruction: ', ['{:.2e}'.format(i) for i in ev[l]])
    print('Relative errors in compressive velocity reconstruction: ', ['{:.2e}'.format(i) for i in ev_comp[l]])
    print('Relative errors in solenoidal velocity reconstruction: ', ['{:.2e}'.format(i) for i in ev_rot[l]])
    print('Fraction of compressive velocity: ', ['{:.2e}'.format(i) for i in frac_vcomp[l]])
    print('Fraction of rotational velocity: ', ['{:.2e}'.format(i) for i in frac_vrot[l]], '\n')

os.chdir(initialpath)

file_dictionary['1d'] = {'percentiles': percentiles,
                         'ev_comp': ev_comp, 'ev_rot': ev_rot, 'ev': ev, 'frac_vcomp': frac_vcomp, 'frac_vrot': frac_vrot}

if verbose:
    print('#####################')
    print('Done with test 1d!')
    print('#####################')

if write:
    with open('test_results_1d.pickle', 'wb') as f:
        pickle.dump(file_dictionary, f)

#### END MESSAGE
if verbose:
    print('Done!')
