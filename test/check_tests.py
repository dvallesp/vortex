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

percentiles = [5, 25, 50, 75, 95]

#######################################
# test1a: constant divergence field
#######################################
if verbose:
    print('#####################')
    print('Starting with test 1a')
    print('#####################')

output_path = os.path.join(initialpath, 'test_files/1a')
maxl = 10
nmax = 128
size = 1

expected_div = 0.03
expected_vmax = 0.01

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
                                   grids_path=os.path.join(output_path, 'simu_masclet'), parameters_path=output_path,
                                   digits=5, are_divrot=True, are_potentials=True, are_velocities=True, verbose=False)

# test divergence
ediv = {}
erot = {}
for l in range(maxl+1):
    Dx = size/nmax/2**l

    if l == 0:
        minipatch = 0
        maxipatch = 0
    else:
        minipatch = 8 * (l-1) + 1
        maxipatch = 8 * l

    divs = np.concatenate([a.flatten() for a in div[minipatch:maxipatch+1]])
    rotsx = np.concatenate([a.flatten() for a in rotx[minipatch:maxipatch+1]])
    rotsy = np.concatenate([a.flatten() for a in roty[minipatch:maxipatch+1]])
    rotsz = np.concatenate([a.flatten() for a in rotz[minipatch:maxipatch+1]])
    rots = np.sqrt(rotsx ** 2 + rotsy ** 2 + rotsz ** 2)

    ediv[l] = [np.percentile(np.abs(divs/expected_div - 1), perc) for perc in percentiles]
    erot[l] = [np.percentile(rots, perc) * Dx / expected_vmax for perc in percentiles]

    print('* At level {}, percentiles: '.format(l), percentiles)
    print('Relative errors in divergence: ', ['{:.2e}'.format(i) for i in ediv[l]])
    print('Errors in rotational: ', ['{:.2e}'.format(i) for i in erot[l]], '\n')

# test velocities
# first assert the written total velocities are the original ones
for ipatch in range(len(vx)):
    if ((vx[ipatch]-truevx[ipatch])**2 + (vy[ipatch]-truevy[ipatch])**2 + (vz[ipatch]-truevz[ipatch])**2).sum() > 0:
        print('Problem with the original velocity in patch {}!'.format(ipatch))

ev = {}
frac_vcomp = {}
frac_vrot = {}
for l in range(maxl+1):
    if l == 0:
        minipatch = 0
        maxipatch = 0
    else:
        minipatch = 8 * (l-1) + 1
        maxipatch = 8 * l

    vx_l = np.concatenate([a.flatten() for a in vx[minipatch:maxipatch+1]])
    vy_l = np.concatenate([a.flatten() for a in vy[minipatch:maxipatch+1]])
    vz_l = np.concatenate([a.flatten() for a in vz[minipatch:maxipatch+1]])
    v_l = np.sqrt(vx_l ** 2 + vy_l ** 2 + vz_l ** 2)

    velcompx_l = np.concatenate([a.flatten() for a in velcompx[minipatch:maxipatch+1]])
    velcompy_l = np.concatenate([a.flatten() for a in velcompy[minipatch:maxipatch+1]])
    velcompz_l = np.concatenate([a.flatten() for a in velcompz[minipatch:maxipatch+1]])
    velcomp_l = np.sqrt(velcompx_l ** 2 + velcompy_l ** 2 + velcompz_l ** 2)

    velrotx_l = np.concatenate([a.flatten() for a in velrotx[minipatch:maxipatch+1]])
    velroty_l = np.concatenate([a.flatten() for a in velroty[minipatch:maxipatch+1]])
    velrotz_l = np.concatenate([a.flatten() for a in velrotz[minipatch:maxipatch+1]])
    velrot_l = np.sqrt(velrotx_l ** 2 + velroty_l ** 2 + velrotz_l ** 2)

    velrecx_l = velcompx_l + velrotx_l
    velrecy_l = velcompy_l + velroty_l
    velrecz_l = velcompz_l + velrotz_l

    ev[l] = [np.mean([np.percentile(np.abs(velrecx_l/vx_l - 1), perc),
                     np.percentile(np.abs(velrecy_l/vy_l - 1), perc),
                     np.percentile(np.abs(velrecz_l/vz_l - 1), perc)]) for perc in percentiles]
    frac_vcomp[l] = [np.percentile(velcomp_l / v_l, perc) for perc in percentiles]
    frac_vrot[l] = [np.percentile(velrot_l / v_l, perc) for perc in percentiles]

    print('* At level {}, percentiles: '.format(l), percentiles)
    print('Relative errors in velocity reconstruction: ', ['{:.2e}'.format(i) for i in ev[l]])
    print('Fraction of compressive velocity: ', ['{:.2e}'.format(i) for i in frac_vcomp[l]])
    print('Fraction of rotational velocity: ', ['{:.2e}'.format(i) for i in frac_vrot[l]], '\n')


os.chdir(initialpath)

if verbose:
    print('#####################')
    print('Done with test 1a!')
    print('#####################')

#######################################
# test1b: constant curl field
#######################################
if verbose:
    print('#####################')
    print('Starting with test 1b')
    print('#####################')

output_path = os.path.join(initialpath, 'test_files/1b')
maxl = 10
nmax = 128
size = 1

expected_rot = -0.02
expected_vmax = 0.01

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
                                   grids_path=os.path.join(output_path, 'simu_masclet'), parameters_path=output_path,
                                   digits=5, are_divrot=True, are_potentials=True, are_velocities=True, verbose=False)

# test rotational
ediv = {}
erotxy = {}
erotz = {}
for l in range(maxl+1):
    Dx = size/nmax/2**l

    if l == 0:
        minipatch = 0
        maxipatch = 0
    else:
        minipatch = 8 * (l-1) + 1
        maxipatch = 8 * l

    divs = np.concatenate([a.flatten() for a in div[minipatch:maxipatch+1]])
    rotsx = np.concatenate([a.flatten() for a in rotx[minipatch:maxipatch+1]])
    rotsy = np.concatenate([a.flatten() for a in roty[minipatch:maxipatch+1]])
    rotsz = np.concatenate([a.flatten() for a in rotz[minipatch:maxipatch+1]])
    rotsxy = np.sqrt(rotsx ** 2 + rotsy ** 2)

    ediv[l] = [np.percentile(divs, perc) * Dx / expected_vmax for perc in percentiles]
    erotz[l] = [np.percentile(np.abs(rotsz/expected_rot - 1), perc) for perc in percentiles]
    erotxy[l] = [np.percentile(divs, perc) * Dx / expected_vmax for perc in percentiles]

    print('* At level {}, percentiles: '.format(l), percentiles)
    print('Errors in divergence: ', ['{:.2e}'.format(i) for i in ediv[l]])
    print('Relative errors in rotational-z: ', ['{:.2e}'.format(i) for i in erotz[l]])
    print('Errors in rotational-xy: ', ['{:.2e}'.format(i) for i in erotxy[l]], '\n')

# test velocities
# first assert the written total velocities are the original ones
for ipatch in range(len(vx)):
    if ((vx[ipatch]-truevx[ipatch])**2 + (vy[ipatch]-truevy[ipatch])**2 + (vz[ipatch]-truevz[ipatch])**2).sum() > 0:
        print('Problem with the original velocity in patch {}!'.format(ipatch))

ev = {}
frac_vcomp = {}
frac_vrot = {}
for l in range(maxl+1):
    if l == 0:
        minipatch = 0
        maxipatch = 0
    else:
        minipatch = 8 * (l-1) + 1
        maxipatch = 8 * l

    vx_l = np.concatenate([a.flatten() for a in vx[minipatch:maxipatch+1]])
    vy_l = np.concatenate([a.flatten() for a in vy[minipatch:maxipatch+1]])
    v_l = np.sqrt(vx_l ** 2 + vy_l ** 2)

    velcompx_l = np.concatenate([a.flatten() for a in velcompx[minipatch:maxipatch+1]])
    velcompy_l = np.concatenate([a.flatten() for a in velcompy[minipatch:maxipatch+1]])
    velcomp_l = np.sqrt(velcompx_l ** 2 + velcompy_l ** 2)

    velrotx_l = np.concatenate([a.flatten() for a in velrotx[minipatch:maxipatch+1]])
    velroty_l = np.concatenate([a.flatten() for a in velroty[minipatch:maxipatch+1]])
    velrot_l = np.sqrt(velrotx_l ** 2 + velroty_l ** 2)

    velrecx_l = velcompx_l + velrotx_l
    velrecy_l = velcompy_l + velroty_l
    velrecz_l = velcompz_l + velrotz_l

    ev[l] = [np.mean([np.percentile(np.abs(velrecx_l/vx_l - 1), perc),
                     np.percentile(np.abs(velrecy_l/vy_l - 1), perc)]) for perc in percentiles]
    frac_vcomp[l] = [np.percentile(velcomp_l / v_l, perc) for perc in percentiles]
    frac_vrot[l] = [np.percentile(velrot_l / v_l, perc) for perc in percentiles]

    print('* At level {}, percentiles: '.format(l), percentiles)
    print('Relative errors in velocity reconstruction: ', ['{:.2e}'.format(i) for i in ev[l]])
    print('Fraction of compressive velocity: ', ['{:.2e}'.format(i) for i in frac_vcomp[l]])
    print('Fraction of rotational velocity: ', ['{:.2e}'.format(i) for i in frac_vrot[l]], '\n')


os.chdir(initialpath)

if verbose:
    print('#####################')
    print('Done with test 1b!')
    print('#####################')

#######################################
# test1c: constant divergence field in 'realistic' grid structure
#######################################
if verbose:
    print('#####################')
    print('Starting with test 1c')
    print('#####################')

output_path = os.path.join(initialpath, 'test_files/1c')
maxl = 9
nmax = 128
size = 1

expected_div = 0.03
expected_vmax = 0.01

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
                                   grids_path=os.path.join(output_path, 'simu_masclet'), parameters_path=output_path,
                                   digits=5, are_divrot=True, are_potentials=True, are_velocities=True, verbose=False)

# test divergence
ediv = {}
erot = {}
for l in range(maxl+1):
    Dx = size/nmax/2**l

    if l == 0:
        minipatch = 0
        maxipatch = 0
    else:
        minipatch = 8 * (l-1) + 1
        maxipatch = 8 * l

    divs = np.concatenate([a.flatten() for a in div[minipatch:maxipatch+1]])
    rotsx = np.concatenate([a.flatten() for a in rotx[minipatch:maxipatch+1]])
    rotsy = np.concatenate([a.flatten() for a in roty[minipatch:maxipatch+1]])
    rotsz = np.concatenate([a.flatten() for a in rotz[minipatch:maxipatch+1]])
    rots = np.sqrt(rotsx ** 2 + rotsy ** 2 + rotsz ** 2)

    ediv[l] = [np.percentile(np.abs(divs/expected_div - 1), perc) for perc in percentiles]
    erot[l] = [np.percentile(rots, perc) * Dx / expected_vmax for perc in percentiles]

    print('* At level {}, percentiles: '.format(l), percentiles)
    print('Relative errors in divergence: ', ['{:.2e}'.format(i) for i in ediv[l]])
    print('Errors in rotational: ', ['{:.2e}'.format(i) for i in erot[l]], '\n')

# test velocities
# first assert the written total velocities are the original ones
for ipatch in range(len(vx)):
    if ((vx[ipatch]-truevx[ipatch])**2 + (vy[ipatch]-truevy[ipatch])**2 + (vz[ipatch]-truevz[ipatch])**2).sum() > 0:
        print('Problem with the original velocity in patch {}!'.format(ipatch))

ev = {}
frac_vcomp = {}
frac_vrot = {}
for l in range(maxl+1):
    if l == 0:
        minipatch = 0
        maxipatch = 0
    else:
        minipatch = npatch[0:l].sum()
        maxipatch = npatch[0:l+1].sum()

    vx_l = np.concatenate([a.flatten() for a in vx[minipatch:maxipatch+1]])
    vy_l = np.concatenate([a.flatten() for a in vy[minipatch:maxipatch+1]])
    vz_l = np.concatenate([a.flatten() for a in vz[minipatch:maxipatch+1]])
    v_l = np.sqrt(vx_l ** 2 + vy_l ** 2 + vz_l ** 2)

    velcompx_l = np.concatenate([a.flatten() for a in velcompx[minipatch:maxipatch+1]])
    velcompy_l = np.concatenate([a.flatten() for a in velcompy[minipatch:maxipatch+1]])
    velcompz_l = np.concatenate([a.flatten() for a in velcompz[minipatch:maxipatch+1]])
    velcomp_l = np.sqrt(velcompx_l ** 2 + velcompy_l ** 2 + velcompz_l ** 2)

    velrotx_l = np.concatenate([a.flatten() for a in velrotx[minipatch:maxipatch+1]])
    velroty_l = np.concatenate([a.flatten() for a in velroty[minipatch:maxipatch+1]])
    velrotz_l = np.concatenate([a.flatten() for a in velrotz[minipatch:maxipatch+1]])
    velrot_l = np.sqrt(velrotx_l ** 2 + velroty_l ** 2 + velrotz_l ** 2)

    velrecx_l = velcompx_l + velrotx_l
    velrecy_l = velcompy_l + velroty_l
    velrecz_l = velcompz_l + velrotz_l

    ev[l] = [np.mean([np.percentile(np.abs(velrecx_l/vx_l - 1), perc),
                     np.percentile(np.abs(velrecy_l/vy_l - 1), perc),
                     np.percentile(np.abs(velrecz_l/vz_l - 1), perc)]) for perc in percentiles]
    frac_vcomp[l] = [np.percentile(velcomp_l / v_l, perc) for perc in percentiles]
    frac_vrot[l] = [np.percentile(velrot_l / v_l, perc) for perc in percentiles]

    print('* At level {}, percentiles: '.format(l), percentiles)
    print('Relative errors in velocity reconstruction: ', ['{:.2e}'.format(i) for i in ev[l]])
    print('Fraction of compressive velocity: ', ['{:.2e}'.format(i) for i in frac_vcomp[l]])
    print('Fraction of rotational velocity: ', ['{:.2e}'.format(i) for i in frac_vrot[l]], '\n')


os.chdir(initialpath)

if verbose:
    print('#####################')
    print('Done with test 1c!')
    print('#####################')

#### END MESSAGE
if verbose:
    print('Done!')
