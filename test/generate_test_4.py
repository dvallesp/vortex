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
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator

mockit = 4000
ncores = 8
verbose = True
smooth = True
if smooth:
    sigma = 2
write = True

#######################################
# test1d: multifrequency field
#######################################
if verbose:
    print('Starting with test 1d')

output_path = os.path.join(initialpath, 'test_files/1d/simu_masclet')
maxl = 4
nmax = 128
size = 1
uni_nmax = nmax * 2**maxl

slope_comp = -2
slope_sol = -5/3
min_n = 2
max_n = 1024
norm_velocities = 0.01
nnorm_velocities = 64

if verbose:
    print('Grids')

npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, \
patchry, patchrz, pare = mock_gridsAprime(nlevels=maxl, Nx=nmax, L=size)

levels = masclet.tools.create_vector_levels(npatch)

if write:
    write_mock_grids(mockit, 100, maxl, 1e-10, 0.1, 100000, npatch, patchnx, patchny, patchnz, patchx, patchy, patchz,
                     patchrx, patchry, patchrz, pare, digits=5, path=output_path)

if verbose:
    print('Velocities')

cellsrx, cellsry, cellsrz = masclet.tools_xyz.compute_position_fields(patchnx, patchny, patchnz, patchrx, patchry,
                                                                      patchrz, npatch, size, nmax, ncores=ncores)

### GENERATE THE VELOCITY FIELDS
uni_cellsize = size / uni_nmax
uni_cellsrx = np.linspace(-size/2 + uni_cellsize/2, size/2 - uni_cellsize/2, uni_nmax)
uni_cellsry = np.linspace(-size/2 + uni_cellsize/2, size/2 - uni_cellsize/2, uni_nmax)
uni_cellsrz = np.linspace(-size/2 + uni_cellsize/2, size/2 - uni_cellsize/2, uni_nmax)
uni_cells = {0: uni_cellsrx, 1: uni_cellsry, 2: uni_cellsrz}

randomphases = {i: [] for i in range(1,10)}
for i in range(1,10):
    randomphases[i] = np.random.rand(max_n-min_n+1) * 2 * np.pi

for mode in ['sol', 'comp']:
    for icomp in range(3):

        print('#################')
        print('Component {:}'.format(icomp), mode)
        print('#################')

        uniform_v = np.zeros((uni_nmax, uni_nmax, uni_nmax))

        addx = np.zeros(uni_cells[0].shape)
        addy = np.zeros(uni_cells[0].shape)
        addz = np.zeros(uni_cells[0].shape)

        for n in range(min_n, max_n+1):
            #print('Term {:}'.format(n))
            Acomp = norm_velocities * (n/nnorm_velocities) ** ((slope_comp+1) / 2)
            Asol = norm_velocities * (n/nnorm_velocities) ** ((slope_sol+1) / 2)
            kn = 2*np.pi*n/size

            for jcomp in range(3):

                if (icomp == jcomp and mode == 'sol') or (icomp != jcomp and mode == 'comp'):
                    continue
                iphase = icomp * 3 + jcomp + 1
                if mode == 'comp':
                    add = Acomp * np.sin(kn * uni_cells[jcomp] + randomphases[iphase][n-min_n])
                else:
                    add = Asol * np.sin(kn * uni_cells[jcomp] + randomphases[iphase][n - min_n])

                if jcomp==0:
                    addx += add
                elif jcomp==1:
                    addy += add
                elif jcomp==2:
                    addz += add

        for ii in range(uniform_v.shape[0]):
            #print('x', ii)
            uniform_v[ii,:,:] += addx[ii]
        for jj in range(uniform_v.shape[1]):
            #print('y', jj)
            uniform_v[:,jj,:] += addy[jj]
        for kk in range(uniform_v.shape[2]):
            #print('z', kk)
            uniform_v[:,:,kk] += addz[kk]

        print('Uniform done.')
        if mode == 'comp':
            if icomp == 0:
                velcompx = []
            elif icomp == 1:
                velcompy = []
            elif icomp == 2:
                velcompz = []
        elif mode == 'sol':
            if icomp == 0:
                velrotx = []
            elif icomp == 1:
                velroty = []
            elif icomp == 2:
                velrotz = []

        for l in range(maxl+1):
            reduction = 2**(maxl-l)
            if l==0:
                minipatch = 0
                maxipatch = 0
            else:
                minipatch = npatch[0:l].sum()+1
                maxipatch = npatch[0:l+1].sum()
            for ipatch in range(minipatch, maxipatch+1):
                print('Patch {:}'.format(ipatch))
                vpatch = np.zeros((patchnx[ipatch], patchny[ipatch], patchnz[ipatch]))

                st_x = 0
                st_y = 0
                st_z = 0
                ip = ipatch
                while ip != 0:
                    reduction_ip = 2**(maxl - levels[ip])
                    st_x += (patchx[ip]-1) * reduction_ip * 2
                    st_y += (patchy[ip]-1) * reduction_ip * 2
                    st_z += (patchz[ip]-1) * reduction_ip * 2
                    ip = pare[ip]
                #print(ipatch, st_x, st_y, st_z)
                if not smooth:
                    for i in range(patchnx[ipatch]):
                        for j in range(patchny[ipatch]):
                            for k in range(patchnz[ipatch]):
                                vpatch[i,j,k] = uniform_v[st_x+i*reduction:st_x+(i+1)*reduction,
                                                          st_y+j*reduction:st_y+(j+1)*reduction,
                                                          st_z+k*reduction:st_z+(k+1)*reduction].mean()
                else:
                    if ipatch == 0:
                        for i in range(patchnx[ipatch]):
                            for j in range(patchny[ipatch]):
                                for k in range(patchnz[ipatch]):
                                    vpatch[i, j, k] = uniform_v[st_x + i * reduction:st_x + (i + 1) * reduction,
                                                      st_y + j * reduction:st_y + (j + 1) * reduction,
                                                      st_z + k * reduction:st_z + (k + 1) * reduction].mean()
                        vpatch = gaussian_filter(vpatch, sigma=sigma, mode='wrap')
                    else:
                        vpatch_ext = np.zeros((patchnx[ipatch]+8*sigma, patchny[ipatch]+8*sigma, patchnz[ipatch]+8*sigma))
                        st_x = st_x - 4 * reduction * sigma
                        st_y = st_y - 4 * reduction * sigma
                        st_z = st_z - 4 * reduction * sigma
                        for i in range(vpatch_ext.shape[0]):
                            for j in range(vpatch_ext.shape[1]):
                                for k in range(vpatch_ext.shape[2]):
                                    vpatch_ext[i,j,k] = \
                                        uniform_v[st_x + i * reduction:st_x + (i + 1) * reduction,
                                                  st_y + j * reduction:st_y + (j + 1) * reduction,
                                                  st_z + k * reduction:st_z + (k + 1) * reduction].mean()
                        vpatch = gaussian_filter(vpatch_ext, sigma=sigma, mode='nearest')[4*sigma:patchnx[ipatch]+4*sigma,
                                                                                      4*sigma:patchny[ipatch]+4*sigma,
                                                                                      4*sigma:patchnz[ipatch]+4*sigma]

                if mode == 'comp':
                    if icomp == 0:
                        velcompx.append(vpatch)
                    elif icomp == 1:
                        velcompy.append(vpatch)
                    elif icomp == 2:
                        velcompz.append(vpatch)
                elif mode == 'sol':
                    if icomp == 0:
                        velrotx.append(vpatch)
                    elif icomp == 1:
                        velroty.append(vpatch)
                    elif icomp == 2:
                        velrotz.append(vpatch)
    if write:
        if mode == 'comp':
            write_mock_clus(velcompx, velcompy, velcompz, mockit, digits=5, path=os.path.join(output_path,'comp'))
        elif mode == 'sol':
            write_mock_clus(velrotx, velroty, velrotz, mockit, digits=5, path=os.path.join(output_path, 'sol'))

vx = [c+s for c,s in zip(velcompx, velrotx)]
vy = [c+s for c,s in zip(velcompy, velroty)]
vz = [c+s for c,s in zip(velcompz, velrotz)]

if write:
    write_mock_clus(vx, vy, vz, mockit, digits=5, path=output_path)

os.chdir(initialpath)

if verbose:
    print('Done!')


#### END MESSAGE
if verbose:
    print('Done!')
