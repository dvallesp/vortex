################################################################################
#                          GENERATE TESTS FOR VORTEX                           #
#                                                                              #
#                          Masclet mock grids module                           #
#                                                                              #
#                            Author: David Vallés                              #
################################################################################

import numpy as np


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
                    
                    patchx[ipatch] = (patchnx[ipatch] / 2 * (i + 1)) % patchnx[ipare] + 1
                    patchy[ipatch] = (patchny[ipatch] / 2 * (j + 1)) % patchny[ipare] + 1
                    patchz[ipatch] = (patchnz[ipatch] / 2 * (k + 1)) % patchnz[ipare] + 1

                    patchrx[ipatch] = patchrx[ipare] + (patchx[ipatch] - 3 / 2) * L / Nx / 2 ** (l - 1)
                    patchry[ipatch] = patchry[ipare] + (patchy[ipatch] - 3 / 2) * L / Nx / 2 ** (l - 1)
                    patchrz[ipatch] = patchrz[ipare] + (patchz[ipatch] - 3 / 2) * L / Nx / 2 ** (l - 1)

    return npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx, patchry, patchrz, pare