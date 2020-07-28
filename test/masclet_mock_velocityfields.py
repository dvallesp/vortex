################################################################################
#                          GENERATE TESTS FOR VORTEX                           #
#                                                                              #
#                                                                              #
#                           Author: David Vallés                               #
################################################################################

import numpy as np


def constant_curl_field(cellsrx, cellsry, cellsrz, constant):
    cellsrho = [np.sqrt(cellrx ** 2 + cellry ** 2) for cellrx, cellry, cellrz in zip(cellsrx, cellsry, cellsrz)]
    v = [constant * rho for rho in cellsrho]
    vx = [vi * y / rho for vi, y, rho in zip(v, cellsry, cellsrho)]
    vy = [- vi * x / rho for vi, x, rho in zip(v, cellsrx, cellsrho)]
    vz = [np.zeros(vxi.shape) for vxi in vx]
    return vx, vy, vz


def constant_div_field(cellsrx, cellsry, cellsrz, constant):
    cellsr = [np.sqrt(cellrx ** 2 + cellry ** 2 + cellrz ** 2) for cellrx, cellry, cellrz in zip(cellsrx, cellsry,
                                                                                                 cellsrz)]
    v = [constant * r for r in cellsr]
    vx = [vi * x / r for vi, x, r in zip(v, cellsrx, cellsr)]
    vy = [vi * y / r for vi, y, r in zip(v, cellsry, cellsr)]
    vz = [vi * z / r for vi, z, r in zip(v, cellsrz, cellsr)]
    return vx, vy, vz
