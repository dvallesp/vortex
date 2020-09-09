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
    vx = [np.nan_to_num(vi * y / rho) for vi, y, rho in zip(v, cellsry, cellsrho)]
    vy = [np.nan_to_num(- vi * x / rho) for vi, x, rho in zip(v, cellsrx, cellsrho)]
    vz = [np.zeros(vxi.shape) for vxi in vx]
    return vx, vy, vz


def constant_div_field(cellsrx, cellsry, cellsrz, constant):
    cellsr = [np.sqrt(cellrx ** 2 + cellry ** 2 + cellrz ** 2) for cellrx, cellry, cellrz in zip(cellsrx, cellsry,
                                                                                                 cellsrz)]
    v = [constant * r for r in cellsr]
    vx = [np.nan_to_num(vi * x / r) for vi, x, r in zip(v, cellsrx, cellsr)]
    vy = [np.nan_to_num(vi * y / r) for vi, y, r in zip(v, cellsry, cellsr)]
    vz = [np.nan_to_num(vi * z / r) for vi, z, r in zip(v, cellsrz, cellsr)]
    return vx, vy, vz


def periodic(cellsrx, cellsry, cellsrz, constant):
    vx = [constant * np.sin(2 * np.pi * x) for x in cellsrx]
    vy = [np.zeros(x.size) for x in cellsrx]
    vz = [np.zeros(x.size) for x in cellsrx]
    return vx, vy, vz


def periodic_compressive(cellsrx, cellsry, cellsrz, constant):
    vx = [constant * np.sin(2 * np.pi * x) for x, y, z in zip(cellsrx, cellsry, cellsrz)]
    vy = [constant * np.sin(2 * np.pi * y) for x, y, z in zip(cellsrx, cellsry, cellsrz)]
    vz = [constant * np.sin(2 * np.pi * z) for x, y, z in zip(cellsrx, cellsry, cellsrz)]
    return vx, vy, vz


def periodic_rotational(cellsrx, cellsry, cellsrz, constant):
    vx = [constant * (np.sin(4 * np.pi * y) + np.sin(6 * np.pi * z)) for x, y, z in zip(cellsrx, cellsry, cellsrz)]
    vy = [constant * (np.sin(4 * np.pi * z) + np.sin(6 * np.pi * x)) for x, y, z in zip(cellsrx, cellsry, cellsrz)]
    vz = [constant * (np.sin(4 * np.pi * x) + np.sin(6 * np.pi * y)) for x, y, z in zip(cellsrx, cellsry, cellsrz)]
    return vx, vy, vz
