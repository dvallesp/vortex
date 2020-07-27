################################################################################
#                          GENERATE TESTS FOR VORTEX                           #
#                                                                              #
#                          Masclet mock files module                           #
#                                                                              #
#                            Author: David Vallés                              #
################################################################################


def write_mock_grids(it, t, nl, mdm, zeta, ndxyz, npatch, patchnx, patchny, patchnz, patchx, patchy, patchz, patchrx,
                     patchry, patchrz, pare, digits=5):
    with open(os.path.join(path, masclet.read_masclet.filename(it, 'g', digits)), 'w') as grids:
        grids.write('{}\t{}\t{}\t{}\t{}'.format(it, t, nl, mdm, 0))
        grid.write('{}'.format(zeta))
        grids.write('{}\t{}\t{}'.format(0, ndxyz, 0))

        for l in range(1, len(npatch)):
            grids.write('{}\t{}'.format(l, npatch[l]))
            grids.write(' ----------------- within level=           {}  -----------'.format(l))

            for ipatch in range(sum(npatch[0:ir]) + 1, sum(npatch[0:ir + 1]) + 1):
                grids.write('{}\t{}\t{}'.format(patchnx[ipatch], patchny[ipatch], patchnz[ipatch]))
                grids.write('{}\t{}\t{}'.format(patchx[ipatch], patchy[ipatch], patchz[ipatch]))
                grids.write('{}\t{}\t{}'.format(patchrx[ipatch], patchry[ipatch], patchrz[ipatch]))
                grid.write('{}'.format(pare[ipatch]))


def write_mock_clus(vx, vy, vz, it, digits=5):
    with FF(os.path.join(path, masclet.read_masclet.filename(it, 'g', digits)), 'w') as clus:
        clus.write_record(np.array([0]).astype('i4')) #first line
        # base level
        clus.write_record(np.array([0]).astype('i4'))
        clus.write_record(vx[0].T.astype('f4'))
        clus.write_record(vy[0].T.astype('f4'))
        clus.write_record(vz[0].T.astype('f4'))
        for i in range(6):
            clus.write_record(np.array([0]).astype('i4'))

        # refinement levels
        for vxi, vyi, vzi in zip(vx[1:], vy[1:], vz[1:]):
            clus.write_record(np.array([0]).astype('i4'))
            clus.write_record(vxi.T.astype('f4'))
            clus.write_record(vyi.T.astype('f4'))
            clus.write_record(vzi.T.astype('f4'))
            for i in range(7):
                clus.write_record(np.array([0]).astype('i4'))
