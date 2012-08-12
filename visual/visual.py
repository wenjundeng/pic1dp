#!/usr/bin/env python
"""PIC1D-PETSc visualization app"""

import os
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import matplotlib.widgets as widgets
import XPetscBinaryIO

class OutputData:
    """class for handling PIC1D-PETSc output data"""

    def __init__(self, datapath):
        #self._datapath = datapath
        io = XPetscBinaryIO.XPetscBinaryIO()
        fdata = open(os.path.join(datapath, 'pic1dp.out'), 'rb')

        self.lx, self.v_max = io.readReal(fdata, 2)
        self.nx, self.nv = io.readInt(fdata, 2)

        ptcldist = io.readBinaryFile(fdata)[0]
        #print ptcldist
        ptcldist.resize(self.nv, self.nx)

        self.ptcldist_xv = np.zeros((self.nv, self.nx + 1))
        self.ptcldist_xv[:, 0 : self.nx] = ptcldist
        self.ptcldist_xv[:, self.nx] = ptcldist[:, 0]
        #print self.ptcldist_xv

        self.x = np.arange(self.nx + 1.0) / self.nx * self.lx
        self.v = (np.arange(self.nv + 0.0) / (self.nv - 1) - 0.5) * 2.0 * self.v_max
        self.xv = np.meshgrid(self.x, self.v)
        print self.x
        print self.v
        print self.xv

    def getptcldist():
        pass


class VisualApp:
    """class of visualization app"""

    def __init__(self, datapath):
        """initialization"""

        # object that handles output data
        self._data = OutputData(datapath)

        # colormap
        cdict = {'red':   [(0.0,  0.0, 0.0),
                           (0.5,  1.0, 1.0),
                           (1.0,  1.0, 1.0)],

                 'green': [(0.0,  0.0, 0.0),
                           (0.5,  1.0, 1.0),
                           (1.0,  0.0, 0.0)],

                 'blue':  [(0.0,  1.0, 1.0),
                           (0.5,  1.0, 1.0),
                           (1.0,  0.0, 0.0)]}

        self._cmap = mpl.colors.LinearSegmentedColormap('BWR', cdict, 256)
        self._levels = (np.arange(64) - 31.5) / 31.5

        # layout
        self._fig = plt.figure(figsize = (22, 11))
        self._ax_xv = self._fig.add_axes([0.035, 0.53, 0.2, 0.44])
        self._ax_xv_colorbar = self._fig.add_axes([0.24, 0.53, 0.01, 0.44])
        self._ax_dist_x = self._fig.add_axes([0.035, 0.045, 0.2, 0.44])
        self._ax_dist_v = self._fig.add_axes([0.3, 0.53, 0.22, 0.44])

        #print self._data.xv[0].shape
        #print self._data.xv[1].shape
        #print self._data.ptcldist_xv.shape
        self._ax_xv.contourf(self._data.xv[0], self._data.xv[1], self._data.ptcldist_xv)



if __name__ == '__main__':
    parser = argparse.ArgumentParser( \
        description = 'Visualization app for PIC1D-PETSc')
    parser.add_argument('data_path', metavar = 'data path', \
        nargs = '?', type = str, default = './')
    args = parser.parse_args()

    app = VisualApp(args.data_path)
    plt.show()

