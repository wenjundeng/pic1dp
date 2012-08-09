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



    def

class VisualApp:
    """class of visualization app"""

    def __init__(self, datapath):
        """initialization"""

        self._datapath = datapath

        # layout
        self._fig = plt.figure(figsize = (22, 11))
        self._ax_xv = self._fig.add_axes([0.035, 0.53, 0.2, 0.44])
        self._ax_xv_colorbar = self._fig.add_axes([0.24, 0.53, 0.01, 0.44])
        self._ax_dist_x = self._fig.add_axes([0.035, 0.045, 0.2, 0.44])
        self._ax_dist_v = self._fig.add_axes([0.3, 0.53, 0.22, 0.44])



if __name__ == '__main__':
    parser = argparse.ArgumentParser( \
        description = 'Visualization app for PIC1D-PETSc')
    parser.add_argument('data_path', metavar = 'data path', \
        nargs = '?', type = str)
    args = parser.parse_args()

    app = VisualApp(args.data_path)
    plt.show()

