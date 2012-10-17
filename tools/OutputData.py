# Copyright 2012 Wenjun Deng <wdeng@wdeng.info>
#
# This file is part of PIC1D-PETSc
#
# PIC1D-PETSc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIC1D-PETSc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PIC1D-PETSc.  If not, see <http://www.gnu.org/licenses/>.


import os
import numpy as np
import XPetscBinaryIO

class OutputData:
    """class for handling PIC1D-PETSc output data"""

    def __init__(self, datapath):
        #self._datapath = datapath
        io = XPetscBinaryIO.XPetscBinaryIO()
        fdata = open(os.path.join(datapath, 'pic1dp.out'), 'rb')

        # read parameters
        self.nspecies, self.nmode, self.nx, self.nv, \
            self.nx_pd, self.nv_pd = io.readInt(fdata, 6)
        self.mode = io.readInt(fdata, self.nmode)
        self.lx, self.v_max = io.readReal(fdata, 2)

        #print self.nspecies, self.nmode, self.nx, self.nv
        #print self.mode
        #print self.lx, self.v_max

        # generate arrays for axises
        self.x = np.arange(self.nx + 1.0) / self.nx * self.lx
        self.x_pd = np.arange(self.nx_pd + 1.0) / self.nx_pd * self.lx
        self.v_pd = (np.arange(self.nv_pd + 0.0) / (self.nv_pd - 1) - 0.5) \
            * 2.0 * self.v_max
        self.xv_pd = np.meshgrid(self.x_pd, self.v_pd)

        # read time dependent data
        self._rawdataset = []
        try: # read until EOF
            while True:
                rawdata = []
                # scalars # index 0
                rawdata.append(io.readReal(fdata, self.nspecies * 3 + 2))
                # electric field Fourier components
                header = io.readInt(fdata, 1)[0]
                rawdata.append(io.readVec(fdata)) # index 1
                header = io.readInt(fdata, 1)[0]
                rawdata.append(io.readVec(fdata)) # index 2
                # electric field and charge density in x space
                header = io.readInt(fdata, 1)[0]
                rawdata.append(io.readVec(fdata)) # index 3
                header = io.readInt(fdata, 1)[0]
                rawdata.append(io.readVec(fdata)) # index 4
                # distribution in x-v
                for i in range(self.nspecies):
                    for i in range(3):
                        rawdata.append( \
                            io.readScalar(fdata, self.nx_pd * self.nv_pd))

                    for i in range(3):
                        rawdata.append(io.readScalar(fdata, self.nv_pd))

                #print rawdata
                self._rawdataset.append(rawdata)
        except (MemoryError, IndexError):
            pass

        fdata.close()

        self.ntime = len(self._rawdataset)
        print '# of time steps read:', self.ntime

    def get_scalar_t(self):
        '''get data of scalar vs t'''
        scalar_t = np.zeros(((self.nspecies + 1) * 3 + 2, self.ntime))
        for itime in range(self.ntime):
            rawdata = self._rawdataset[itime]
            scalar_t[0 : self.nspecies * 3 + 2, itime] = rawdata[0]
            # calculate summation from all species
            for ispecies in range(self.nspecies):
                scalar_t[self.nspecies * 3 + 2, itime] += rawdata[0][ispecies * 3 + 2]
                scalar_t[self.nspecies * 3 + 3, itime] += rawdata[0][ispecies * 3 + 3]
                scalar_t[self.nspecies * 3 + 4, itime] += rawdata[0][ispecies * 3 + 4]

        return scalar_t

    def get_mode_t(self):
        '''get data of electric field Fourier mode vs t'''
        mode_t = np.zeros((self.nmode * 2, self.ntime))
        for itime in range(self.ntime):
            rawdata = self._rawdataset[itime]
            mode_t[0 : self.nmode, itime] = rawdata[1]
            mode_t[self.nmode : self.nmode * 2, itime] = rawdata[2]

        return mode_t

    def get_field_x(self, itime):
        '''get data of electric field and charge density vs x'''
        field_x = np.zeros((2, self.nx + 1))
        rawdata = self._rawdataset[itime]
        field_x[0, : self.nx] = rawdata[3]
        field_x[1, : self.nx] = rawdata[4]

        # boundary condition
        for i in range(2):
            field_x[i, self.nx] = field_x[i, 0]

        return field_x

    def get_ptcldist_xv(self, itime, ispecies, iptcldist, periodicbound = True):
        '''get data of particle distribution in x-v plane'''
        if periodicbound:
            ptcldist_xv = np.zeros((self.nv_pd, self.nx_pd + 1))
        else:
            ptcldist_xv = np.zeros((self.nv_pd, self.nx_pd))

        rawdata = self._rawdataset[itime]
        if ispecies < self.nspecies:
            ptcldist_xv[:, 0 : self.nx_pd] \
                = rawdata[5 + ispecies * 6 + iptcldist].reshape( \
                (self.nv_pd, self.nx_pd))
        else:
            for i in range(self.nspecies):
                ptcldist_xv[:, 0 : self.nx_pd] \
                    += rawdata[5 + i * 6 + iptcldist].reshape((self.nv_pd, self.nx_pd))
        if periodicbound:
            ptcldist_xv[:, self.nx_pd] = ptcldist_xv[:, 0] # boundary condition

        return ptcldist_xv

    def get_ptcldist_v(self, itime, ispecies, iptcldist):
        '''get data of particle distribution in v space'''
        rawdata = self._rawdataset[itime]
        if ispecies < self.nspecies:
            return rawdata[8 + ispecies * 6 + iptcldist]
        else:
            ret = np.zeros((self.nv_pd))
            for i in range(self.nspecies):
                ret += rawdata[8 + i * 6 + iptcldist]
            return ret

    def growthrate_energy_fit(self, time1, time2):
        '''calculate growth rate using data between time1 and time2'''
        scalar_t = self.get_scalar_t()
        itime1 = np.searchsorted(scalar_t[0], time1) - 1
        itime2 = np.searchsorted(scalar_t[0], time2)
        t = scalar_t[0, itime1 : itime2]
        energy = scalar_t[1, itime1 : itime2]

        # fit growthrate
        n = itime2 - itime1
        lnenergy = np.log(energy)
        sum_t = np.sum(t)
        sum_lnenergy = np.sum(lnenergy)
        sum_tlnenergy = np.sum(t * lnenergy)
        sum_t2 = np.sum(t * t)
        gamma = (n * sum_tlnenergy - sum_t * sum_lnenergy) \
            / (n * sum_t2 - sum_t * sum_t)
        return gamma

    def findpeak_energy(self, time1, time2):
        '''find peak using data between time1 and time2'''
        scalar_t = self.get_scalar_t()
        itime1 = np.searchsorted(scalar_t[0], time1) - 1
        itime2 = np.searchsorted(scalar_t[0], time2)
        t = scalar_t[0, itime1 : itime2]
        energy = scalar_t[1, itime1 : itime2]
        peakindex = np.argmax(energy)
        return [t[peakindex], energy[peakindex]]

