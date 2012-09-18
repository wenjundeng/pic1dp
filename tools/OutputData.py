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
        self.nspecies, self.nmode, self.nx, self.nv = io.readInt(fdata, 4)
        self.mode = io.readInt(fdata, self.nmode)
        self.lx, self.v_max = io.readReal(fdata, 2)

        #print self.nspecies, self.nmode, self.nx, self.nv
        #print self.mode
        #print self.lx, self.v_max

        # generate arrays for axises
        self.x = np.arange(self.nx + 1.0) / self.nx * self.lx
        self.v = (np.arange(self.nv + 0.0) / (self.nv - 1) - 0.5) \
            * 2.0 * self.v_max
        self.xv = np.meshgrid(self.x, self.v)

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
                        rawdata.append(io.readScalar(fdata, self.nx * self.nv))

                    for i in range(3):
                        rawdata.append(io.readScalar(fdata, self.nv))

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
            ptcldist_xv = np.zeros((self.nv, self.nx + 1))
        else:
            ptcldist_xv = np.zeros((self.nv, self.nx))

        rawdata = self._rawdataset[itime]
        if ispecies < self.nspecies:
            ptcldist_xv[:, 0 : self.nx] \
                = rawdata[5 + ispecies * 6 + iptcldist].reshape((self.nv, self.nx))
        else:
            for i in range(self.nspecies):
                ptcldist_xv[:, 0 : self.nx] \
                    += rawdata[5 + i * 6 + iptcldist].reshape((self.nv, self.nx))
        if periodicbound:
            ptcldist_xv[:, self.nx] = ptcldist_xv[:, 0] # boundary condition

        return ptcldist_xv

    def get_ptcldist_v(self, itime, ispecies, iptcldist):
        '''get data of particle distribution in v space'''
        rawdata = self._rawdataset[itime]
        if ispecies < self.nspecies:
            return rawdata[8 + ispecies * 6 + iptcldist]
        else:
            ret = np.zeros((self.nv))
            for i in range(self.nspecies):
                ret += rawdata[8 + i * 6 + iptcldist]
            return ret

    def growthrate_energe_fit(self, time1, time2):
        '''calculate growth rate using data between time1 and time2'''
        scalar_t = self.get_scalar_t()
        itime1 = np.searchsorted(scalar_t[0], time1) - 1
        itime2 = np.searchsorted(scalar_t[0], time2)
        t = scalar_t[0, itime1 : itime2]
        energe = scalar_t[1, itime1 : itime2]

        # fit growthrate
        n = itime2 - itime1
        lnenerge = np.log(energe)
        sum_t = np.sum(t)
        sum_lnenerge = np.sum(lnenerge)
        sum_tlnenerge = np.sum(t * lnenerge)
        sum_t2 = np.sum(t * t)
        gamma = (n * sum_tlnenerge - sum_t * sum_lnenerge) \
            / (n * sum_t2 - sum_t * sum_t)
        return gamma

    def findpeak_energe(self, time1, time2):
        '''find peak using data between time1 and time2'''
        scalar_t = self.get_scalar_t()
        itime1 = np.searchsorted(scalar_t[0], time1) - 1
        itime2 = np.searchsorted(scalar_t[0], time2)
        t = scalar_t[0, itime1 : itime2]
        energe = scalar_t[1, itime1 : itime2]
        peakindex = np.argmax(energe)
        return [t[peakindex], energe[peakindex]]

