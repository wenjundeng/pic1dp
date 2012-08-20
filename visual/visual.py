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

class XScalarFormatter(mpl.ticker.ScalarFormatter):
    """eXtended ScalarFormatter"""

    def __init__(self, useOffset = True, useMathText = False, \
        useLocale = None, precision = None):

        mpl.ticker.ScalarFormatter.__init__(\
            self, useOffset, useMathText, useLocale)
        self._precision = precision

    def set_precision(self, precision):
        """use this function to set precision"""
        self._precision = precision

    def get_offset(self):
        """Return scientific notation, plus offset"""
        if len(self.locs)==0: return ''
        s = ''
        if self.orderOfMagnitude or self.offset:
            offsetStr = ''
            sciNotStr = ''
            if self.offset:
                offsetStr = self.format_data(self.offset)
                if self.offset > 0: offsetStr = '+' + offsetStr
            if self.orderOfMagnitude:
                if self._usetex or self._useMathText:
                    sciNotStr = self.format_data(10**self.orderOfMagnitude)
                else:
                    sciNotStr = '1e%d'% self.orderOfMagnitude
            if self._usetex or self._useMathText:
                if sciNotStr != '':
                    sciNotStr = r'\times%s' % sciNotStr
                s =  ''.join(('$',sciNotStr,offsetStr,'$'))
            else:
                s =  ''.join((sciNotStr,offsetStr))

        return self.fix_minus(s)

    def _set_format(self):
        # set the format string to format all the ticklabels
        # The floating point black magic (adding 1e-15 and formatting
        # to 8 digits) may warrant review and cleanup.
        locs = (np.asarray(self.locs)-self.offset) / 10**self.orderOfMagnitude+1e-15
        sigfigs = [len(str('%1.8f'% loc).split('.')[1].rstrip('0')) \
                   for loc in locs]
        sigfigs.sort()
        if self._precision is None:
            self.format = '%1.' + str(sigfigs[-1]) + 'f'
        else:
            self.format = '%1.' + str(self._precision) + 'f'

        if self._usetex or self._useMathText:
            self.format = '$%s$' % self.format


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
        self.v = (np.arange(self.nv + 0.0) / (self.nv - 1) - 0.5) * 2.0 * self.v_max
        self.xv = np.meshgrid(self.x, self.v)

        # read time dependent data
        self._rawdataset = []
        try: # read until EOF
            while True:
                rawdata = []
                # scalars
                rawdata.append(io.readReal(fdata, self.nspecies * 3 + 2)) # index 0
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
        scalar_t = np.zeros((self.nspecies * 3 + 2, self.ntime))
        for itime in range(self.ntime):
            rawdata = self._rawdataset[itime]
            #print rawdata[0]
            #print 'aaa'
            scalar_t[:, itime] = rawdata[0]

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

    def get_ptcldist_xv(self, itime, ispecies, iptcldist):
        '''get data of particle distribution in x-v plane'''
        ptcldist_xv = np.zeros((self.nv, self.nx + 1))
        rawdata = self._rawdataset[itime]
        ptcldist_xv[:, 0 : self.nx] \
            = rawdata[5 + ispecies * 6 + iptcldist].reshape((self.nv, self.nx))
        ptcldist_xv[:, self.nx] = ptcldist_xv[:, 0] # boundary condition
        return ptcldist_xv

    def get_ptcldist_v(self, itime, ispecies, iptcldist):
        '''get data of particle distribution in v space'''
        rawdata = self._rawdataset[itime]
        return rawdata[8 + ispecies * 6 + iptcldist]


class VisualApp:
    """class of visualization app"""

    def __init__(self, datapath):
        """initialization"""

        # object that handles output data
        self._data = OutputData(datapath)

        # time dependent data for plotting
        self._scalar_t = self._data.get_scalar_t()
        self._mode_t = self._data.get_mode_t()

        # initial parameters
        self._iscalar = 0 # scalar index
        self._imode = 0 # mode index
        self._itime = 0 # time index
        self._itime1 = 0 # time range index 1
        self._itime2 = self._data.ntime # time range index 2
        self._ispecies = 0 # species index
        self._iptcldist = 0 # 0: g; 1: f; 2: delta f
        self._ani_playing = False # whether animation is playing

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

        # formatters
        self._scalar_t_formatter = XScalarFormatter( \
            useOffset = True, useMathText = True, precision = 1)
        self._scalar_t_formatter.set_powerlimits((-2, 3))
        self._ptcldist_xv_colorbar_formatter = XScalarFormatter( \
            useOffset = True, useMathText = True, precision = 1)
        self._ptcldist_xv_colorbar_formatter.set_powerlimits((-2, 3))

        # layout
        self._fig = plt.figure(figsize = (22, 11))
        # plots
        self._ax_scalar_t = self._fig.add_axes([0.04, 0.53, 0.18, 0.44])
        self._ax_mode_t = self._fig.add_axes([0.04, 0.045, 0.18, 0.44])
        self._ax_mode_amp_t = self._fig.add_axes([0.26, 0.53, 0.18, 0.44])
        self._ax_mode_norm_t = self._fig.add_axes([0.26, 0.045, 0.18, 0.44])
        self._ax_field_x = self._fig.add_axes([0.58, 0.53, 0.18, 0.44])
        #self._ax_chargeden_x = self._fig.add_axes([0.54, 0.045, 0.22, 0.44])
        self._ax_ptcldist_xv = self._fig.add_axes([0.8, 0.53, 0.16, 0.44])
        self._ax_ptcldist_xv_colorbar = self._fig.add_axes([0.965, 0.53, 0.01, 0.44])
        self._ax_ptcldist_v = self._fig.add_axes([0.8, 0.045, 0.18, 0.44])
        # widgets
        self._ax_scalar_chooser = self._fig.add_axes([0.45, 0.85, 0.08, 0.1])
        self._ax_mode_chooser = self._fig.add_axes([0.45, 0.7, 0.08, 0.1])
        self._ax_ptcldist_chooser = self._fig.add_axes([0.45, 0.55, 0.08, 0.1])
        self._ax_species_chooser = self._fig.add_axes([0.45, 0.4, 0.08, 0.1])
        self._ax_ani_playpause = self._fig.add_axes([0.45, 0.25, 0.08, 0.04])

        # clicking on window: time chooser and time range chooser
        self._fig.canvas.mpl_connect('button_press_event', self._canvas_on_press)
        self._fig.canvas.mpl_connect('button_release_event', self._canvas_on_release)

        # scalar chooser
        self._scalar_labels = ['$\int E^2 \mathrm{d} x$', \
            '$\int \, (g|f|\delta f) v^2 \mathrm{d} v \mathrm{d} x$']
        self._scalar_chooser = widgets.RadioButtons(self._ax_scalar_chooser, \
            self._scalar_labels, active = self._iscalar)
        self._scalar_chooser.on_clicked(self._scalar_chooser_on_click)

        # mode chooser
        self._mode_labels = []
        for i in range(self._data.nmode):
            self._mode_labels.append(str(self._data.mode[i]))
        self._mode_chooser = widgets.RadioButtons(self._ax_mode_chooser, \
            self._mode_labels, active = self._imode)
        self._mode_chooser.on_clicked(self._mode_chooser_on_click)

        # particle distribution chooser
        self._ptcldist_labels = ['$g$', '$f$', '$\delta f$']
        self._ptcldist_chooser = widgets.RadioButtons( \
            self._ax_ptcldist_chooser, self._ptcldist_labels, \
            active = self._iptcldist)
        self._ptcldist_chooser.on_clicked(self._ptcldist_chooser_on_click)

        # species chooser
        self._species_labels = []
        for i in range(self._data.nspecies):
            self._species_labels.append(str(i + 1))
        self._species_chooser = widgets.RadioButtons( \
            self._ax_species_chooser, self._species_labels, \
            active = self._ispecies)
        self._species_chooser.on_clicked(self._species_chooser_on_click)

        # animation play/pause button
        self._ani_playpause = widgets.Button(self._ax_ani_playpause, 'Play animation')
        self._ani_playpause.on_clicked(self._ani_playpause_on_click)

        # timer for animation
        self._timer = self._fig.canvas.new_timer(interval = 200)
        self._timer.add_callback(self._ani_advance)
        self._timer.start()

        # update all plots
        self.update_plot_all()

    def _canvas_on_press(self, event):
        # time range chooser
        if event.inaxes == self._ax_mode_t:
            self._itime1 = np.searchsorted(self._scalar_t[0], event.xdata) - 1
            if self._itime1 < 0:
                self._itime1 = 0
            if self._itime1 > self._data.ntime - 2:
                self._itime1 = self._data.ntime - 2
            t1 = self._scalar_t[0][self._itime1]
            self._ax_mode_t_tln1.set_xdata(np.array([t1, t1]))
            plt.draw()

    def _canvas_on_release(self, event):
        # time range chooser
        if event.inaxes == self._ax_mode_t:
            self._itime2 = np.searchsorted(self._scalar_t[0], event.xdata)
            if self._itime2 < self._itime1 + 2:
                self._itime2 = self._itime1 + 2
            if self._itime2 > self._data.ntime:
                self._itime2 = self._data.ntime
            #print self._itime2, self._itime1

            t2 = self._scalar_t[0][self._itime2 - 1]
            self._ax_mode_t_tln2.set_xdata(np.array([t2, t2]))
            self.update_plot_mode_ampnorm_t()
            plt.draw()

        # time chooser
        if (event.inaxes == self._ax_scalar_t) \
            or (event.inaxes == self._ax_mode_amp_t) \
            or (event.inaxes == self._ax_mode_norm_t):

            self._itime = np.searchsorted(self._scalar_t[0], event.xdata) - 1
            if self._itime < 0:
                self._itime = 0
            if self._itime > self._data.ntime - 1:
                self._itime = self._data.ntime - 1
            #print event.xdata, self._scalar_t[0], self._itime

            t = self._scalar_t[0][self._itime]
            self._ax_scalar_t_tln.set_xdata(np.array([t, t]))
            self._ax_mode_t_tln.set_xdata(np.array([t, t]))
            self._ax_mode_amp_t_tln.set_xdata(np.array([t, t]))
            self._ax_mode_norm_t_tln.set_xdata(np.array([t, t]))
            self.update_plot_field_x()
            self.update_plot_ptcldist_xv()
            self.update_plot_ptcldist_v()
            plt.draw()

    def _scalar_chooser_on_click(self, label):
        '''Deal with clicking scalar chooser'''
        #for i, l in zip(range(2), self._scalar_labels):
        self._iscalar = self._scalar_labels.index(label)
        self.update_plot_scalar_t()
        plt.draw()

    def _mode_chooser_on_click(self, label):
        '''Deal with clicking mode chooser'''
        self._imode = self._mode_labels.index(label)
        self.update_plot_mode_t()
        self.update_plot_mode_ampnorm_t()
        plt.draw()

    def _ptcldist_chooser_on_click(self, label):
        '''Deal with clicking particle distribution chooser'''
        self._iptcldist = self._ptcldist_labels.index(label)
        if self._iscalar > 0:
            self.update_plot_scalar_t()
        self.update_plot_ptcldist_xv()
        self.update_plot_ptcldist_v()
        plt.draw()

    def _species_chooser_on_click(self, label):
        '''Deal with clicking species chooser'''
        self._ispecies = self._species_labels.index(label)
        if self._iscalar > 0:
            self.update_plot_scalar_t()
        self.update_plot_ptcldist_xv()
        self.update_plot_ptcldist_v()
        plt.draw()

    def _ani_playpause_on_click(self, event):
        '''Deal with clicking animation play/pause button'''
        if self._ani_playing:
            self._ani_playing = False
            self._ani_playpause.label.set_text('Play animation')
            plt.draw()
            #self._timer.stop() # for some backends, timer can't stop
        else:
            self._ani_playing = True
            self._ani_playpause.label.set_text('Pause animation')
            #self._timer.start()

    def _ani_advance(self):
        '''Advance animation frame'''
        if not self._ani_playing:
            return
        self._itime = self._itime + 1
        if self._itime > self._data.ntime - 1:
            self._itime = self._data.ntime - 1
            self._ani_playing = False
            #self._ani_playpause_on_click(None)

        t = self._scalar_t[0][self._itime]
        self._ax_scalar_t_tln.set_xdata(np.array([t, t]))
        self._ax_mode_t_tln.set_xdata(np.array([t, t]))
        self._ax_mode_amp_t_tln.set_xdata(np.array([t, t]))
        self._ax_mode_norm_t_tln.set_xdata(np.array([t, t]))
        self.update_plot_field_x()
        self.update_plot_ptcldist_xv()
        self.update_plot_ptcldist_v()
        plt.draw()

    def update_plot_scalar_t(self):
        """Update plot for scalar vs t"""
        if self._iscalar == 0:
            iscalar = 1
        else:
            iscalar = 2 + self._ispecies * 3 + self._iptcldist

        self._ax_scalar_t.clear()
        self._ax_scalar_t.set_xlabel('$\omega_{\mathrm{pe}} t$')
        self._ax_scalar_t.set_ylabel(self._scalar_labels[self._iscalar])
        self._ax_scalar_t.yaxis.set_major_formatter(self._scalar_t_formatter)
        self._ax_scalar_t.plot( \
            self._scalar_t[0], self._scalar_t[iscalar])
        # time chooser indicator
        self._ax_scalar_t_tln = \
            self._ax_scalar_t.axvline(x = self._scalar_t[0][self._itime], \
            color = '0.5', ls = '--')

    def update_plot_mode_t(self):
        """Update plot for mode vs t"""
        self._ax_mode_t.clear()
        self._ax_mode_t.set_xlabel('$\omega_{\mathrm{pe}} t$')
        self._ax_mode_t.set_ylabel('$E$ Fourier mode, Re and Im')
        self._ax_mode_t.plot( \
            self._scalar_t[0], self._mode_t[self._imode], '-', \
            self._scalar_t[0], self._mode_t[self._data.nmode + self._imode], \
            '--')
        # time range chooser indicator
        self._ax_mode_t_tln1 = \
            self._ax_mode_t.axvline(x = self._scalar_t[0][self._itime1], \
            color = 'y', ls = '--')
        self._ax_mode_t_tln2 = \
            self._ax_mode_t.axvline(x = self._scalar_t[0][self._itime2 - 1], \
            color = 'y', ls = '--')
        # time chooser indicator
        self._ax_mode_t_tln = \
            self._ax_mode_t.axvline(x = self._scalar_t[0][self._itime], \
            color = '0.5', ls = '--')

    def update_plot_mode_ampnorm_t(self):
        """Update plot for mode amplitude and normalized vs t"""
        t = self._scalar_t[0][self._itime1 : self._itime2]
        re = self._mode_t[self._imode][self._itime1 : self._itime2]
        im = self._mode_t[self._data.nmode + self._imode] \
            [self._itime1 : self._itime2]
        amp = np.sqrt(re**2 + im**2)
        gamma = np.log(amp[self._itime2 - self._itime1 - 1] / amp[0]) \
            / (t[self._itime2 - self._itime1 - 1] - t[0])
        print 'gamma = ', gamma

        re_norm = re / np.exp(gamma * (t - t[0]))
        im_norm = im / np.exp(gamma * (t - t[0]))

        self._ax_mode_amp_t.clear()
        self._ax_mode_amp_t.set_xlabel('$\omega_{\mathrm{pe}} t$')
        self._ax_mode_amp_t.set_ylabel('Amplitude of $E$ Fourier mode')
        self._ax_mode_amp_t.plot(t, amp)
        self._ax_mode_amp_t.set_yscale('log')

        self._ax_mode_norm_t.clear()
        self._ax_mode_norm_t.set_xlabel('$\omega_{\mathrm{pe}} t$')
        self._ax_mode_norm_t.set_ylabel('Normalized $E$ Fourier mode')
        self._ax_mode_norm_t.plot( \
            t, re_norm, '-', \
            t, im_norm, '--')
        # time chooser indicator
        self._ax_mode_amp_t_tln = \
            self._ax_mode_amp_t.axvline(x = self._scalar_t[0][self._itime], \
            color = '0.5', ls = '--')
        self._ax_mode_norm_t_tln = \
            self._ax_mode_norm_t.axvline(x = self._scalar_t[0][self._itime], \
            color = '0.5', ls = '--')

    def update_plot_field_x(self):
        """Update plot for electric field and charge density vs x"""
        field_x = self._data.get_field_x(self._itime)

        self._ax_field_x.clear()
        self._ax_field_x.set_xlabel('$x$')
        self._ax_field_x.set_ylabel('$E$ and charge density')
        self._ax_field_x.plot( \
            self._data.x, field_x[0], '-', \
            self._data.x, field_x[1], '--')

    def update_plot_ptcldist_xv(self):
        """Update plot for particle distribution in x-v plane"""
        ptcldist_xv = self._data.get_ptcldist_xv( \
            self._itime, self._ispecies, self._iptcldist)

        self._ax_ptcldist_xv.clear()
        self._ax_ptcldist_xv_colorbar.clear()
        self._ax_ptcldist_xv.set_xlabel('$x$')
        self._ax_ptcldist_xv.set_ylabel('$v$')
        cf = self._ax_ptcldist_xv.contourf( \
            self._data.xv[0], self._data.xv[1], ptcldist_xv)
        plt.colorbar(cf, cax = self._ax_ptcldist_xv_colorbar, \
            format = self._ptcldist_xv_colorbar_formatter)

    def update_plot_ptcldist_v(self):
        '''Update plot for particle distribution in v space'''
        ptcldist_v = self._data.get_ptcldist_v( \
            self._itime, self._ispecies, self._iptcldist)

        self._ax_ptcldist_v.clear()
        self._ax_ptcldist_v.set_xlabel('$v$')
        self._ax_ptcldist_v.set_ylabel('particle distribution')
        self._ax_ptcldist_v.plot(self._data.v, ptcldist_v)

    def update_plot_all(self):
        """Update all plots"""
        self.update_plot_scalar_t()
        self.update_plot_mode_t()
        self.update_plot_mode_ampnorm_t()
        self.update_plot_field_x()
        self.update_plot_ptcldist_xv()
        self.update_plot_ptcldist_v()
        plt.draw()


if __name__ == '__main__':
    parser = argparse.ArgumentParser( \
        description = 'Visualization app for PIC1D-PETSc')
    parser.add_argument('data_path', metavar = 'data path', \
        nargs = '?', type = str, default = './')
    args = parser.parse_args()

    app = VisualApp(args.data_path)
    plt.show()

