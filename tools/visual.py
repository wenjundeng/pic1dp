#!/usr/bin/env python

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


"""PIC1D-PETSc visualization app"""

import argparse
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
import matplotlib.widgets as widgets
import OutputData
import XScalarFormatter

class VisualApp:
    """Class of visualization app"""

    def __init__(self, datapath):
        """Initialization"""

        # object that handles output data
        self._data = OutputData.OutputData(datapath)

        # time dependent data for plotting
        self._scalar_t = self._data.get_scalar_t()
        self._mode_t = self._data.get_mode_t()

        # initial parameters
        self._iscalar = 0 # scalar index
        self._imode = 0 # mode index
        self._itime = 0 # time index
        self._itime1 = 0 # time range index 1
        self._itime2 = self._data.ntime # time range index 2
        self._ispecies = self._data.nspecies # species index
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
        self._scalar_t_formatter = XScalarFormatter.XScalarFormatter( \
            useOffset = True, useMathText = True, precision = 2)
        self._scalar_t_formatter.set_powerlimits((-2, 3))
        self._ptcldist_xv_colorbar_formatter \
            = XScalarFormatter.XScalarFormatter( \
            useOffset = True, useMathText = True, precision = 2)
        self._ptcldist_xv_colorbar_formatter.set_powerlimits((-2, 3))

        # layout
        self._fig = plt.figure(figsize = (22, 11))
        self._fig.canvas.set_window_title(os.path.abspath(datapath))
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
        self._ax_scalar_chooser.set_title('Scalar')
        self._ax_mode_chooser = self._fig.add_axes([0.45, 0.7, 0.08, 0.1])
        self._ax_mode_chooser.set_title('Fourier mode')
        self._ax_ptcldist_chooser = self._fig.add_axes([0.45, 0.55, 0.08, 0.1])
        self._ax_ptcldist_chooser.set_title('Distribution')
        self._ax_species_chooser = self._fig.add_axes([0.45, 0.4, 0.08, 0.1])
        self._ax_species_chooser.set_title('Species')
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
        for ispecies in range(self._data.nspecies):
            self._species_labels.append(str(ispecies + 1))
        self._species_labels.append('Sum')
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

        # not playing, return directly
        if not self._ani_playing:
            return

        # temporarily turn off playing status to avoid double advancing
        self._ani_playing = False

        self._itime = self._itime + 1
        if self._itime > self._data.ntime - 1:
            self._itime = self._data.ntime - 1

        t = self._scalar_t[0][self._itime]
        self._ax_scalar_t_tln.set_xdata(np.array([t, t]))
        self._ax_mode_t_tln.set_xdata(np.array([t, t]))
        self._ax_mode_amp_t_tln.set_xdata(np.array([t, t]))
        self._ax_mode_norm_t_tln.set_xdata(np.array([t, t]))
        self.update_plot_field_x()
        self.update_plot_ptcldist_xv()
        self.update_plot_ptcldist_v()

        if self._itime >= self._data.ntime - 1:
            self._ani_playing = False
            self._ani_playpause.label.set_text('Play animation')
            plt.draw()
        else:
            # keep on playing
            self._ani_playing = True
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
        nlevel = 64
        fmax = np.max(ptcldist_xv)
        fmin = np.min(ptcldist_xv)
        levels = fmin + (fmax - fmin) * np.arange(nlevel) / (nlevel - 1.0)
        cf = self._ax_ptcldist_xv.contourf( \
            self._data.xv[0], self._data.xv[1], ptcldist_xv, levels)
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
# end of class VisualApp


if __name__ == '__main__':
    parser = argparse.ArgumentParser( \
        description = 'Visualization app for PIC1D-PETSc')
    parser.add_argument('data_path', metavar = 'data path', \
        nargs = '?', type = str, default = './')
    args = parser.parse_args()

    visapp = VisualApp(args.data_path)
    plt.show()

