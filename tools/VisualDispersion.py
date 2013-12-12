# Copyright 2012, 2013 Wenjun Deng <wdeng@wdeng.info>
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


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import XScalarFormatter

class VisualDispersion:
    '''Class of visualization of dispersion'''

    def __init__(self, disp, arrk, arromega):
        self._disp = disp
        self._arrk = arrk
        self._arromega = arromega
        self._ispecies = disp.nspecies
        # formatters
        self._modestruct_colorbar_formatter \
            = XScalarFormatter.XScalarFormatter( \
                useOffset = True, useMathText = True, precision = 2)
        self._modestruct_colorbar_formatter.set_powerlimits((-2, 3))
        # layout
        self._fig = plt.figure(figsize = (13, 5))
        self._fig.canvas.set_window_title('Visual Dispersion')
        self._ax_omega_re = self._fig.add_axes([0.07, 0.1, 0.3, 0.82])
        self._ax_omega_im = self._ax_omega_re.twinx()
        self._ax_modestruct = self._fig.add_axes([0.6, 0.1, 0.3, 0.82])
        self._ax_modestruct_colorbar \
            = self._fig.add_axes([0.92, 0.1, 0.02, 0.82])
        # widgets
        self._ax_species_chooser = self._fig.add_axes([0.43, 0.62, 0.1, 0.3])
        self._ax_species_chooser.set_title('Species')

        # clicking on window, k chooser
        self._fig.canvas.mpl_connect('button_release_event', self._canvas_on_release)

        # species chooser
        self._species_labels = []
        for ispecies in range(self._disp.nspecies):
            self._species_labels.append(str(ispecies + 1))
        self._species_labels.append('Sum')
        self._species_chooser = widgets.RadioButtons( \
            self._ax_species_chooser, self._species_labels, \
            active = self._ispecies)
        self._species_chooser.on_clicked(self._species_chooser_on_click)

        self.update_plot_omega()
        self.update_plot_modestruct()
        plt.draw()

    def _canvas_on_release(self, event):
        # k chooser
        if event.inaxes == self._ax_omega_im:
            self._disp.set_k(event.xdata)
            ik = np.searchsorted(self._arrk, event.xdata)
            self._disp.append_guess(self._arromega[ik - 2 : ik + 2])
            self._disp.solveomega()
            self._disp.print_komega()

            self.update_plot_omega()
            self.update_plot_modestruct()
            plt.draw()

    def _species_chooser_on_click(self, label):
        '''Deal with clicking species chooser'''
        self._ispecies = self._species_labels.index(label)
        self.update_plot_modestruct()
        plt.draw()

    def update_plot_omega(self):
        '''Update omega plot'''
        self._ax_omega_re.clear()
        self._ax_omega_im.clear()
        self._ax_omega_re.set_xlabel('$k$')
        self._ax_omega_re.set_ylabel('$\omega_{\mathrm{r}}$', color = 'b')
        self._ax_omega_im.set_ylabel('$\gamma$', color = 'r')

        self._ax_omega_re.spines['left'].set_color('b')
        self._ax_omega_re.spines['right'].set_color('r')
        self._ax_omega_re.tick_params(axis = 'y', colors = 'b')
        self._ax_omega_im.tick_params(axis = 'y', colors = 'r')
        #self._ax_omega_re.yaxis.label.set_color('b')
        #self._ax_omega_im.yaxis.label.set_color('r')

        self._ax_omega_re.plot(self._arrk, self._arromega.real, 'b')
        self._ax_omega_im.plot(self._arrk, self._arromega.imag, 'r')
        self._ax_omega_im.axhline(color = '0.5', ls = '--')
        self._ax_omega_re.axvline(x = self._disp.k, color = '0.5', ls = '--')

    def update_plot_modestruct(self):
        '''Update mode structure plot'''
        modestruct = self._disp.get_modestruct(self._ispecies)

        self._ax_modestruct.clear()
        self._ax_modestruct.set_title('Mode structure of $\delta f$')
        self._ax_modestruct_colorbar.clear()
        self._ax_modestruct.set_xlabel('$x$')
        self._ax_modestruct.set_ylabel('$v$')
        nlevel = 64
        fmax = np.max(modestruct[2])
        fmin = np.min(modestruct[2])
        x_mg, v_mg = np.meshgrid(modestruct[0], modestruct[1])
        levels = fmin + (fmax - fmin) * np.arange(nlevel) / (nlevel - 1.0)
        #print modestruct[0].shape, modestruct[1].shape, modestruct[2].shape
        cf = self._ax_modestruct.contourf( \
            x_mg, v_mg, modestruct[2], levels)
        plt.colorbar(cf, cax = self._ax_modestruct_colorbar, \
            format = self._modestruct_colorbar_formatter)

# end of class VisualDispersion

