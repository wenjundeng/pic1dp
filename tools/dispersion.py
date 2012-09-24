#!/usr/bin/env python
'''Numerically solve for omega (real frequency and growth rate)
using analytic dispersion of 1D electrostatic Vlasov-Poisson plamsa
consisted of (shifted) Maxwellian species'''

import argparse
import numpy as np
import scipy.special as special
import sys

def plasma_z(x):
    '''plasma dispersion Z function'''
    #print x
    return 1.0j * np.sqrt(np.pi) * np.exp(-x**2) \
        * (1.0 - special.erf(-1.0j * x))

def muller(func, x0, x1, x2, functol = 1e-14, xtol = 1e-14, niter_max = 100):
    '''Find complex root using
    [Muller method](http://en.wikipedia.org/wiki/Muller's_method)'''
    iiter = 0
    #print 'abs(x2 - x1) =', abs(x2 - x1)
    #print 'abs(func(x2)) =', abs(func(x2))
    while abs(func(x2)) > functol and abs(x2 - x1) > xtol \
        and iiter < niter_max:
        #print abs(func(x2)), abs(x2 - x1), iiter
        w = (func(x2) - func(x1)) / (x2 - x1) \
            + (func(x2) - func(x0)) / (x2 - x0) \
            - (func(x1) - func(x0)) / (x1 - x0)
        sqrtdelta = np.sqrt(w**2 - 4.0 * func(x2) * ( \
            (((func(x2) - func(x1)) / (x2 - x1)) \
            - (func(x1) - func(x0)) / (x1 - x0)) \
            / (x2 - x0) + 0j))
        if np.abs(w + sqrtdelta) > np.abs(w - sqrtdelta):
            denorm = w + sqrtdelta
        else:
            denorm = w - sqrtdelta
        x0 = x1
        x1 = x2
        x2 = x1 - 2.0 * func(x1) / denorm
        iiter += 1
    #print abs(func(x2)), abs(x2 - x1), iiter
    return x2


class Dispersion:
    '''dispersion relation'''

    #npara_nonspecies = 1
    npara_species = 5

    def __init__(self, species_paras, k):
        self._guess0 = 0.4739 + 0.153j
        self._guess1 = 1.793 + 0.491j
        self._guess2 = 0.9371 + 0.287j

        self.k = k
        self._omega = self._guess2
        self._solved = False

        self.nspecies = len(species_paras) / self.npara_species
        self.species_charge = np.zeros(self.nspecies)
        self.species_mass = np.zeros(self.nspecies)
        self.species_temperature = np.zeros(self.nspecies)
        self.species_density = np.zeros(self.nspecies)
        self.species_v0 = np.zeros(self.nspecies)
        for ispecies in range(self.nspecies):
            self.species_charge[ispecies] \
                = species_paras[ispecies * self.npara_species]
                #+ self.npara_nonspecies]
            self.species_mass[ispecies] \
                = species_paras[ispecies * self.npara_species + 1]
                #+ self.npara_nonspecies + 1]
            self.species_temperature[ispecies] \
                = species_paras[ispecies * self.npara_species + 2]
                #+ self.npara_nonspecies + 2]
            self.species_density[ispecies] \
                = species_paras[ispecies * self.npara_species + 3]
                #+ self.npara_nonspecies + 3]
            self.species_v0[ispecies] \
                = species_paras[ispecies * self.npara_species + 4]
                #+ self.npara_nonspecies + 4]

    def set_k(self, k):
        '''change k'''
        if k != self.k:
            self.k = k
            self._solved = False

    def set_guess(self, guess0 = None, guess1 = None, guess2 = None):
        '''change initial guess'''
        if guess0 is not None:
            self._guess0 = guess0
        if guess1 is not None:
            self._guess1 = guess1
        if guess2 is not None:
            self._guess2 = guess2

    def append_guess(self, arrguess):
        '''append to initial guess list'''
        for guess in arrguess:
            if guess == self._guess0:
                self._guess0, self._guess1, self._guess2 \
                    = self._guess1, self._guess2, self._guess0
            elif guess == self._guess1:
                self._guess1, self._guess2 \
                    = self._guess2, self._guess1
            elif guess == self._guess2:
                pass
            else:
                self._guess0, self._guess1, self._guess2 \
                    = self._guess1, self._guess2, guess

    def dispfunc(self, omega):
        '''dispersion function'''
        D = 1.0
        for ispecies in range(self.nspecies):
            vth2 = self.species_temperature[ispecies] \
                / self.species_mass[ispecies]
            #print 'k =', self.k
            #print 'omega / k / sqrt(2) =', omega / self.k / np.sqrt(2.0)
            zeta = (omega / self.k - self.species_v0[ispecies]) \
                / np.sqrt(2.0 * vth2)
            #print 'zeta =', zeta
            #print 'Z(zeta) =', plasma_z(zeta)
            D += self.species_density[ispecies] \
                * self.species_charge[ispecies]**2 \
                / self.species_mass[ispecies] \
                / (self.k**2 * vth2) \
                * (1.0 + zeta * plasma_z(zeta))
        return D

    def solveomega(self):
        '''solve dispersion relation for omega'''
        if self._solved:
            return self._omega
        else:
            self._omega = muller(self.dispfunc, self._guess0, self._guess1, self._guess2)
            self._solved = True
            self.append_guess([self._omega])
            return self._omega

    def get_modestruct(self, ispecies, \
        v_max = 8.0, nx = 64, nv = 128):
        self.solveomega()
        modestruct = np.zeros((nv, nx + 1))
        for iv in range(nv):
            v = (v_max * 2.0) / (nv - 1.0) * iv - v_max
            for ix in range(nx):
                x = (2.0 * np.pi / self.k) / nx * ix
                if ispecies < self.nspecies:
                    ms_species = self.species_charge[ispecies] \
                        / self.species_temperature[ispecies] \
                        * (v - self.species_v0[ispecies]) \
                        / np.sqrt(2.0 * np.pi \
                            * self.species_temperature[ispecies] \
                            / self.species_mass[ispecies]) \
                        * np.exp(-(v - self.species_v0[ispecies])**2 \
                            / (2.0 * self.species_temperature[ispecies] \
                            / self.species_mass[ispecies]))
                else:
                    ms_species = 0.0
                    for ispecies2 in range(self.nspecies):
                        ms_species += self.species_density[ispecies2] \
                            * self.species_charge[ispecies2] \
                            / self.species_temperature[ispecies2] \
                            * (v - self.species_v0[ispecies2]) \
                            / np.sqrt(2.0 * np.pi \
                                * self.species_temperature[ispecies2] \
                                / self.species_mass[ispecies2]) \
                            * np.exp(-(v - self.species_v0[ispecies2])**2 \
                                / (2.0 * self.species_temperature[ispecies2] \
                                / self.species_mass[ispecies2]))
                #end of else of if ispecies < self.nspecies
                ms_harmonic = 1j / (self._omega - self.k * v) \
                    * np.exp(1j * self.k * x)
                modestruct[iv, ix] = ms_species * ms_harmonic.real * 2.0
            # end of for ix in range(nx)
        # end of for iv in range(nv)

        # periodic boundary condition
        modestruct[:, nx] = modestruct[:, 0]
        mg = np.meshgrid(\
            (2.0 * np.pi / self.k) / nx * np.arange(nx + 1.0), \
            (v_max * 2.0) / (nv - 1.0) * np.arange(nv) - v_max)

        return [mg[0], mg[1], modestruct]
    # end of def get_modestruct()

    def print_komega(self, k = None, omega = None):
        '''print k and omega information'''
        if k is None:
            k = self.k
        if omega is None:
            self.solveomega()
            omega = self._omega
        print 'k = ', k, ': omega =', omega,
        print ' (gamma / omega_r =', omega.imag / omega.real * 100.0, '%)'

# end of class Dispersion

if __name__ == '__main__':
    parser = argparse.ArgumentParser( \
        description = 'Numerically solve for omega ' \
        + '(real frequency and growth rate) using dispersion of ' \
        + '1D electrostatic Vlasov-Poisson plasma consisted of '\
        + '(shifted) Maxwellian species')

    parser.add_argument('para', metavar = 'dispersion parameters', \
        help = 'for each species, input charge Z, ' \
        + 'mass m, temperature T, density n, flow v0 in sequence', \
        nargs = '*', type = float)
    parser.add_argument('-ig', metavar = '<initial guess>', \
        help = 'specify at most three different initial guesses', \
        nargs = '+', type = complex)
    parser.add_argument('-k', \
        help = 'specify one value to calculate for single k value; ' \
        + 'specify two values to calculate for a range of k values; ' \
        + 'specify three values to calculate for a range of k values ' \
        + 'given by the last two values, ' \
        + 'and start scanning from the first value', \
        nargs = '+', type = float, default = [0.5])
    parser.add_argument('-sks', metavar = '<k step size for scanning>', \
        help = 'specify step size of k for scanning', \
        nargs = 1, type = float, default = [0.005])
    parser.add_argument('-vis', \
        help = 'visualization, make plots of omega(k) and mode structure', \
        action = 'store_true')
    args = parser.parse_args()
    #print args

    if len(args.para) < Dispersion.npara_species:
        sys.exit('Error: not enough parameters for at least one species.')

    disp = Dispersion(args.para, args.k[0])
    if args.ig is not None:
        disp.append_guess(args.ig)

    #print disp.dispfunc(1.28505698-0.06612800j)
    #print disp.dispfunc(args.ig[2])
    #omega = muller(disp.dispfunc, args.ig[0], args.ig[1], args.ig[2])
    omega = disp.solveomega()
    disp.print_komega()

    nk = len(args.k)
    nk1 = 0
    if nk > 1:
        if nk == 2:
            nk1 = int((args.k[1] - args.k[0]) / args.sks[0]) + 2
        else:
            nk1 = int((args.k[2] - args.k[0]) / args.sks[0]) + 2
        if nk1 > 0:
            arrk = args.k[0] + np.arange(nk1) * args.sks[0]
            arromega = np.zeros(nk1, dtype = complex)
            arromega[0] = omega
            for ik in range(1, nk1):
                disp.set_k(arrk[ik])
                arromega[ik] = disp.solveomega()
            if nk > 2:
                nk2 = int((args.k[0] - args.k[1]) / args.sks[0]) + 1
                if nk2 > 0:
                    arrk2 = args.k[0] + np.arange(-nk2, 0.0) * args.sks[0]
                    arromega2 = np.zeros(nk2, dtype = complex)
                    #print disp._guess0, disp._guess1, disp._guess2
                    disp.append_guess(arromega[3 : : -1])
                    #print disp._guess0, disp._guess1, disp._guess2
                    for ik in range(nk2 - 1, -1, -1):
                        disp.set_k(arrk2[ik])
                        arromega2[ik] = disp.solveomega()
                    arrk = np.append(arrk2, arrk)
                    arromega = np.append(arromega2, arromega)
            # end of if nk > 2
            if not args.vis:
                print
                for ik in range(len(arrk)):
                    disp.print_komega(arrk[ik], arromega[ik])
        # end of if nk1 > 0
    #end of if nk > 1
    if args.vis:
        import matplotlib.pyplot as plt
        import VisualDispersion
        if nk1 < 1:
            arrk = [args.k[0]]
            arromega = [omega]
        disp.set_k(args.k[0])
        visdisp = VisualDispersion.VisualDispersion( \
            disp, arrk, arromega)
        plt.show()

# end of if __name__ == '__main__'

