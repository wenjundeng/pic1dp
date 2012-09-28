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


'''Get various information from run(s),
such as growth rate and saturation level'''

import argparse
import numpy as np

import OutputData


def intfdt(t, f):
    '''calculate [int f(t) dt]'''
    nt = len(t)
    integral = (t[1] - t[0]) * f[0] + (t[nt - 1] - t[nt - 2]) * f[nt - 1]
    for it in range(1, nt - 1):
        integral += (t[it + 1] - t[it - 1]) * f[it]
    integral /= 2.0
    return integral


def printvalref(desc, value, value_ref):
    '''print value and difference from reference value'''
    print desc, value, ' diff with ref:', value - value_ref,
    if value_ref != 0.0:
        print '(', (value - value_ref) / value_ref * 100.0, '%)'
    else:
        print


parser = argparse.ArgumentParser( \
    description = 'Get various information from run(s)')
parser.add_argument('-g', \
    metavar = '<# of runs in group>', \
    help = 'get information from a group of runs', \
    nargs = '+', type = int)
parser.add_argument('-wg', \
    metavar = '<data file>',
    help = 'write group results to a data file', \
    type = str)
parser.add_argument('-gr', \
    metavar = ('<lower bound>', '<upper bound>'), \
    help = 'time boundaries for growth rate calculation', \
    nargs = 2, type = float)
parser.add_argument('-gref', metavar = '<reference growth rate>', \
    help = 'specify reference growth rate instead of using value from' \
    + ' first run and group', nargs = 1, type = float)
parser.add_argument('-sr', \
    metavar = ('<lower bound>', '<upper bound>'), \
    help = 'time boundaries for saturation level calculation', \
    nargs = 2, type = float)
parser.add_argument('datapaths', metavar = 'data path', \
    help = 'data path for each run', \
    nargs = '*', type = str, default = ['./'])
args = parser.parse_args()
#print args

data = []
if args.g is not None:
    igroup = 0
    nrun_group = args.g[igroup]
    irun_group = 0
    if args.gr is not None:
        gamma_group = np.zeros(nrun_group)
    if args.sr is not None:
        peak_group = np.zeros((nrun_group, 2))

for irun in range(len(args.datapaths)):
    print
    print 'run', irun,
    if irun == 0:
        print '(ref)',
    print ': ', args.datapaths[irun]
    data.append(OutputData.OutputData(args.datapaths[irun]))

    # calculate integral of energy over time
    scalar_t = data[irun].get_scalar_t()
    t = scalar_t[0, :]
    eng = scalar_t[1, :]
    if irun == 0:
        t_ref = t
        eng_ref = eng

    intengdt = intfdt(t, eng)
    if irun == 0:
        intengdt_ref = intengdt

    printvalref('int energy dt =', intengdt, intengdt_ref)

    if data[irun].ntime == data[0].ntime:
        intengdiffdt = intfdt(t_ref, np.abs(eng - eng_ref))
        print 'int |energy - energy_ref| dt =', intengdiffdt
        print 'int |energy - energy_ref| dt / int energy_ref dt =', \
            intengdiffdt / intengdt_ref * 100.0, '%'

    # calculate growth rate
    if args.gr is not None:
        gamma = data[irun].growthrate_energy_fit(args.gr[0], args.gr[1]) / 2.0
        if irun == 0:
            if args.gref is not None:
                gamma_ref = args.gref[0]
            else:
                gamma_ref = gamma
        printvalref('growth rate =', gamma, gamma_ref)
        if args.g is not None:
            gamma_group[irun_group] = gamma

    # calculate saturation level
    if args.sr is not None:
        peak = data[irun].findpeak_energy(args.sr[0], args.sr[1])
        if irun == 0:
            peak_ref = peak
        printvalref('saturation level (energy) =', peak[1], peak_ref[1])
        printvalref('saturation time =', peak[0], peak_ref[0])
        if args.g is not None:
            peak_group[irun_group, :] = peak[:]

    # group operation
    if args.g is not None:
        if irun_group == 0:
            t_group = t
            eng_group = np.zeros((nrun_group, data[irun].ntime))
        # end of if irun_group == 0
        if data[irun].ntime == np.shape(eng_group)[1]:
            eng_group[irun_group, :] = eng
        else:
            print 'Warning: # of time steps of ', args.datapaths[irun], \
                'does not fit group head''s. Following energy integral would' \
                + ' not make sense.'
        irun_group += 1
        if irun_group == nrun_group:
            # report group result
            print
            print 'group', igroup,
            if igroup == 0:
                print '(ref)'
            else:
                print

            # report group energy integral
            eng_groupavg = np.mean(eng_group, axis = 0)
            eng_groupstd = np.std(eng_group, axis = 0, ddof = 1)
            intenggavgdt = intfdt(t_group, eng_groupavg)
            intenggstddt = intfdt(t_group, eng_groupstd)
            if igroup == 0:
                intenggavgdt_ref = intenggavgdt
                intenggstddt_ref = intenggstddt
            printvalref('int energy_avg dt =', intenggavgdt, intenggavgdt_ref)
            printvalref('int energy_std dt =', intenggstddt, intenggstddt_ref)
            print 'int energy_std dt / int energy_avg dt =', \
                    intenggstddt / intenggavgdt * 100.0, '%'

            group_result = np.array([igroup, intenggavgdt, intenggstddt])

            # report group growth rate
            if args.gr is not None:
                gamma_groupavg = np.mean(gamma_group)
                gamma_groupstd = np.std(gamma_group, ddof = 1)
                if igroup == 0:
                    if args.gref is not None:
                        gamma_groupavg_ref = args.gref[0]
                    else:
                        gamma_groupavg_ref = gamma_groupavg
                    gamma_groupstd_ref = gamma_groupstd
                printvalref('growth rate avaerage =', gamma_groupavg, \
                    gamma_groupavg_ref)
                printvalref('growth rate std =', gamma_groupstd, \
                    gamma_groupstd_ref)
                print 'growth rate std / average =', \
                    gamma_groupstd / gamma_groupavg * 100.0, '%'
                group_result = np.append(group_result, \
                    [gamma_groupavg, gamma_groupstd])
            # report group saturation level
            if args.sr is not None:
                peak_groupavg = np.mean(peak_group, axis = 0)
                peak_groupstd = np.std(peak_group, axis = 0, ddof = 1)
                if igroup == 0:
                    peak_groupavg_ref = peak_groupavg
                    peak_groupstd_ref = peak_groupstd
                printvalref('saturation level average =', peak_groupavg[1], \
                    peak_groupavg_ref[1])
                printvalref('saturation level std =', peak_groupstd[1], \
                    peak_groupstd_ref[1])
                print 'saturation level std / average =', \
                    peak_groupstd[1] / peak_groupavg[1] * 100.0, '%'
                printvalref('saturation time average =', peak_groupavg[0], \
                    peak_groupavg_ref[0])
                printvalref('saturation time std =', peak_groupstd[0], \
                    peak_groupstd_ref[0])
                print 'saturation time std / average =', \
                    peak_groupstd[0] / peak_groupavg[0] * 100.0, '%'
                group_result = np.append(group_result, \
                    [peak_groupavg[1], peak_groupstd[1]])

            if igroup == 0:
                group_result_all = np.array([group_result])
            else:
                group_result_all = np.append(group_result_all, \
                    np.array([group_result]), axis = 0)

            print '****************************************'
            # move on to next group
            irun_group = 0
            igroup += 1
            if igroup < len(args.g):
                nrun_group = args.g[igroup]
            else:
                nrun_group = args.g[-1]
            if args.gr is not None:
                gamma_group = np.zeros(nrun_group)
            if args.sr is not None:
                peak_group = np.zeros((nrun_group, 2))
        # end of if irun_group = nrun_group:
    # end of if args.g is not None:
# end of for irun in range(len(args.datapaths)):

if args.wg is not None:
    np.savetxt(args.wg, group_result_all)

