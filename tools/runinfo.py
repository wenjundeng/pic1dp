#!/usr/bin/env python
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
parser.add_argument('-gr', \
    metavar = ('<lower bound>', '<upper bound>'), \
    help = 'time boundaries for growth rate calculation', \
    nargs = 2, type = float)
parser.add_argument('-sr', \
    metavar = ('<lower bound>', '<upper bound>'), \
    help = 'time boundaries for saturation level calculation', \
    nargs = 2, type = float)
parser.add_argument('-g', \
    metavar = '<# of runs in group>', \
    help = 'get information from a group of runs', \
    nargs = '+', type = int)
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

    # calculate integral of energe over time
    scalar_t = data[irun].get_scalar_t()
    t = scalar_t[0, :]
    eng = scalar_t[1, :]
    if irun == 0:
        t_ref = t
        eng_ref = eng

    intengdt = intfdt(t, eng)
    if data[irun].ntime == data[0].ntime:
        intengdiffdt = intfdt(t_ref, np.abs(eng - eng_ref))
    if irun == 0:
        intengdt_ref = intengdt

    printvalref('int energe dt =', intengdt, intengdt_ref)

    if data[irun].ntime == data[0].ntime:
        print 'int |energe - energe_ref| dt =', intengdiffdt
        print 'int |energe - energe_ref| dt / int energe_ref dt =', \
            intengdiffdt / intengdt_ref * 100.0, '%'

    # calculate growth rate
    if args.gr is not None:
        gamma = data[irun].growthrate_energe_fit(args.gr[0], args.gr[1]) / 2.0
        if irun == 0:
            gamma_ref = gamma
        printvalref('growth rate =', gamma, gamma_ref)
        if args.g is not None:
            gamma_group[irun_group] = gamma

    # calculate saturation level
    if args.sr is not None:
        peak = data[irun].findpeak_energe(args.sr[0], args.sr[1])
        if irun == 0:
            peak_ref = peak
        printvalref('saturation level (energe) =', peak[1], peak_ref[1])
        printvalref('saturation time =', peak[0], peak_ref[0])
        if args.g is not None:
            peak_group[irun_group, :] = peak[:]

    # group operation
    if args.g is not None:
        irun_group += 1
        if irun_group == nrun_group:
            # report group result
            print
            print 'group', igroup,
            if igroup == 0:
                print '(ref)'
            else:
                print

            if args.gr is not None:
                gamma_groupavg = np.mean(gamma_group)
                gamma_groupstd = np.std(gamma_group)
                if igroup == 0:
                    gamma_groupavg_ref = gamma_groupavg
                    gamma_groupstd_ref = gamma_groupstd
                printvalref('growth rate avaerage =', gamma_groupavg, \
                    gamma_groupavg_ref)
                printvalref('growth rate std =', gamma_groupstd, \
                    gamma_groupstd_ref)
            if args.sr is not None:
                peak_groupavg = np.mean(peak_group, axis = 0)
                peak_groupstd = np.std(peak_group, axis = 0)
                if igroup == 0:
                    peak_groupavg_ref = peak_groupavg
                    peak_groupstd_ref = peak_groupstd
                printvalref('saturation level average =', peak_groupavg[1], \
                    peak_groupavg_ref[1])
                printvalref('saturation level std =', peak_groupstd[1], \
                    peak_groupstd_ref[1])
                printvalref('saturation time average =', peak_groupavg[0], \
                    peak_groupavg_ref[0])
                printvalref('saturation time std =', peak_groupstd[0], \
                    peak_groupstd_ref[0])

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

