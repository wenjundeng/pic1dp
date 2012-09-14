#!/usr/bin/env python
'''Get various information from run(s),
such as growth rate and saturation level'''

import argparse
import numpy as np

import OutputData

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
parser.add_argument('datapaths', metavar = 'data path', \
    help = 'data path for each run', \
    nargs = '*', type = str, default = ['./'])
args = parser.parse_args()
#print args

data = []
for irun in range(len(args.datapaths)):
    print
    print 'run', irun,
    if irun == 0:
        print '(ref)',
    print ': ', args.datapaths[irun]
    data.append(OutputData.OutputData(args.datapaths[irun]))

    # calculate integral of energe over time
    scalar_t = data[irun].get_scalar_t()
    if irun == 0:
        scalar_t_ref = scalar_t

    intengdt = (scalar_t[0, 1] - scalar_t[0, 0]) * scalar_t[1, 0] \
        + (scalar_t[0, data[irun].ntime - 1] \
            - scalar_t[0, data[irun].ntime - 2]) \
        * scalar_t[1, data[irun].ntime - 1]
    if data[irun].ntime == data[0].ntime:
        intengdiffdt = (scalar_t_ref[0, 1] - scalar_t_ref[0, 0]) \
            * abs(scalar_t[1, 0] - scalar_t_ref[1, 0]) \
            + (scalar_t_ref[0, data[0].ntime - 1] \
                - scalar_t_ref[0, data[0].ntime - 2]) \
            * abs(scalar_t[1, data[irun].ntime - 1] \
                - scalar_t_ref[1, data[0].ntime - 1])
    for itime in range(1, data[irun].ntime - 1):
        intengdt += (scalar_t[0, itime + 1] - scalar_t[0, itime - 1]) \
            * scalar_t[1, itime]
        if data[irun].ntime == data[0].ntime:
            intengdiffdt += (scalar_t_ref[0, itime + 1] \
                    - scalar_t_ref[0, itime - 1]) \
                * abs(scalar_t[1, itime] - scalar_t_ref[1, itime])

    intengdt /= 2.0
    if irun == 0:
        intengdt_ref = intengdt

    print 'int energe dt = ', intengdt, \
        ' diff with ref: ', intengdt - intengdt_ref, \
        '(', (intengdt - intengdt_ref) / intengdt_ref * 100.0, '%)'

    if data[irun].ntime == data[0].ntime:
        intengdiffdt /= 2.0
        print 'int |energe - energe_ref| dt = ', intengdiffdt
        print 'int |energe - energe_ref| dt / int energe_ref dt = ', \
            intengdiffdt / intengdt_ref

    # calculate growth rate
    if args.gr is not None:
        gamma = data[irun].growthrate_energe_fit(args.gr[0], args.gr[1]) / 2.0
        if irun == 0:
            gamma_ref = gamma
        print 'growth rate = ', gamma, \
            ' diff with ref: ', gamma - gamma_ref, \
            '(', (gamma - gamma_ref) / gamma_ref * 100.0, '%)'

    # calculate saturation level
    if args.sr is not None:
        peak = data[irun].findpeak_energe(args.sr[0], args.sr[1])
        if irun == 0:
            peak_ref = peak
        print 'saturation level (energe) = ', peak[1], \
            ' diff with ref: ', peak[1] - peak_ref[1], \
            '(', (peak[1] - peak_ref[1]) / peak_ref[1] * 100.0, '%)'
        print 'saturation time = ', peak[0], \
            ' diff with ref: ', peak[0] - peak_ref[0], \
            '(', (peak[0] - peak_ref[0]) / peak_ref[0] * 100.0, '%)'

