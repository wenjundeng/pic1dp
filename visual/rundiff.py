#!/usr/bin/env python
'''Calculate difference of delta f of two runs

The difference is given by: integral |(df_2 - df_1)| dv dx
df_1: delta f of first run
df_2: delta f of second run

Two runs must have the same x-v grid structure in output'''

import argparse
import numpy as np

import visual

parser = argparse.ArgumentParser( \
    description = 'Calculate difference of delta f of two runs')
parser.add_argument('data_path1', metavar = 'data path 1', type = str)
parser.add_argument('data_path2', metavar = 'data path 2', type = str)
args = parser.parse_args()

data1 = visual.OutputData(args.data_path1)
data2 = visual.OutputData(args.data_path2)

ptcldist_xv_1 = data1.get_ptcldist_xv(data1.ntime - 1, 0, 2, False)
ptcldist_xv_2 = data2.get_ptcldist_xv(data2.ntime - 1, 0, 2, False)

dfdiff = np.sum(np.abs(ptcldist_xv_2 - ptcldist_xv_1)) \
    * (data1.lx / data1.nx) * (data1.v_max * 2.0 / (data1.nv - 1))
print 'integral |(df_2 - df_1)| dv dx = ', dfdiff

