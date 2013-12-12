#!/usr/bin/env python

# Copyright 2013 Wenjun Deng <wdeng@wdeng.info>
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


'''PIC1D-PETSc distribution data file generator'''

import numpy as np
import argparse
import os
import OutputData

parser = argparse.ArgumentParser( \
    description = 'Generate 2D distribution data file')
parser.add_argument('data_path', metavar = 'data path', \
    nargs = '?', type = str, default = './')
parser.add_argument('-xv', metavar = '<coordinate type index>', \
    nargs = 1, type = int, default = [0], help = '0: xv; 1: v; default: 0')
parser.add_argument('-t', metavar = '<time index>', \
    nargs = 1, type = int, default = [0], help = 'default: 0')
parser.add_argument('-s', metavar = '<species index>', \
    nargs = 1, type = int, default = [0], help = 'default: 0')
parser.add_argument('-d', metavar = '<distribution index>', \
    nargs = 1, type = int, default = [0], help = '0: g; 1: f; 2: delta f; default: 0')
parser.add_argument('-vis', action = 'store_true', \
    help = 'visualization')
args = parser.parse_args()
#print args

fn_pd_xv = 'ptcldist_xv'
fn_pd_v = 'ptcldist_v'
fn_x = 'x'
fn_v = 'v'

fn_ext = '_' + str(args.t[0]) + '_' + str(args.s[0]) + '_' + str(args.d[0])
if os.path.normpath(args.data_path) != '.':
    fn_ext += '_' + os.path.basename(os.path.normpath(args.data_path))
fn_ext += '.dat'

fn_pd_xv += fn_ext
fn_pd_v += fn_ext
fn_x += fn_ext
fn_v += fn_ext

outdat = OutputData.OutputData(args.data_path)

if args.xv[0] == 0:
    pd_xv = outdat.get_ptcldist_xv(args.t[0], args.s[0], args.d[0])
    np.savetxt(fn_pd_xv, pd_xv)
    np.savetxt(fn_x, outdat.x_pd)
    np.savetxt(fn_v, outdat.v_pd)

else:
    pd_v = np.zeros((outdat.nv_pd, 2))
    pd_v[:, 0] = outdat.v_pd
    pd_v[:, 1] = outdat.get_ptcldist_v(args.t[0], args.s[0], args.d[0])

    np.savetxt(fn_pd_v, pd_v)

if args.vis:
    import matplotlib.pyplot as plt
    if args.xv[0] == 0:
        x_mg, v_mg = np.meshgrid(outdat.x_pd, outdat.v_pd)
        plt.contourf(x_mg, v_mg, pd_xv, 64)
        plt.show()

