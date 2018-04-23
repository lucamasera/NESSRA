#!/usr/bin/env python3


__author__ = 'Luca Masera'
__version__ = '0.1'
__date__ = 'April 23. 2018'


import os
import sys
import shutil
import argparse as ap
import subprocess as sb

from pcim import pcim, info, warning, error


p = ap.ArgumentParser(description='NAME and VERSION: pcim.py, ver. '+__version__+' ('+__date__+')'+
    'AUTHOR: '+__author__+
    'DESCRIPTION:')
arg = p.add_argument
arg('-l', '--lgn', type=str, required=True, help='')
arg('-d', '--data', type=str, required=True, help='')
arg('-t', '--tilesizes', nargs='+', type=int, required=True, help='')
arg('-i', '--iterations', nargs='+', type=int, required=True, help='')
arg('-a', '--alphas', nargs='+', type=float, required=True, help='')
arg('-k', '--mc4_list_lenght', type=int, required=True, help='')
arg('-o', '--experiment_name', type=str, required=True, help='')
arg('-p', '--pcpp', type=str, required=True, help='')
arg('-n', '--ncpu', type=int, required=False, default=1, help='')

args = p.parse_args()

if os.path.isdir(args.experiment_name):
    warning('"{}" directory already exists, contents will be removed!'.format(args.experiment_name))
    shutil.rmtree(args.experiment_name)

os.makedirs(args.experiment_name)

for i in args.iterations:
    for a in args.alphas:
        for t in args.tilesizes:
            info('Computing i={} a={} t={}'.format(i,a,t))
            pcim(args.lgn, args.data, t, i, a, os.path.join(args.experiment_name, 'i{}_t{}_a{}'.format(i,t,a)), args.pcpp, args.ncpu)

cmd = ['Rscript', 'rank_expansions.R', args.experiment_name, str(args.mc4_list_lenght)]
t = sb.check_call(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL)