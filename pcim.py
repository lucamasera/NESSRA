#!/usr/bin/env python3


__author__ = 'Francesco Asnicar'
__version__ = '0.3'
__date__ = '24 July 2017'


import os
import sys
import shutil
import time
import math
import random
import argparse as ap
import subprocess as sb
import multiprocessing as mp


def info(s, init_new_line=False):
    if s:
        nfo = '\n' if init_new_line else ''
        nfo += '[i] '
        sys.stdout.write(nfo + str(s) + '\n')
        sys.stdout.flush()


def warning(s, init_new_line=False):
    if s:
        wrn = '\n' if init_new_line else ''
        wrn += '[w] '
        sys.stdout.write(wrn + str(s) + '\n')
        sys.stdout.flush()


def error(s, init_new_line=False):
    if s:
        err = '\n' if init_new_line else ''
        err += '[e] '
        sys.stderr.write(err + str(s) + '\n')
        sys.stderr.flush()


def read_params():
    p = ap.ArgumentParser(description='NAME and VERSION: pcim.py, ver. '+__version__+' ('+__date__+')'+
        'AUTHOR: '+__author__+
        'DESCRIPTION:')
    arg = p.add_argument
    arg('-l', '--lgn', type=str, required=True, help='')
    arg('-d', '--data', type=str, required=True, help='')
    arg('-t', '--tilesize', type=int, required=True, help='')
    arg('-i', '--iterations', type=int, required=True, help='')
    arg('-a', '--alpha', type=float, required=True, help='')
    arg('-o', '--workingdirectory', type=str, required=False, help='')
    arg('-p', '--pcpp', type=str, required=True, help='')
    arg('-n', '--ncpu', type=int, required=False, default=1, help='')

    args = p.parse_args()

    # checks
    if not os.path.isfile(args.lgn):
        error('File not found: '+args.lgn)
        sys.exit(1)

    if not os.path.isfile(args.data):
        error('File not found: '+args.data)
        sys.exit(1)

    if not args.workingdirectory.endswith('/'):
        warning(args.workingdirectory+' added missing /')
        args.workingdirectory += '/'

    if not os.path.isfile(args.pcpp):
        error('File not found: '+args.pcpp)
        sys.exit(1)

    return vars(p.parse_args())


def getLgnProbes(lgn_path, obs_path):
    lgn = []

    with open(lgn_path) as f:
        for row in f.readlines()[1:]:
            lgn += row.strip().split(',')

    lgnNames = set(lgn) # filters repetitions on probe names
    lgnIndices = []

    with open(obs_path) as f:
        rowNumber = 0

        for row in f.readlines()[1:]:
            probeId = row.split(',')[0]

            if probeId in lgnNames:    # the row correspond to a LGN node
                lgnIndices.append(rowNumber)

            rowNumber += 1

    return (lgnIndices, rowNumber)


def subsetToString(subList):
    if len(subList) > 0:
        s = str(subList[0])

        for x in subList[1:]:
            s += ' '+str(x)

        return s
    else:
        warning('Empty "subList"')
        return ''


def fix_lgn_string(lgn_name):
    return lgn_name[:lgn_name.rfind('.')].replace('_', '-') # find and remove the extension and replace '_' with '-'


def pcim(expression_data, lgn_path, alpha, iterations, tile_size, executions_path, ncpu):
    lgnProbes, nProbes = getLgnProbes(lgn_path, expression_data) # gets the LGN indices from the LGN network file, also returns the total number of probes from the 'complete' file
    subsetSize = tile_size - len(lgnProbes) # get the number of genes to add to the LGN

    # check size of subsets
    if(subsetSize <= 0):
        error('Tile dimension minor or equal to lgn size')
        return 1

    allProbes = range(nProbes)
    otherProbes = list(set(allProbes) - set(lgnProbes)) # indices of non-LGN nodes
    workingDirectory = executions_path

    # print parameters file
    with open(os.path.join(workingDirectory, 'parameters.txt'), 'w') as f:
        f.write('\n'.join(['{}: {}'.format(lbl, val) for lbl, val in zip(['PC-IM', 'LGN', 'alpha', 'iterations', 'tile size'], [lgn_path, lgn_path, alpha, iterations, tile_size])]) + '\n')

    lgn_name = lgn_path[(lgn_path.rfind('/') + 1):]
    lgn_name = fix_lgn_string(lgn_name)

    counter = dict() # create a dictionary to store frequencies
    lst_lst = []

    for i in range(iterations):
        random.shuffle(otherProbes)
        probeBag = list(otherProbes) # probes not in the LGN to be extracted

        while probeBag != []:
            subset, probeBag = probeBag[:subsetSize], probeBag[subsetSize:] # takes a subset of non-LGN indices
            toAdd = subsetSize-len(subset)

            if toAdd > 0:
                extra = list(set(otherProbes).difference(set(subset)))
                random.shuffle(extra)
                exxtra = extra[:toAdd] # get exactly extra probes to add
                subset += exxtra

                for key in exxtra:
                    if key in counter:
                        counter[key] += 1
                    else:
                        counter[key] = 1

            subset = subset+lgnProbes # adds the LGN to the subset
            random.shuffle(subset)
            lst_lst.append(subsetToString(subset))

    # write the PC executions
    input_files = []
    npcs = round(len(lst_lst)/float(ncpu))

    for i in range(ncpu):
        output_file = os.path.join(workingDirectory, 'inputs/pc_exe_'+str(i)+'.txt')
        input_files.append(output_file)

        if i == ncpu-1:
            with open(output_file, 'w') as f:
                f.write('\n'.join(lst_lst[i*int(npcs):]) + '\n')
        else:
            with open(output_file, 'w') as f:
                f.write('\n'.join(lst_lst[i*int(npcs):(i+1)*int(npcs)]) + '\n')

    # write frequencies
    with open(os.path.join(workingDirectory, 'frequency.txt'), 'w') as freqFile:
        freqFile.write('\n'.join(['{}\t{}'.format(c, k) for k, c in counter.items()]))

    with open(os.path.join(workingDirectory, 'lgn_nodes.txt'),'w') as lgnFile:
        lgnFile.write('\n'.join([str(p) for p in lgnProbes]))

    return input_files


def pcpp_exec(cmd):
    if not terminating.is_set():
        t0 = time.time()
        # info('Starting "{}"'.format(cmd[2]))
 
        try:
            t = sb.check_call(cmd, stdout=sb.DEVNULL, stderr=sb.DEVNULL)
        except Exception as e:
            terminating.set()
            error('\n    '.join([str(a) for a in [e, type(e), e.args]]), init_new_line=True)
            error('cannot execute command: {}'.format(cmd), init_new_line=True)
            raise

        # info('"{}" (output: "{}") done in {} s'.format(cmd[2], cmd[3], int(time.time()-t0)))
    else:
        terminating.set()


def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def pcpp_execs(inputs, pcpp, expr_data, out_folder, ncpu, alpha, cutoff=0):
    # tasks = ([pcpp, '-i', expression_data, '-t', tile, '-o',  os.path.splitext(os.path.basename(tile))[0], '-a', str(alpha), '-d', working_directory, '-f', '0'] for tile in inputs) # stand-alone PC++
    tasks = ([os.path.abspath(pcpp), '-i', expr_data, '-t',  tile, '-o', os.path.abspath(out_folder)+'/results/th'+str(i), '-a', str(alpha), '-f', '0'] for i, tile in enumerate(inputs)) # BOINC pc++
    terminating = mp.Event()

    with mp.Pool(initializer=initt, initargs=(terminating, ), processes=ncpu) as pool:
        try:
            [_ for _ in pool.imap_unordered(pcpp_exec, tasks, chunksize=1)]
        except Exception as e:
            error('\n    '.join([str(a) for a in [e, type(e), e.args]]), init_new_line=True)
            error('pcpp_execs crashed', init_new_line=True)
            sys.exit(1)


if __name__ == '__main__':
    t0 = time.time()
    pars = read_params()

    # create output directory
    if os.path.isdir(pars['workingdirectory']):
        warning('"{}" directory already exists, contents will be overwritten!'.format(pars['workingdirectory']))
        shutil.rmtree(pars['workingdirectory'])

    os.makedirs(pars['workingdirectory'])
    os.makedirs(pars['workingdirectory']+'/inputs')
    os.makedirs(pars['workingdirectory']+'/results')
 
    with open(pars['data']) as f:
        for line in f:
            col_number = len(line.split(','))-1
            break

    num_exp = 1
    input_files = pcim(os.path.abspath(pars['data']), os.path.abspath(pars['lgn']), pars['alpha'], pars['iterations'], pars['tilesize'], os.path.abspath(pars['workingdirectory']), pars['ncpu'])

    # pcpp_execs(input_files, pars['pcpp'], pars['ncpu'], pars['data'], pars['workingdirectory'], pars['alpha'])
    pcpp_execs(input_files, pars['pcpp'], pars['data'], pars['workingdirectory'], pars['ncpu'], pars['alpha'])
    t1 = time.time()
    info('{}'.format(int(t1-t0)))
    sys.exit(0)
