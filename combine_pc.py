#!/usr/bin/env python3

#import numpy as np
import os, sys, shutil
from collections import Counter

pcim_folder = sys.argv[1]
results_path = os.path.join(pcim_folder, 'results')
extra_count = {int(p):int(f) for f,p in [l.strip().split() for l in open(os.path.join(pcim_folder, 'frequency.txt')).readlines()]}
lgn_probes = {int(p.strip()) for p in open(os.path.join(pcim_folder, 'lgn_nodes.txt')).readlines()}

iterations = 5

pc_results = os.listdir(results_path)

edge_count = Counter()

for pc in pc_results:
    with open(os.path.join(results_path, pc)) as f:
        edges = [tuple(sorted(map(int, e.strip().split(','))))  for e in f.readlines()]
    edge_count.update(edges)

#print(sorted(edge_count.items()))

#print(extra_count)

expansion_count = Counter()

for e, c in edge_count.items():
    x = set(e).difference(lgn_probes)
    if len(x) == 1:
        p = x.pop()
        expansion_count.update({p:c,})

expansion = []
for p, c in expansion_count.items():
    expansion.append([p, c, c/((iterations + extra_count.get(p, 0) * len(lgn_probes)))])

with open(os.path.join(pcim_folder, 'nessra.expansion'), 'w') as expansion_file:
    for p, f_abs, f_rel in sorted(expansion, key=lambda x: (x[2], x[1], x[0]), reverse=True):
        expansion_file.write('{},{},{:.6f}\n'.format(p, f_abs, f_rel))
