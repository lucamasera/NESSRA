#!/bin/bash

N_CPUS=4

cd pc && make && cd ..

wget https://gene.disi.unitn.it/test/download/2b9/vv_experiments.csv -O frontiers_paper_data/vv_experiments.csv

./NESSRA.py -l frontiers_paper_data/input/aba.csv -d frontiers_paper_data/vv_experiments.csv -t 500 1000 -i 1000 2000 -a 0.05 -p pc/bin/pc++ -o results/aba -n $N_CPUS -k 1000 -c
./NESSRA.py -l frontiers_paper_data/input/anthocyanins.csv -d frontiers_paper_data/vv_experiments.csv -t 500 1000 -i 1000 2000 -a 0.05 -p pc/bin/pc++ -o results/anthocyanins -n $N_CPUS -k 1000 -c
./NESSRA.py -l frontiers_paper_data/input/erf.csv -d frontiers_paper_data/vv_experiments.csv -t 500 1000 -i 1000 2000 -a 0.05 -p pc/bin/pc++ -o results/erf -n $N_CPUS -k 1000 -c
./NESSRA.py -l frontiers_paper_data/input/stilbenoids.csv -d frontiers_paper_data/vv_experiments.csv -t 500 1000 -i 1000 2000 -a 0.05 -p pc/bin/pc++ -o results/stilbenoids -n $N_CPUS -k 1000 -c
