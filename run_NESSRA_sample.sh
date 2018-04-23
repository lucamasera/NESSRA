#!/bin/bash

cd pc && make && cd ..

./NESSRA.py -l sample_data/lgn_input -d sample_data/experiments.csv -t 500 1000 -i 10 20 -a 0.0001 0.0005 -p pc/bin/pc++ -o test  -n 4 -k 100