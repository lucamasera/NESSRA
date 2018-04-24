#!/bin/bash

cd pc && make && cd ..

./NESSRA.py -l example_data/lgn_input -d example_data/experiments.csv -t 500 1000 -i 100 200 -a 0.01 0.05 -p pc/bin/pc++ -o test  -n 4 -k 100
