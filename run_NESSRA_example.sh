#!/bin/bash

cd pc && make && cd ..

./NESSRA.py -l example_data/lgn_input -d example_data/experiments.csv -t 1500 1000 -i 1 2 -a 0.0001 0.0005 -p pc/bin/pc++ -o test  -n 4 -k 100
