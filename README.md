# NESSRA: Network Expansion by Stratified Subsetting and Ranking Aggregation
[NESSRA](http://journals.sagepub.com/doi/abs/10.1177/1094342016662508) is a pipeline that finds candidate genes to expand a known gene-network.
NESSRA applies the `skeleton` function of the [`PC-algorithm`](https://www.researchgate.net/publication/242448131_Causation_Prediction_and_Search) on all the disjoint subsets (of size `t`) of the genes in a given organism.
This step is repeated `i` times shuffling the input genes, and the results are combined in an expansion list.
Several expansion lists are produced changing the subset-size, number of iterations, and alpha parameters.
The resulting expansion lists are combined with the MC4 aggregation method.

This code is the standalone implementation of the [gene@home](https://gene.disi.unitn.it/test/) project.

## Requirements
The NESSRA pipeline is implemented in `Python3`, `C++` and `R`.
It requires the [`markovchain`](https://cran.r-project.org/web/packages/markovchain/index.html) `R` package.

The code has been tested with:
- `Ubuntu 18.04`
- `Python 3.6.5`
- `R 3.4.4`
- `gcc 7.3.0`
- `g++ 7.3.0`


## Usage
Compile the source code of the `pc++` applicative.
```
cd pc/ && make && cd ..
```
Here are listed the parameters of NESSRA.
```
-h, --help            show this help message and exit
-l --lgn              file containing the edges of the gene-network to
                      expand;
-d --data             gene-expression experiment data matrix;
-t --tilesizes        one or more values for the subsets-size;
-i --iterations       one or more values for the number of iterations;
-a --alphas           one or more values for the confidence threshold alpha;
-k --mc4_list_lenght  lenght of the list to aggregate with MC4;
-o --experiment_name  output folder containing the experiments;
-p --pcpp             path to the pc++ executable;
-n  --ncpu            number of CPUs to use in parallel.
```

We provide an example script to run an expansion on a subset of vitis gene-expression data.
```
./run_NESSRA_example.sh
```
This may take several minutes.

## Results reprodution
Download the input [dataset](https://gene.disi.unitn.it/test/download/2b9/vv_experiments.csv) for Vitis vinifera.
```
wget https://gene.disi.unitn.it/test/download/2b9/vv_experiments.csv
```
