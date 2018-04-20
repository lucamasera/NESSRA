/**
	Authors: Francesco Asnicar, Luca Masera, Paolo Morettin, Nadir Sella, Thomas Tolio.
	Copyright (C) 2013, all rights reserved

	This file (utility.cpp) is part of the PC++ project.

	You can NOT redistribute it.

	PC++ is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#ifdef _MSC_VER
#define isnan(x) _isnan(x)  // VC++ uses _isnan() instead of isnan()
#endif
#define M_SQRT1_2 0.70710678118654757273731092936941422522068023681640625 // definition for M_SQRT1_2
#include <iostream> // cerr
#include <sstream> // stringstream, stream
#include <fstream> // ifstream, ofstream
#include <algorithm> // atoi, sort, atof
#include <cmath> // M_SQRT1_2, pow
#ifndef _UTILITY
#include "utility.hpp"
#endif
#ifndef _ERF
#include "erf.h"
#endif

using namespace std;

/** Compares two pairs of the type <probe identifier, lookup index>.
	
	@param const intpair &l
	First intpair (formally <p1,l1>).
	
	@param const intpair &r
	Second intpair (formally <p1,l1>).

	@return
	TRUE if p1 is less than p2 (formally p1 < p2). Otherwise FALSE.
*/
bool comparator (const intpair &l, const intpair &r) {
	return l.first < r.first; 
}

/** Reads the file TILE modifying both sizes and data.

	@param const string tilePath
    Path of the file containing the selected tiles.

    @param int* &tilesDim
    Reference of the array containing the lengths of the selected tiles.

	@param intpair** &tiles
    Reference of the matrix of pairs <probe identifier, lookup index>.

    @param int &tileRows
    Reference of the number of selected tiles.

    @return nothing.
*/
bool readTile(const string tilePath, int* &tilesDim, intpair** &tiles, int &tileRows) {
 	ifstream tile(tilePath.c_str());
 	string line;
 	tileRows = 0;

	if(!tile.is_open()) {
		return false;
	}
	
 	//count the numbers of rows
 	while (getline(tile, line)) {
 		tileRows++;
 	}

 	tile.close();

	// create tileRows arrays
 	tilesDim = new int[tileRows];

	// intialize tilesDim structure
 	for (int i = 0; i < tileRows; i++) {
 		tilesDim[i] = 0;
 	}

 	tiles = new intpair*[tileRows];

    tile.open(tilePath.c_str());

    // For each line (thus, tile), reads the number of words in order to know how many nodes there are in the subgraph.
    for (int  i = 0; getline(tile, line); i++) {
	    stringstream stream(line);

	 	while (getline(stream, line, ' ')) {
	 		tilesDim[i]++;
	 	}

	 	// Instantiate a new array of intpair as long as the previous read.
	 	tiles[i] = new intpair[tilesDim[i]];
	}

 	tile.close();

 	tile.open(tilePath.c_str());

 	// For each tile, extracts the index of involved probes and associates it to a (increasing) lookup index.
 	for (int r = 0; r < tileRows; r++) {
 		getline(tile, line);
 		stringstream stream(line);

	 	for (int c = 0; c < tilesDim[r]; c++) {
	 		getline(stream, line, ' ');
	        int temp = atoi(line.c_str());
 			tiles[r][c] = make_pair(temp, c);
	 	}
	}

 	tile.close();

    //sort the indexes of the subgraphs
    for (int i = 0; i < tileRows; i++) {
 		sort(tiles[i], tiles[i] + tilesDim[i]);
 	}
	
	return true;
 }

/** Reads the file CGN saving the biological data that will be used to compute the correlation coefficients.

	@param const string cgnPath
	Path of the file containing the complete gene network.

	@param const intpair* nodesPermutation
	Permutation of the nodes taken into consideration. 

	@param Graph* &g
	The reference of the Graph object representing the gene network.

    @return nothing.
*/
bool readCGN(const string cgnPath, const intpair* nodesPermutation, Graph* &g) {
 	ifstream cgn(cgnPath.c_str());
 	string line;
 	int cols = 0;
	
	if(!cgn.is_open()) {
		return false;
	}

 	getline(cgn, line);
 	cgn.close();
	
 	stringstream stream(line);

    //count the numbers of columns
 	while (getline(stream, line, ',')) {
 		cols++;
 	}

    //decrement cols beacause we discard the first column that contains the probe ids
 	cols--;

    // Initializes the bioData matrix
 	g->nCols = cols;
 	g->bioData = new double*[g->nRows];

 	for (int i = 0; i < g->nRows; i++) {
 		g->bioData[i] = new double[g->nCols];
 	}

    //read the file and save the values in bioData
    cgn.open(cgnPath.c_str());

    // trash the header
    getline(cgn, line);

	//counter for the subset graph (aka, tile)
    int c = 0;

    //i start from 1, as the number of rows shows in the text editors
    for (int i = 0; (c < g->nRows) && getline(cgn, line); i++) {
        //check if it is the rows chosen in readTile()
    	if (i == nodesPermutation[c].first) {
    		stringstream stream(line);

    		for (int j = 0; ((j - 1) < g->nCols) && getline(stream, line, ','); j++) {
    			if (j == 0) { //i'm reading the first token that is the probe id
    				g->probeIDs[nodesPermutation[c].second] = line;
    			} else {
    				g->bioData[nodesPermutation[c].second][j - 1] = atof(line.c_str());
    				g->means[nodesPermutation[c].second] += g->bioData[nodesPermutation[c].second][j - 1];
    			}
    		}

    		g->means[nodesPermutation[c].second] /= (double) g->nCols;
    		c++;
    	}
    }

    cgn.close();
	return true;
}

/** Computes the continous density function.
	M_SQRT1_2 takes value 1/sqrt(2).
	
	@param const double value
	Value for which it will be computed its cumulative normal distribution.

	@return The cumulative normal distribution decimal value for the passed parameter.
*/
double comulativeNormalDistribution(const double value) {
	return 0.5 * a_erfc(-value * M_SQRT1_2);
}

/** Finds the correlation coefficient when l is greater than 1 (formally, when l > 1).
	
	@param const int a
	Index of the selected row, that also represents the "departure" node i of an edge i->j.

	@param const int b
	Index of the selected column, that also represents the "arrival" node j of an edge i->j.

	@param const Graph* g
	The Graph object representing the gene network.

	@param const int* neighbours
	Set of the nearby nodes.

	@param const int* subset
	Subset of the lookup indexes of the nearby nodes.
	Note that this subset has cardinality l.

	@param const int l
	Reached dimension actually taken into account for the subset cardinality.
	In this case l is greater than 1 (formally, l > 1).

	@param double ** p
	Rho value.
	
	@return The decimal value of the computed correlation for the edge a->b depending on the given neighbours' subset.
*/
double correlations(const int a, const int b, const Graph* g, const int* neighbours, const int* subset, const int l, double ** p) {
	int dim = l + 2;

	//initialization of p (looks like rho)
	for (int i = 0; i < dim - 1; i++) {
		for (int j = i + 1; j < dim; j++) {
			int first, second;

			if (i == 0) {
				first = a;
			} else if (i == 1) {
				first = b;
			} else {
				first = neighbours[subset[i - 2]];
			}

			if (j == 1) {
				second = b;
			} else {
				second = neighbours[subset[j - 2]];
			}

			p[i][j] = p[j][i] = g->rho[first][second];
		}
	}

	for (int k = 1; k <= l; k++) {
		for (int i = 0; i <= (l - k); i++) {
			for (int j = i + 1; j < (dim - k); j++) {
				p[i][j] = p[j][i] = (p[i][j] - p[i][dim - k] * p[j][dim - k]) / (sqrt((1 - pow(p[i][dim - k], 2)) * (1 - pow(p[j][dim - k], 2))));
			}
		}
	}

	return p[0][1];
}

/** Checks if a given string (of the form array of chars) whether representing a float number or not.
	
	@param const char* number
	String (or, rather, array of characters) that should represent a decimal number.

	@return TRUE if the string follows the correct format for representing a float number. FALSE otherwise.
*/
bool isFloat(const char* number) {
	bool noMorePoint = false;

	for (int i = 0; number[i] != '\0'; i++) {
		if ((number[i] < '0') || (number[i] > '9')) {
			if (number[i] == '.') {
				if (!noMorePoint) {
					noMorePoint = true;
				} else {
					return false;
				}
			} else {
				return false;
			}
		}
	}

	return true;
}

/** Prints the uncutted edges of the graph in a .csv file.

	@param Graph* g
	The Graph object representing the gene network.

	@param const intpair* nodesPermutation
	Permutation of the nodes taken into consideration. 

	@param const string fileName
	Name of the output file.

	@return nothing.
*/
void fprintEdgesCSV(Graph* g, const intpair* nodesPermutation, const string fileName, const int mode) {
 	ofstream out(fileName.c_str(), ios_base::out);
 	ostringstream strs;
 	
 	// print the "header" (the probeIDs of the tile)
 	strs << "# ";

 	for (int i = 0; i < g->nRows; i++) {
 		strs << g->probeIDs[nodesPermutation[i].second];

 		if (i < (g->nRows - 1)) {
 			strs << ",";
 		}
 	}

	strs << endl;

 	if (mode == 0) {
	 	// print undirected edges
		for (int i = 0; i < g->nRows-1; i++) {
			for (int j = i+1; j < g->nRows; j++) {
				if (g->matrix[nodesPermutation[i].second][nodesPermutation[j].second]) {
					strs << g->probeIDs[nodesPermutation[i].second] << "," << g->probeIDs[nodesPermutation[j].second] << endl;
				}
			}

			out << strs.str();
			strs.str(std::string()); // cleans the stringstream
		}
	} else if (mode == 1) {
		// print directed edges
		for (int i = 0; i < g->nRows; i++) {
			for (int j = i + 1; j < g->nRows; j++) {
				// undirect
				if (g->matrix[nodesPermutation[i].second][nodesPermutation[j].second] &&
					g->matrix[nodesPermutation[j].second][nodesPermutation[i].second]) {
					strs << g->probeIDs[nodesPermutation[i].second] << "\t<->\t" << g->probeIDs[nodesPermutation[j].second] << endl;
				}

				// right direct
				if (g->matrix[nodesPermutation[i].second][nodesPermutation[j].second] &&
					!g->matrix[nodesPermutation[j].second][nodesPermutation[i].second]) {
					strs << g->probeIDs[nodesPermutation[i].second] << "\t-->\t" << g->probeIDs[nodesPermutation[j].second] << endl;
				}

				// left direct
				if (!g->matrix[nodesPermutation[i].second][nodesPermutation[j].second] &&
					g->matrix[nodesPermutation[j].second][nodesPermutation[i].second]) {
					strs << g->probeIDs[nodesPermutation[i].second] << "\t<--\t" << g->probeIDs[nodesPermutation[j].second] << endl;
				}
			}

			out << strs.str();
			strs.str(std::string()); // cleans the stringstream
		}
	}  else if (mode == 2) {
		strs.str(std::string()); // cleans the header, not valid for igraph format

		// print directed edges in igraph format
		for (int i = 0; i < g->nRows; i++) {
			for (int j = 0; j < g->nRows; j++) {
				if (g->matrix[nodesPermutation[i].second][nodesPermutation[j].second]) {
					strs << g->probeIDs[nodesPermutation[i].second] << "\t" << g->probeIDs[nodesPermutation[j].second] << endl;
				}
			}

			out << strs.str();
			strs.str(std::string()); // cleans the stringstream
		}
	} else if (mode == 3) {
		// print directed edges
		for (int i = 0; i < g->nRows; i++) {
			for (int j = i + 1; j < g->nRows; j++) {
				// undirect
				if (g->matrix[nodesPermutation[i].second][nodesPermutation[j].second] &&
					g->matrix[nodesPermutation[j].second][nodesPermutation[i].second]) {

					if(g->double_directed[nodesPermutation[i].second][nodesPermutation[j].second]) {
						strs << g->probeIDs[nodesPermutation[i].second] << "\t<->\t" << g->probeIDs[nodesPermutation[j].second] << endl;
					} else {
						strs << g->probeIDs[nodesPermutation[i].second] << "\t---\t" << g->probeIDs[nodesPermutation[j].second] << endl;
					}
				}

				// right direct
				if (g->matrix[nodesPermutation[i].second][nodesPermutation[j].second] &&
					!g->matrix[nodesPermutation[j].second][nodesPermutation[i].second]) {
					strs << g->probeIDs[nodesPermutation[i].second] << "\t-->\t" << g->probeIDs[nodesPermutation[j].second] << endl;
				}

				// left direct
				if (!g->matrix[nodesPermutation[i].second][nodesPermutation[j].second] &&
					g->matrix[nodesPermutation[j].second][nodesPermutation[i].second]) {
					strs << g->probeIDs[nodesPermutation[i].second] << "\t<--\t" << g->probeIDs[nodesPermutation[j].second] << endl;
				}
			} 
		}
	} else {
		cerr << "[E] Invalid given mode: " << mode << ". No output file written!" << endl;
	}

	out.close();
 }

/** Counts the number of (uncutted) edges in the graph.

	@param bool** matrix
	Matrix of booleans representing a tabular form of the presence and absence of all the edges in the graph.
	The boolean located in the cell of i-th row and j-th column represents the presence/absence of the edge i->j.

	@param const int rows
	The number of rows in the matrix.

	@param const int cols
	The number of rows in the matrix.

	@return The number of uncutted edges.
*/
int countArcs(bool** matrix, const int rows, const int cols, const int mode) {
	int counter = -1;

	if (mode == 0) {
		counter = 0;

		for (int i = 0; i < rows; i++) {
			for (int j = i + 1; j < cols; j++) {
				if (matrix[i][j]) {
					counter++;
				}
			}
		}	
	} else {
		counter = 0;
		
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (matrix[i][j]) {
					counter++;
				}
			}
		}
	}

	return counter;
}
