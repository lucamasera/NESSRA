/**
	Authors: Francesco Asnicar, Luca Masera, Paolo Morettin, Nadir Sella, Thomas Tolio.
	Copyright (C) 2013, all rights reserved

	This file (skeleton.cpp) is part of the PC++ project.

	You can NOT redistribute it.

	PC++ is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#include <cmath> // isnan, sqrt, log, abs, pow
#include <iostream> // cerr, cout
#include <cstdlib> // atof
#ifdef _MSC_VER
#define isnan(x) _isnan(x)  // VC++ uses _isnan instead of isnan
#endif
#ifndef _GRAPH
#include "Graph.hpp"
#endif
#ifndef _UTILITY
#include "utility.hpp"
#endif
#ifndef _SKELETON
#include "skeleton.hpp"
#endif

using namespace std;

/** 
 *
 */
void testAndRemove(const int* neighbours, const int* subset, double correlationCoefficient, Graph* &g, const int r, const int c, const int l,
	const double alpha, const bool star, const unsigned short int undirect) {
	double pVal;
	double const cutAt = 0.9999999;
	bool NAdelete = true;

	if (isnan(correlationCoefficient)) {
		correlationCoefficient = 0;
	}

	correlationCoefficient = min(cutAt, max(-cutAt, correlationCoefficient));

	pVal = sqrt(g->nCols - l - 3.0) * 0.5 * log((1 + correlationCoefficient) / (1 - correlationCoefficient));

	if (isnan(pVal)) {
		pVal = 0.0;
	}

	pVal = 2 * (1 - comulativeNormalDistribution(abs(pVal)));

	if (isnan(pVal)) {
		pVal = NAdelete ? 1.0 : 0.0;
	}

	//test d-separation
	if (pVal >= alpha) {
		
		// output_format = 0 means undirect, hence no separation set
		if (undirect) {
			// save the size of the sepset at position r, c
			g->lenSepSet[r][c] = l;

			if (l > 0) {
				//create the vector for the sepset at position r, c
				g->sepSet[r][c] = new int[l];
				
				//save the sepset at position r, c
				for (int pos = 0; pos < l; pos++){
					//do not save on position c, r for saving memory
					g->sepSet[r][c][pos] = neighbours[subset[pos]];
				}
			}
		}

		if (star){
			// mark edges
			g->cutMap[r][c] = g->cutMap[c][r] = true;
		} else {
			//remove edges
			g->matrix[r][c] = g->matrix[c][r] = false;

			//decrement neighbours
			g->numNeighbours[r]--;
			g->numNeighbours[c]--;
		}
		
	}
}

/** Cuts all the edge we marked as to cut.
	
	@param Graph* &g
	The reference of the Graph object representing the gene network.
 */
void remove(Graph* &g) {
    for (int r = 0; r < g->nRows - 1; r++) {
        for (int c = r + 1; c < g->nRows; c++) {
            if (g->cutMap[r][c]) {
                g->numNeighbours[r]--;
                g->numNeighbours[c]--;
                g->matrix[r][c] = g->matrix[c][r] = false;
            }
        }
    }
}

/**
 *
 */
double getCorrelationCoefficient(const int* neighbours, const int* subset, const int l, Graph* &g, const int r, const int c, double** p) {
	double correlationCoefficient;

	if (l == 2) {
		double rhoRC, rhoRH1, rhoCH1;
		int h1 = neighbours[subset[0]];
		int h2 = neighbours[subset[1]];

		rhoRC = (g->rho[r][c] - (g->rho[r][h2] * g->rho[c][h2])) / sqrt((1 - pow(g->rho[r][h2], 2)) * (1 - pow(g->rho[c][h2], 2)));
		rhoRH1 = (g->rho[r][h1] - (g->rho[r][h2] * g->rho[h1][h2])) / sqrt((1 - pow(g->rho[r][h2], 2)) * (1 - pow(g->rho[h1][h2], 2)));
		rhoCH1 = (g->rho[c][h1] - (g->rho[c][h2] * g->rho[h1][h2])) / sqrt((1 - pow(g->rho[c][h2], 2)) * (1 - pow(g->rho[h1][h2], 2)));

		correlationCoefficient = (rhoRC - (rhoRH1 * rhoCH1)) / sqrt((1 - pow(rhoRH1, 2)) * (1 - pow(rhoCH1, 2)));
	} else if (l == 3) {
		double rhoRC_H2H3, rhoRH1_H2H3, rhoCH1_H2H3, rhoRC_H3, rhoRH1_H3, rhoRH2_H3, rhoCH1_H3, rhoCH2_H3, rhoH1H2_H3;
		int h1 = neighbours[subset[0]];
		int h2 = neighbours[subset[1]];
		int h3 = neighbours[subset[2]];

		rhoRC_H3 = (g->rho[r][c] - (g->rho[r][h3] * g->rho[c][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[c][h3], 2)));
		rhoRH1_H3 = (g->rho[r][h1] - (g->rho[r][h3] * g->rho[h1][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[h1][h3], 2)));
		rhoRH2_H3 = (g->rho[r][h2] - (g->rho[r][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));
		rhoCH1_H3 = (g->rho[c][h1] - (g->rho[c][h3] * g->rho[h1][h3])) / sqrt((1 - pow(g->rho[c][h3], 2)) * (1 - pow(g->rho[h1][h3], 2)));
		rhoCH2_H3 = (g->rho[c][h2] - (g->rho[c][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[c][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));
		rhoH1H2_H3 = (g->rho[h1][h2] - (g->rho[h1][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[h1][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));

		rhoRC_H2H3 = (rhoRC_H3 - (rhoRH2_H3 * rhoCH2_H3)) / sqrt((1 - pow(rhoRH2_H3, 2)) * (1 - pow(rhoCH2_H3, 2)));
		rhoRH1_H2H3 = (rhoRH1_H3 - (rhoRH2_H3 * rhoH1H2_H3)) / sqrt((1 - pow(rhoRH2_H3, 2)) * (1 - pow(rhoH1H2_H3, 2)));
		rhoCH1_H2H3 = (rhoCH1_H3 - (rhoCH2_H3 * rhoH1H2_H3)) / sqrt((1 - pow(rhoCH2_H3, 2)) * (1 - pow(rhoH1H2_H3, 2)));

		correlationCoefficient = (rhoRC_H2H3 - (rhoRH1_H2H3 * rhoCH1_H2H3)) / sqrt((1 - pow(rhoRH1_H2H3, 2)) * (1 - pow(rhoCH1_H2H3, 2)));
	} else {
		correlationCoefficient = correlations(r, c, g, neighbours, subset, l, p);
	}

	return correlationCoefficient;
}

/**
 *
 * Thanks to Ilya Bursov, http://stackoverflow.com/questions/19327847/n-choose-k-for-large-n-and-k
 */
void iterativeComb(int* neighbours, const int neighboursDim, const int l, Graph* &g, const int r, const int c, const double alpha, double** p,
	int* currentCombination, const bool star, const unsigned short int undirect) {
	double coeff;

	for (int i = 0; i < neighboursDim; i++) {
		currentCombination[i] = i;
	}

	currentCombination[l - 1] = l - 1 - 1;

	do {
		if (currentCombination[l - 1] == (neighboursDim - 1)) {
			int i = l - 1 - 1;

			while (currentCombination[i] == (neighboursDim - l + i)) {
				i--;
			}

			currentCombination[i]++;

			for (int j = i + 1; j < l; j++) {
				currentCombination[j] = currentCombination[i] + j - i;
			}
		} else {
			currentCombination[l - 1]++;
		}

		coeff = getCorrelationCoefficient(neighbours, currentCombination, l, g, r, c, p);
		testAndRemove(neighbours, currentCombination, coeff, g, r, c, l, alpha, star, undirect);
	} while (!g->cutMap[r][c] && g->matrix[r][c] && !((currentCombination[0] == (neighboursDim - l)) && (currentCombination[l - 1] == (neighboursDim - 1))));
}

/** Finds all the subsets adj(i)\{j} with cardinality equals to l (formally, |adj(i)\{j}| = l).
	
	@param Graph* &g
	The reference of the Graph object representing the gene network.

	@param const int i
	Index of the selected row, that also represents the "departure" node i of an edge i->j.

	@param const int j
	Index of the selected column, that also represents the "arrival" node j of an edge i->j.

	@param const int l
	Reached dimension actually taken into account for the subset cardinality.
	In this case l is greater than 1 (formally, l > 1).

	@param bool &hasWorked


	@param const double alpha
	Complement of the confidence interval

	@param int* neighbours
	Set of the indexes of the nearby nodes of the given edge i->j

	@param double** p
	Rhos.

	@return nothing.
*/
void findAllSubsets(Graph* &g, const int i, const int j, const int l, const double alpha, int* neighbours, double** p,
	int* currentCombination, const bool star, const unsigned short int undirect) {
	int neighboursDim = 0;

	if (l == 0) {
		testAndRemove(NULL, NULL, g->rho[i][j], g, i, j, l, alpha, star, undirect);
	} else if (l == 1) {
		//find neighbours (when l > 0) of i without considering j
		for (int k = 0; (g->matrix[i][j]) && !(g->cutMap[i][j]) && (k < g->nRows) && (neighboursDim < g->numNeighbours[i]); k++) {
			if (g->matrix[i][k] && (k != j)) {
				neighboursDim++;
				double correlationCoefficient = (g->rho[i][j] - (g->rho[i][k] * g->rho[j][k])) / sqrt((1 - pow(g->rho[i][k], 2)) * (1 - pow(g->rho[j][k], 2)));
				int pos[1];
				currentCombination[0] = k;
				pos[0] = 0;
				testAndRemove(currentCombination, pos, correlationCoefficient, g, i, j, l, alpha, star, undirect);
			}
		}
	} else {
		//find neighbours (when l > 1) of i without considering j
		for (int k = 0; (k < g->nRows) && (neighboursDim < g->numNeighbours[i]); k++) {
			if ((g->matrix[i][k]) && (k != j)) {
				neighbours[neighboursDim++] = k;
			}
		}

		//look for all subset of length l
		iterativeComb(neighbours, neighboursDim, l, g, i, j, alpha, p, currentCombination, star, undirect);
	}
}

/**
 *
 */
void skeleton(Graph* &g, const double alpha, const bool star, const unsigned short int undirect) {
	int l = -1;
	bool hasWorked = true; //boolean to see that there is at least an arc i,j s.t. |ad j(C, i)\{j}| >= l TO CHECK FROM THE TEXT
	int* neighbours = new int[g->nRows]; //alloc the neighbours array (save time)
	int* currentCombination = new int[g->nRows]; //alloc an array for saving time
	double** p = new double*[g->nRows]; //alloc the rho for correlations() (save time)

	for (int i = 0; i < g->nRows; i++) {
		p[i] = new double[g->nRows];
	}

	//PC-algorithm
	while ((hasWorked) && (l < g->nRows)) {
		l++;
		hasWorked = false;

		for (int i = 0; i < g->nRows; i++) {
			for (int j = 0; j < g->nRows; j++) {
				//check if exists the arc between i and j
				if (g->matrix[i][j] && (g->numNeighbours[i] > l)) {
					hasWorked = true;
					findAllSubsets(g, i, j, l, alpha, neighbours, p, currentCombination, star, undirect);
				}
			}
		}

		if (star) {
            remove(g);
            g->initializeCutMap();
        }
	}

	cout << "\tl: " << l;

	// free the memory		
	delete[] neighbours;

	for (int i = 0; i < g->nRows; i++) {
		delete[] p[i];
	}
	delete[] p;

	delete[] currentCombination;
}
