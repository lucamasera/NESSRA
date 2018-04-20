/**
	Authors: Francesco Asnicar, Luca Masera, Paolo Morettin, Nadir Sella, Thomas Tolio.
	Copyright (C) 2013, all rights reserved

	This file (graph.cpp) is part of the PC++ project.

	You can NOT redistribute it.

	PC++ is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#include <cmath> // pow, sqrt
#include <cstdlib> // srand
#include <algorithm> // random_shufle
#ifndef _GRAPH
#include "Graph.hpp"
#endif

/**
 *  Construction that take as parameter the dimension with which the matrix is build
 */
Graph::Graph(const int dim, const bool undirect) {
	nRows = dim;

	numNeighbours = new int[nRows];

	//create the bool matrix
	matrix = new bool*[nRows];

	for (int i = 0; i < nRows; i++) {
		matrix[i] = new bool[nRows];
	}

	cutMap = new bool*[nRows];

	for (int i = 0; i < nRows; i++) {
		cutMap[i] = new bool[nRows];
	}

	means = new double[nRows];
	probeIDs = new std::string[nRows];
	standardDeviations = new double[nRows];

	//create rho matrix
	rho = new double*[nRows];

	for (int i = 0; i < nRows; i++) {
		rho[i] = new double[nRows];
	}

	// create the sepSet basic structure
	if (!undirect) {
		sepSet = new int**[nRows];

		for (int i = 0; i < nRows; i++) {
			sepSet[i] = new int*[nRows];
		}

		// create the lenSepSet matrix
		lenSepSet = new int*[nRows];

		for (int i = 0; i < nRows; i++) {
			lenSepSet[i] = new int[nRows];
		}

		initializeMatrixNull(sepSet, nRows);
		initializeMatrixZero(lenSepSet, nRows);
	}

	//initialize matrix, numNeighbours and l
	initializeMatrix(matrix, nRows);
	initializeCutMap();
	initializeNeighbours(numNeighbours, nRows);
	initializeZero(means, nRows);
	initializeZero(standardDeviations, nRows);
	
}

/**
 *
 */
Graph::~Graph(void) {
	//empty the memory for matrix
	for (int i = 0; i < nRows; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;

	//empty the memory for cutMap
	for (int i = 0; i < nRows; i++) {
		delete[] cutMap[i];
	}
	delete[] cutMap;

	//empty the memory for bioData
	for (int i = 0; i < nRows; i++) {
		delete[] bioData[i];
	}
	delete[] bioData;

	//empty the memory for means
	delete[] means;

	//empty the memory for probeIDs
	delete[] probeIDs;

	//empty the memory for standardDeviations
	delete[] standardDeviations;

	//empty the memory for numNeighbours
	delete[] numNeighbours;

	//empty the memory for rho
	for (int i = 0; i < nRows; i++) {
		delete[] rho[i];
	}
	delete[] rho;
}


/**
 *  Compute the standard deviations for each node in the graph
 */
void Graph::computeStandardDeviations(void) {
	for (int r = 0; r < nRows; r++) {
		for (int c = 0; c < nCols; c++) {
			standardDeviations[r] += pow((bioData[r][c] - means[r]), 2);
		}

		standardDeviations[r] /= (double) nCols;
		standardDeviations[r] = sqrt(standardDeviations[r]);
	}
}

/**
 *  Compute the correlation coefficient of the base case, and store it in rho
 */
void Graph::computeCorrelations(void) {
	double covariance = 0.0;

	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nRows; j++) {
			covariance = 0.0;

			for (int k = 0; k < nCols; k++) {
				covariance += (bioData[i][k] - means[i]) * (bioData[j][k] - means[j]);
			}

			//divide covariance by nCols
			covariance /= nCols;

			//covariance(i, j) / sigma(i) * sigma(j)
			rho[i][j] = covariance / (standardDeviations[i] * standardDeviations[j]);
		}
	}
}

/**
 *  Initialize the boolean matrix to 'true', but the diagonal, setted to 'false'
 */
void Graph::initializeMatrix(bool** matrix, const int dim) {
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (i != j) {
				matrix[i][j] = true;
			} else {
				matrix[i][j] = false;
			}
		}
	}
}

/**
 *  Initialize the boolean matrix to 'true', but the diagonal, setted to 'false'
 */
void Graph::initializeCutMap() {
	for (int i = 0; i < Graph::nRows; i++) {
		for (int j = 0; j < Graph::nRows; j++) {
			Graph::cutMap[i][j] = false;
		}
	}
}

/**
 *  Initialize the array numNeighbours with the value dim - 1, since the initial graph is connected
 */
void Graph::initializeNeighbours(int* numNeighbours, const int dim) {
	for (int i = 0; i < dim; i++) {
		numNeighbours[i] = dim - 1;
	}
}

/**
 *  Initialize the given array till dim to 0.0
 */
void Graph::initializeZero(double* array, const int dim) {
	for (int i = 0; i < dim; i++) {
		array[i] = 0.0;
	}
}

/**
 *  Initialize the int matrix to 0
 */
void Graph::initializeMatrixZero(int** matrix, const int dim) {
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			matrix[i][j] = 0;
		}
	}
}

/**
 *  Initialize the int 3D matrix to NULL
 */
void Graph::initializeMatrixNull(int*** matrix, const int dim) {
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			matrix[i][j] = NULL;
		}
	}
}
