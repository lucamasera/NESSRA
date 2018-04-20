/**
	Authors: Francesco Asnicar, Luca Masera, Paolo Morettin, Nadir Sella, Thomas Tolio.
	Copyright (C) 2013, all rights reserved

	This file (utility.hpp) is part of the PC++ project.

	You can NOT redistribute it.

	PC++ is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#ifndef _GRAPH
#include "Graph.hpp"
#endif
#ifndef _SKELETON
#include "skeleton.hpp"
#endif

#ifndef _UTILITY
#define _UTILITY

/**
 *  This library will contain all the functions used by the PC-algorithm
 */

// Compares two pairs of the type <probe identifier, lookup index>.
bool comparator (const intpair &, const intpair &);

// Reads the file TILE modifying both sizes and data.
bool readTile(const std::string, int* &, intpair** &, int &);

// Reads the file CGN saving the biological data that will be used to compute the correlation coefficients.
bool readCGN(const std::string, const intpair*, Graph* &);

// Computes the continous density function.
double comulativeNormalDistribution(const double);

// Finds the correlation coefficient when l is greater than 1 (formally, when l > 1).
double correlations(const int, const int, const Graph*, const int*, const int*, const int, double **);

// Checks if a given string (of the form array of chars) whether representing a float number or not.
bool isFloat(const char*);

// Prints the uncutted edges of the graph in a .csv file.
void fprintEdgesCSV(Graph*, const intpair*, const std::string, const int);

// Counts the number of (uncutted) edges in the graph.
int countArcs(bool**, const int, const int, const int);

#ifdef _MSC_VER
//erfc function for VC++
double erfcWindows(double);
#endif

#endif //_UTILITY
