/**
	Authors: Francesco Asnicar, Luca Masera, Paolo Morettin, Nadir Sella, Thomas Tolio.
	Copyright (C) 2013, all rights reserved

	This file (skeleton.hpp) is part of the PC++ project.

	You can NOT redistribute it.

	PC++ is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#ifndef _GRAPH
#include "Graph.hpp" // Graph
#endif

#ifndef _SKELETON
#define _SKELETON

// Type definition for the pair of integers.
typedef std::pair<int,int> intpair;

// 
void testAndRemove(const int*, const int*, double, Graph* &, const int, const int, const int, const double, const unsigned short int);

// 
void remove(Graph* &);

// 
double getCorrelationCoefficient(const int*, const int*, const int, Graph* &, const int, const int, double**);

// 
void iterativeComb(int*, const int, const int, Graph* &, const int, const int, const double, double**, int*, const bool, const unsigned short int);

// 
void findAllSubsets(Graph* &, const int, const int, const int, const double, int*, double**, int*, const bool, const unsigned short int);

//
void skeleton(Graph* &, const double, const bool, const unsigned short int);

#endif //_SKELETON
