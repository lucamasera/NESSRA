/**
    Authors: Francesco Asnicar, Luca Masera, Paolo Morettin, Nadir Sella, Thomas Tolio.
    Copyright (C) 2013, all rights reserved

    This file (graph.hpp) is part of the PC++ project.

    You can NOT redistribute it.

    PC++ is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#include <string> // string

#ifndef _GRAPH
#define _GRAPH

/**
 *  This class collect all the information to represent the graph
 */

class Graph {
public:
    Graph(const int, const bool); //constructor : initialize the matrix and the numNeighbours
    ~Graph(void); //destructor : free the allocated memory
    bool** matrix; //matrix that represent the graph
    bool** cutMap; //matrix that represent the graph
	bool** double_directed; //keep track of double directed arcs
    int nRows; //represent the number of rows of matrix (and columns), bioData, means, standardDeviations and numNeighbours
    int nCols; //represent the number of columns of bioData
    double** bioData; //matrix that will contains the data to compute Pearson coefficient for the d-separation test
    double* means; //array that contains the means for each node in the graph
    std::string* probeIDs; //array that contains the name of each probe taken in account
    double* standardDeviations; //array that contains the standard deviations for each node in the graph
    int* numNeighbours; //represents the number of adjacents for each node (thought as column vector)
    double** rho; //represents the correlation matrix
    int*** sepSet; // contains the separation set for each pairs of nodes
    int** lenSepSet; // contains the lenght of the separation set for each pairs of nodes
    void computeStandardDeviations(void); //compute the standard deviation for each node in the graph
    void computeCorrelations(void); //compute the correlation coefficient of the base case, and store it in rho
    void initializeCutMap(); //initialize the boolean matrix to 'true', but the diagonal, setted to 'false'

private:
    void initializeMatrix(bool**, const int); //initialize the boolean matrix to 'true', but the diagonal, setted to 'false'
    void initializeNeighbours(int*, const int); //initialize the array numNeighbours with the value dim-1, since the initial graph is connected
    void initializeZero(double*, const int); //initialize the given array till dim to 0.0
    void initializeMatrixZero(int**, const int); //initialize the int matrix to 0
    void initializeMatrixNull(int***, const int); //Initialize the int 3D matrix to NULL
};
#endif //_GRAPH
