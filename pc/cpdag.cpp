/**
	Authors: Francesco Asnicar, Luca Masera, Paolo Morettin, Nadir Sella, Thomas Tolio.
	Copyright (C) 2013, all rights reserved

	This file (cpdag.cpp) is part of the PC++ project.

	You can NOT redistribute it.

	PC++ is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#include <iostream> // cout
#ifndef _CPDAG
#include "cpdag.hpp"
#endif

using namespace std;

void extimateCPDAG(Graph* &g, bool mode) {
	bool change;
	bool** adj_copy;

	// create a copy of the given graph
	adj_copy = new bool*[g->nRows];

	for (int i = 0; i < g->nRows; i++) {
		adj_copy[i] = new bool[g->nRows];

		for (int j = 0; j < g->nRows; j++) {
			adj_copy[i][j] = g->matrix[i][j];
			
		}
	}

	//Blanzieri double direction
	if(mode){
		g->double_directed = new bool*[g->nRows];

		for (int i = 0; i < g->nRows; i++) {
			g->double_directed[i] = new bool[g->nRows];

			for (int j = 0; j < g->nRows; j++) {
				g->double_directed[i][j] = false;
				g->matrix[i][j] = false;
			}
		}
	}

	// R0 : all pairs (i,j) nonadj with a common node k
	for (int r = 0; r < g->nRows; r++) {
		for (int c = 0; c < g->nRows; c++) {
			// r and c are adjacent, but not the same node
			if ((c != r) && (adj_copy[r][c])) {
				for (int k = 0; k < g->nRows; k++) {
					// r is a common neighbour of c and k (c -- r -- k)
					if (adj_copy[r][k] && !adj_copy[c][k] && c != k) {

						bool belong = false;

						// check if r belong to sepSet[k][c], 
						for (int i = 0; i < g->lenSepSet[k][c]; i++) {
							if (g->sepSet[k][c][i] == r) {
								belong = true;
							}
						}

						// or if r belong to sepSet[c][k] 
						for (int i = 0; i < g->lenSepSet[c][k]; i++) {
							if (g->sepSet[c][k][i] == r) {
								belong = true;
							}
						}

						// r does not belong to sepSet[k][c] or sepSet[c][k]
						if (!belong) { // hence direct the arcs c -> r <- k
							if(mode){
								//never reverse and ordered arc
								// c -> r
								g->matrix[c][r] = true;
								// k -> r
								g->matrix[k][r] = true;
							} else {
								// c -> r
								g->matrix[c][r] = true;
								g->matrix[r][c] = false;
								// k -> r
								g->matrix[k][r] = true;
								g->matrix[r][k] = false;
							}

							// cout << "[R0] " << c+1 << " -> " << r+1 << " <- " << k+1 << endl;
						}
					}
				}
			}
		}
	}

	// save the double directed
	if(mode){
		for (int i = 0; i < g->nRows; i++) {
			for (int j = 0; j < g->nRows; j++) {
				
				if(g->matrix[i][j] && g->matrix[j][i]){
					g->double_directed[i][j] = g->double_directed[j][i] = true;
				}
			}
		}
	}

	// save the changes
	for (int i = 0; i < g->nRows; i++) {
		for (int j = 0; j < g->nRows; j++) {
			
			if(mode){			
				//save the arcs double directed not considered by V structures
				if(adj_copy[i][j] && (!g->matrix[i][j] && !g->matrix[j][i])) {
					g->matrix[i][j] = g->matrix[j][i] = true;
				}
			}

			adj_copy[i][j] = g->matrix[i][j];
		}
	}


	// cout << "V step:" << endl;

	// for (int i = 0; i < g->nRows; i++) {
	// 	for (int j = 0; j < g->nRows; j++) {
	// 		cout << g->matrix[i][j] << "\t";
	// 	}

	// 	cout << endl;
	// }

	// cout << endl;

	 do {
	 	change = false;

		// R1 : if there is an arrow r -> c s.t. r and i are nonadj, then c -> i
		for (int c = 0; c < g->nRows; c++) {
			for (int r = 0; r < g->nRows; r++) {
				// look for a directed arc r -> c
				if (adj_copy[r][c] && !adj_copy[c][r]) {
					// look for an undirected arc c -- i
					for (int i = 0; i < g->nRows; i++) {
						//               ( c -- i )             and             ( r nonadj i )
						if ((g->matrix[c][i] && g->matrix[i][c]) && (!g->matrix[r][i] && !g->matrix[i][r])) {
							// hence direct the arc c -> i
							if(mode) {
								if(!g->double_directed[c][i]){
									g->double_directed[c][i] = true;
									change = true;
									// cout << "[R1] " << c+1 <<  " -> " << i+1 << endl;
								}
							} else {
								g->matrix[c][i] = true;
								g->matrix[i][c] = false;
								change = true;
							}
							

							// cout << "[R1] " << c+1 <<  " -> " << i+1 << endl;
						}
					}
				}
			}
		}

		// save the changes
		if (change) {
			for (int i = 0; i < g->nRows; i++) {
				for (int j = 0; j < g->nRows; j++) {
					
					if(mode) {
						//save the double directed
						if((g->double_directed[i][j] && !g->double_directed[j][i]) || (!g->double_directed[i][j] && g->double_directed[j][i])){

							if(!g->double_directed[i][j]){
								g->matrix[i][j] = false;
							} else {
								g->matrix[j][i] = false;
							}

							g->double_directed[i][j] = g->double_directed[j][i] = false;
						}
					}

					adj_copy[i][j] = g->matrix[i][j];
				}
			}
		}


		// R2 : if there is a chain r -> k -> c, then r -> c
		for (int c = 0; c < g->nRows; c++) {
			for (int r = 0; r < g->nRows; r++) {
				// if the arc between r and c in undirected
				if (adj_copy[r][c] && adj_copy[c][r]) {
					// look for all neighbours of r
					for (int k = 0; k < g->nRows; k++) {
						//               ( r -> k )              and              ( k -> c )
						if ((g->matrix[r][k] && !g->matrix[k][r]) && (g->matrix[k][c] && !g->matrix[c][k])) {
							// hence direct the arc r -> c
							if(mode){
								if(!g->double_directed[r][c]) {
									g->double_directed[r][c] = true;
									change = true;
									// cout << "[R2] " << r+1 << " -> " << c+1 << endl;
								} 
							} else {
								g->matrix[r][c] = true;
								g->matrix[c][r] = false;
								change = true;
							}
							

							// cout << "[R2] " << r+1 << " -> " << c+1 << endl;
						}
					}
				}
			}
		}

		// save the changes
		if (change) {
			for (int i = 0; i < g->nRows; i++) {
				for (int j = 0; j < g->nRows; j++) {

					if(mode) {
						//save the double directed
						if((g->double_directed[i][j] && !g->double_directed[j][i]) || (!g->double_directed[i][j] && g->double_directed[j][i])) {

							if(!g->double_directed[i][j]) {
								g->matrix[i][j] = false;
							} else {
								g->matrix[j][i] = false;
							}

							g->double_directed[i][j] = g->double_directed[j][i] = false;
						}
					}

					adj_copy[i][j] = g->matrix[i][j];
				}
			}
		}

		// R3 : if there are two chains r -- k -> c and r -- l -> c with k nonadj l, then r -> c
		for (int c = 0; c < g->nRows; c++) {
			for (int r = 0; r < g->nRows; r++) {
				// if the arc between r and c in undirected
				if (adj_copy[r][c] && adj_copy[c][r]) {
					// look for a chain that looks like r -- k -> c
					for (int k = 0; k < g->nRows; k++) {
						//               ( r -- k )             and              ( k -> c )
						if ((g->matrix[r][k] && g->matrix[k][r]) && (g->matrix[k][c] && !g->matrix[c][k])) {
							// look for another a chain that looks like r -- l -> c
							for (int l = 0; l < g->nRows; l++) {
								// k and l must be nonadj
								if ((k != l) && (!g->matrix[k][l] && !g->matrix[l][k])) {
									//               ( r -- l )             and              ( l -> c )
									if ((g->matrix[r][l] && g->matrix[l][r]) && (g->matrix[l][c] && !g->matrix[c][l])) {
										// hence direct the arc r -> c
										if(mode) {
											if(!g->double_directed[r][c]) {
												g->double_directed[r][c] = true;
												change = true;
												// cout << "[R3] " << r+1 << " -- " << k+1 << " -> " << c+1 << ", " << r+1 << " -- " << l+1 << " -> " << c+1 << ", " << r+1 << " -> " << c+1 << endl;
											}
										} else {
											g->matrix[r][c] = true;
											g->matrix[c][r] = false;
											change = true;
										}
										// cout << "[R3] " << r+1 << " -- " << k+1 << " -> " << c+1 << ", " << r+1 << " -- " << l+1 << " -> " << c+1 << ", " << r+1 << " -> " << c+1 << endl;
									}
								}
							}
						}
					}
				}
			}
		}

		// save the changes
		if (change) {
			for (int i = 0; i < g->nRows; i++) {
				for (int j = 0; j < g->nRows; j++) {

					if(mode) {
						//save the double directed
						if((g->double_directed[i][j] && !g->double_directed[j][i]) || (!g->double_directed[i][j] && g->double_directed[j][i])) {

							if(!g->double_directed[i][j]) {
								g->matrix[i][j] = false;
							} else {
								g->matrix[j][i] = false;
							}

							g->double_directed[i][j] = g->double_directed[j][i] = false;
						}
					}

					adj_copy[i][j] = g->matrix[i][j];
				}
			}
		}

		// R4 : 

		// save the changes
		if (change) {
			for (int i = 0; i < g->nRows; i++) {
				for (int j = 0; j < g->nRows; j++) {
					adj_copy[i][j] = g->matrix[i][j];
				}
			}
		}
	
	} while (change);

	// int dd = 0;
	// cout << "Double directed: " << endl;

	// for (int i = 0; i < g->nRows; i++) {
	// 	for (int j = 0; j < g->nRows; j++) {
	// 		if(g->double_directed[i][j]) 
	// 			dd++;
	// 	}

	// }

	// cout << dd/2 << endl;

	// deallocate adj_copy
	for (int i = 0; i < g->nRows; i++) {
		delete[] adj_copy[i];
	}
	delete[] adj_copy;
}

// //
// // MUST BE REMOVED
// //
// int main(int argc, char* argv[]) {
// 	int dim = 7;
// 	Graph* g = new Graph(dim);

// 	for (int i = 0; i < dim; i++) {
// 		for (int j = 0; j < dim; j++) {
// 			g->matrix[i][j] = false;
// 		}
// 	}

// 	g->matrix[0][1] = true;
// 	g->matrix[0][3] = true;

// 	g->matrix[1][0] = true;
// 	g->matrix[1][2] = true;
// 	g->matrix[1][3] = true;
// 	g->matrix[1][4] = true;
// 	g->matrix[1][6] = true;

// 	g->matrix[2][1] = true;
// 	g->matrix[2][4] = true;
// 	g->matrix[2][5] = true;

// 	g->matrix[3][0] = true;
// 	g->matrix[3][1] = true;

// 	g->matrix[4][1] = true;
// 	g->matrix[4][2] = true;
// 	g->matrix[4][5] = true;

// 	g->matrix[5][2] = true;
// 	g->matrix[5][4] = true;
// 	g->matrix[5][6] = true;

// 	g->matrix[6][1] = true;
// 	g->matrix[6][5] = true;



// 	 // g->matrix[0][1] = true;
// 	 // g->matrix[1][0] = true;

// 	 // g->matrix[1][2] = true;
// 	 // g->matrix[2][1] = true;

// 	 // g->matrix[2][0] = true;
// 	 // g->matrix[0][2] = true;

// 	// cout << "Matrix:" << endl;

// 	// for (int i = 0; i < dim; i++) {
// 	// 	for (int j = 0; j < dim; j++) {
// 	// 		cout << g->matrix[i][j] << "\t";
// 	// 	}

// 	// 	cout << endl;
// 	// }

// 	// cout << endl;

// 	extimateCPDAG(g, 1);

// 	// cout << "Matrix cpdag:" << endl;

// 	// for (int i = 0; i < dim; i++) {
// 	// 	for (int j = 0; j < dim; j++) {
// 	// 		cout << g->matrix[i][j] << "\t";
// 	// 	}

// 	// 	cout << endl;
// 	// }

// 	return 0;
// }
