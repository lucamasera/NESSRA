/**
	Authors: Francesco Asnicar, Luca Masera, Paolo Morettin, Nadir Sella, Thomas Tolio.
	Copyright (C) 2013, all rights reserved

	This file (pc.cpp) is part of the PC++ project.

	You can NOT redistribute it.

	PC++ is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
**/

#include <iostream> // cerr, cout
#include <unistd.h> // getopt
#include <cstdlib> // atof, atoi
#include <stdio.h> // EOF windows
#include <sstream> // stringstream
#ifndef _SKELETON
#include "skeleton.hpp"
#endif
#ifndef _UTILITY
#include "utility.hpp"
#endif
#ifndef _CPDAG
#include "cpdag.hpp"
#endif
#ifndef _PC
#include "pc.hpp"
#endif

using namespace std;

/**
 *
 */
void printUsage(void) {
	cout << "PC++" << endl;
	cout << endl;
	cout << "NAME" << endl;
    cout << "\tpc++, perform the PC algorithm" << endl;
    cout << endl;
	cout << "SYNOPSIS" << endl;
	cout << "\tpc++ -i complete_gene_network -t tiles -o output_file [-a alpha" << endl;
	cout << "\t     -d output_directory -f output_format]" << endl;
	cout << endl;
	cout << "DESCRIPTION" << endl;
	cout << "\tBla bla bla" << endl;
    cout << endl;
	cout << "\tThe options are as follows:" << endl;
	cout << endl;
	cout << "\t-i\tSpecifies the input file, in biological case is the observation" << endl;
	cout << "\t\tfile." << endl;
	cout << endl;
	cout << "\t-t\tRepresents the subgraph to which apply the PC algorithm. In" << endl;
	cout << "\t\tparticular it contains the index rows that will be extracted" << endl;
	cout << "\t\tfrom the global observation file. If this file contains more" << endl;
	cout << "\t\tthan one row, then it iterates through all rows, executing a PC" << endl;
	cout << "\t\tfor each one." << endl;
	cout << endl;
	cout << "\t-o\tThe first part of the filename to use for the output file." << endl;
	cout << endl;
	cout << "\t-a\tAlpha value to use in the statistical test. If not provide it" << endl;
	cout << "\t\twill be used the default value 0.05." << endl;
	cout << endl;
	cout << "\t-d\tMust be a folder in which store all the computed results. It" << endl;
	cout << "\t\tmust contains the \'/\' at the end." << endl;
	cout << endl;
	cout << "\t-s\tFlag for the order independent skeleton function. The default" << endl;
	cout << "\t\tvalue is false" << endl;
	cout << endl;
	cout << "\t-f\tSpecified the output format. Values are: 0 for the undirected" << endl;
	cout << "\t\tversion (saved into a .csv file), 1 for the directed version" << endl;
	cout << "\t\t(saved into a .ncol file), 2 for the igraph standard (compatible" << endl;
	cout << "\t\twith R, python and Cytoscape and saved into a .ncol file)." << endl;
	cout << "\t-b\tBlanzieri's orientation method" << endl;
}

/**
 *
 */
int main(int argc, char* argv[]) {
	string* cgn = NULL;
	string* tile = NULL;
	string* outfile = NULL;
	string* directory = NULL;
	stringstream ss;
	double alpha = 0.05;
	int option;
	int tileRows; //number of rows of the tile.txt file
	int* tilesDim = NULL; //array that contains the length of each tiles read
	intpair** tiles = NULL; //array that contains the numbers of the rows that will be extracted and represent the subgraph
	Graph* g = NULL;
	bool blanz = false;
	bool star = false;
	unsigned short int output_format = 1;
	const string undirected_str = "_u.csv";
	const string directed_str = "_d.ncol";

	// check the input parameters
	while ((option = getopt(argc, argv, "i:t:o:a:d:f:bs?")) != EOF) {
		switch (option) {
			case 'i' : // saved the complete gene network file
				cgn = new string(optarg);
				break;
			case 't' : // saved the tiles file
				tile = new string(optarg);
				break;
			case 'o' : // saved the output file name
				outfile = new string(optarg);
				break;
			case 'a' : // check the given alpha value
				if (isFloat(optarg)) {
					alpha = atof(optarg);
				} else {
					alpha = 0.05;
					cerr << "[E] alpha is not a float. Use the default value: " << alpha << endl;
				}

				if (alpha < 0.0 || alpha > 1.0) {
					alpha = 0.05;
					cerr << "[E] alpha is out of the valid range: [0.0, 1.0]. Use the default value: " << alpha << endl;
				}

				break;
			case 'd' : // the user specified a directory in which store the results
				directory = new string(optarg);
				break;
			case 'f' : // the user want to print the directed arcs in the igraph format
				if (isFloat(optarg)) {
					output_format = atoi(optarg);
				} else {
					output_format = 1;
					cerr << "[E] output format is not a number. Use the default value: " << output_format << endl;
				}

				if (output_format > 2) {
					output_format = 1;
					cerr << "[E] output format is out of the valid range: 0 (undirect), 1 (our direct version) and 2 (igraph compatible direct output). Use the default value: " << output_format << endl;
				}

				break;
			case 'b' : // Blanzieri's orientation method
				blanz = true;
				break;
			case 's' : // Blanzieri's orientation method
				star = true;
				break;
			case '?' :
			default :
				printUsage();
				return -1;
		}
	}

	if ((cgn == NULL) || (tile == NULL) || (outfile == NULL)) {
		printUsage(); 
		return -1;
	}

	// read the file with the tiles
	if(!readTile(*tile, tilesDim, tiles, tileRows)) {
		cerr << "[E] File \"" << *tile << "\" does not exist!" << endl;
		return -1;
	}
	
	// loop through the given tiles and execute them
	for (int c = 0; c < tileRows; c++) {
		cout << "#" << (c + 1);

		//create the graph
		g = new Graph(tilesDim[c], ((output_format == 0) ? true : false));

		//extract the subgraph from the complete genes file
		if(!readCGN(*cgn, tiles[c], g)){
			cerr << "[E] File \"" << *cgn << "\" does not exist!" << endl;
			return -1;
		}

		//cout << endl << endl;
		//for (int i=0; i < g->nRows; i++)
		//	cout << g->probeIDs[i] << endl;
		//cout << endl << endl;

		//compute the standard deviations
		g->computeStandardDeviations();
		
		//compute the correlations coefficients
		g->computeCorrelations();

		// invoke the skeleton function (output_format = 0 means undirect, hence no separation set)
		skeleton(g, alpha, star, output_format);

		cout << "\tnodes: " << g->nRows << "\tarcs: ";

		ss.str(std::string()); // cleans the stringstream

		if (directory != NULL) {
			// create the filename as output_dir + outfile + '_u.csv'
			ss << *(directory);
		}

		ss << *(outfile) << "_" << c;

		// output the undirected results in PC++ format
		if (output_format == 0) {
			ss << undirected_str; // create the filename as outfile + '_u.csv'

			cout << countArcs(g->matrix, g->nRows, g->nRows, output_format) << endl;
		} else {
			ss << directed_str; // create the filename as outfile + '_d.ncol'

			// direct the arcs
			extimateCPDAG(g, blanz);

			cout << countArcs(g->matrix, g->nRows, g->nRows, output_format) << endl;

			// sepset
			// for (int r = 0; r < g->nRows; r++) {
			// 	for (int c = 0; c < g->nRows; c++) {

			// 		cout << "[" << r+1 << "] [" << c+1 << "] : "; 

			// 		for (int i = 0; i < g->lenSepSet[r][c]; i++) {
			// 			cout << g->sepSet[r][c][i]+1 << ", ";
			// 		}

			// 		cout << endl;
			// 	}
			// }
		}

		// print the CPDAG
		fprintEdgesCSV(g, tiles[c], ss.str(), output_format);

		if (output_format) {
			//empty the memory for sepSet
			for (int i = 0; i < g->nRows; i++) {
				for (int j = 0; j < g->nRows; j++) {
					delete[] g->sepSet[i][j];
				}
				delete[] g->sepSet[i];
			}
			delete[] g->sepSet;

			//empty the memory for lenSepSet
			for (int i = 0; i < g->nRows; i++) {
				delete[] g->lenSepSet[i];
			}
			delete[] g->lenSepSet;
		}

		delete g;
	}

	//free the memory
	for (int i = 0; i < tileRows; i++) {
		delete[] tiles[i];
	}
	delete[] tiles;

	delete[] tilesDim;

	return 0;
}
