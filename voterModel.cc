#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {
    /* Arguments list:
       [0] : n, number of vertices
       [1] : h/lambda, avg vertex degree
       [2.5] : number of opinions
       [2] : initial distribution of opinions
       [3] : a/alpha, model parameter
       [4] : rewiring scheme ("same" = 0, "random" = 1)
       [5] : maximum allowable iterations
       [6] : graph data collected every [6] number of iterations
       [7] : output filename
    */
    int n = atoi(argv[0]);
    double h = atof(argv[1]);
    int k = atoi(argv[2]);
    int i, j;
    double kInitDist[k];
    //Populate array of initial distributions
    for(i = 0; i < k; i++) {
	kInitDist[i] = atof(argv[i+2]);
    }
    double a = atof(argv[k+2]);
    double rewireTo = atoi(argv[k+3]);
    int maxIter = atoi(argv[k+4]);
    int collectionInterval = atoi(argv[k+5]);
    char fileName[128];
    strcpy(fileName, argv[argc]);
    out << fileName << "\n";
}
