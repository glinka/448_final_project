#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <random>
#include <Eigen/Sparse>
using namespace Eigen;
using namespace std;

int main(int argc, char *argv[]) {
    /* Arguments list:
       [1] : n, number of vertices
       [2] : h/lambda, avg vertex degree
       [3] : number of opinions
       [4] : initial distribution of opinions
       [5] : a/alpha, model parameter
       [6] : rewiring scheme ("same" = 0, "random" = 1)
       [7] : maximum allowable iterations
       [8] : graph data collected every [6] number of iterations
       [9] : output filename
    */
    int n = atoi(argv[1]);
    double avgDeg = atof(argv[2]);
    int k = atoi(argv[3]);
    int i, j;
    double initDist[k];
    double sum = 0;
    //Populate array of initial distributions
    for(i = 0; i < k; i++) {
	initDist[i] = atof(argv[i+4]);
	sum += initDist[i];
    }
    if(sum != 1) {
	cerr << "sum(u) != 1, improper initial distribution" << "\n";
	return 1;
    }
    double a = atof(argv[k+4]);
    string rewireTo = string(argv[k+5]);
    int maxIter = atoi(argv[k+6]);
    int collectionInterval = atoi(argv[k+7]);
    string fileName = string(argv[k+8]);
    cout << fileName << rewireTo << "\n";
    //init matrices
    MatrixXd A = MatrixXd::Zero(n,n);
    MatrixXi Opns = MatrixXi::Zero(n,1);
    double p = avgDeg/(n-1);
    int v = 1;
    int w = -1;
    default_random_engine generator;
    uniform_real_distribution<double> distr(0.0, 1.0);
    while(v <= n) {
	double rv1 = distr(generator);
	w += (int) 1 + floor(log(1-rv1)/log(1-p));
	while((w >= v) && (v <= n)) {
	    double rv2 = distr(generator);
	    double partialSum = 0.0;
	    int indexCounter = 0;
	    while(partialSum < rv2) {
		partialSum += initDist[indexCounter];
		indexCounter++;
	    }
	    Opns(v-1,0) = indexCounter;
	    w = w - v;
	    v++;
	}
	if(v < n) {
	    A(v-1,w-1) = 1;
	    A(w-1,v-1) = 1;
	}
    }
    int totalEdges = 0;
    MatrixXi opnCounts = MatrixXi::Zero(k,1);

    int currentOpn;
    for(int i = 0; i < n; i++) {
	currentOpn = Opns(i,0);
	opnCounts(currentOpn-1,0) += 1;
	for(j = i+1; j < n; j++) {
	    totalEdges += A(i,j);
	}
    }
    int conflicts = 0;
    for(i = 0; i < n; i++) {
	currentOpn = Opns(i,0);
	for(j = i+1; j < n; j++) {
	    if((A(i,j) != 0) && (Opns(j,0) != currentOpn)) {
		conflicts++;
	    }
	}
    }
    int iters = 0;
    MatrixXi N10timeCourse = MatrixXi::Zero((int) maxIter/collectionInterval,0);
    MatrixXi N1timeCourse = MatrixXi::Zero((int) maxIter/collectionInterval, 0);
    MatrixXi stepTimeCourse = MatrixXi::Zero((int) maxIter/collectionInterval, 0);
    int step = 0;
}
