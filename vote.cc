#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <Eigen/Sparse>
#include <chrono>
#include <iomanip>
#include "vote.h"

using namespace Eigen;
using namespace std;
  
votingModel::votingModel(int n, int k, int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, string rewireTo, string fileName): ROUND_CONST(0.01), n(n), A(MatrixXi::Zero(n,n)), Opns(MatrixXi::Zero(n,1)), k(k), maxIter(maxIter), collectionInterval(collectionInterval), a(a), avgDeg(avgDeg), initDist(initDist), rewireTo(rewireTo), fileName(fileName) {};

/**
   initializes and runs voting model to frozen state or maxIter
**/  
int votingModel::vote() {
    /* 
       n: number of vertices
       k: number of opinions
       avgDeg: lambda/avg vertex degree
       initDist: initial distribution of opinions
       a: alpha, model parameter
       rewireTo: rewiring scheme ("same" = 0, "random" = 1)
       maxIter: maximum allowable iterations
       collectionInterval: graph data collected every "collectionInterval" number of iterations
       fileName: output filename
    */
    int i, j;
    double sum = 0;
    //check that initial distribution has been properly specified
    for(i = 0; i < k; i++) {
	sum += initDist[i];
    }
    if(sum != 1) {
	cerr << "sum(u) != 1, improper initial distribution" << "\n";
	return 1;
    }
    //init graph and uniform random number generator (mersenne twister)
    initGraph();
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 mt1(seed);
    //adding one to ensure the distribution lies in [0,1)
    double normalization = (double) (mt1.max()+1);
    //init matrices
    int totalEdges = 0;
    int currentOpn;
    MatrixXi opnCounts = MatrixXi::Zero(k,1);
    //counts initial numbers of each opinion
    for(int i = 0; i < n; i++) {
	currentOpn = Opns(i,0);
	opnCounts(currentOpn-1,0) = opnCounts(currentOpn-1,0) + 1;
	for(j = i+1; j < n; j++) {
	    totalEdges += A(i,j);
	}
    }
    //count initial number of conflicts
    /**
       Initially the code calculated the number of conflicts from scratch at each step (i.e from
       the adjacency matrix and opinion vector), but this was prohibitively expensive
    **/
    int conflicts = 0;
    for(i = 0; i < n; i++) {
	currentOpn = Opns(i,0);
	for(j = i+1; j < n; j++) {
	    if((A(i,j) != 0) && (Opns(j,0) != currentOpn)) {
		conflicts++;
	    }
	}
    }
    /**
       iters: number of steps the simulation has executed
       N10timeCourse: stores the number of discordant edges
       minorityOpnTimeCourse: stores the number of vertices holding the minority opinion
       V: during each step, stores the index of potential neighbors to chosenVertex (only used in rewire-to-same)
       ***********
       chosenVertex: randomly chosen vertex, could be the source of discrepancies given I actually do not choose an
       edge, but rather a vertex, which result in a different distribution of selected vertices
       ***********
       neighborIndex: temporary array index of the chosenVertex's neighbor in V
       nNeighbors: deg(chosenVertex)
       neighborNumber: similar function to neighborIndex, specifies chosenVertex's neighbor
       neighbor: index of chosenVertex's neighbor
       conflictCounter: tracks changes number of conflicts in system at each step
       newNeighbor: index of new neighbor of chosenVertex (if rewiring)
       collectionStep: number of times data has been collected
       actionToPerform: uniform random number on [0,1), determines what to do with the edge
    **/
    int iters = 0;
    int chosenVertex, chosenEdge, neighborIndex, neighbor, conflictCounter, newNeighbor, collectionStep, edgeCount;
    double actionToPerform;
    MatrixXi N10timeCourse = MatrixXi::Zero((int) maxIter/collectionInterval, 1);
    MatrixXi minorityOpnTimeCourse = MatrixXi::Zero((int) maxIter/collectionInterval, 1);
    MatrixXi stepTimeCourse = MatrixXi::Zero((int) maxIter/collectionInterval, 1);
    MatrixXi V;
    if(rewireTo == "random") {
	while((conflicts > 0) && (iters < maxIter)) {
	    //chose vertex at random, ensure it's degree is nonzero
	    //is choosing an element of A(i,j) equivalent to choosing an edge?
	    //I believe yes, as it simply scales everthing by 2, doesn't change
	    //probability to choose a given edge
	    chosenEdge = ((int) floor(2*totalEdges*(mt1()/normalization))) + 1;
	    i = 0;
	    edgeCount = 0;
	    while(edgeCount < chosenEdge) {
		edgeCount += degs[i++];
	    }
	    edgeCount -= degs[--i];
	    j = 0;
	    while(edgeCount < chosenEdge) {
		edgeCount += A(i,j++);
	    }
	    j--;
	    if(mt1()/normalization > 0.5) {
		chosenVertex = i;
		neighbor = j;
	    }
	    else {
		chosenVertex = j;
		neighbor = i;
	    }
	    //test edges existence
	    if(A(i,j) == 0) {
		cout << "nonexistent edge, error" << endl;
	    }
	    actionToPerform = mt1()/normalization;
	    conflictCounter = 0;
	    //case one, simply change opinion of chosenVertex to match neighbor
	    if(actionToPerform > a) {
		//do nothing if opinions already match
		if(Opns(chosenVertex,0) != Opns(neighbor,0)) {
		    //adjust opinions and opinion counts
		    opnCounts(Opns(chosenVertex,0) - 1,0) = opnCounts(Opns(chosenVertex,0) - 1,0) - 1;
		    opnCounts(Opns(neighbor,0) - 1,0) = opnCounts(Opns(neighbor,0) - 1,0) + 1;
		    Opns(chosenVertex,0) = Opns(neighbor,0);
		    //adjust conflicts
		    for(j = 0; j < n; j++) {
			if((A(chosenVertex,j) != 0) && (Opns(j,0) != Opns(chosenVertex,0))) {
			    conflictCounter++;
			}
			else if(A(chosenVertex,j) != 0) {
			    conflictCounter--;
			}
		    }
		    conflicts += conflictCounter;
		}
	    }
	    //case two, rewire edge
	    else {
		if(Opns(chosenVertex,0) != Opns(neighbor,0)) {
		    conflictCounter--;
		}
		//disconnect edge and find new neighbor from all vertices, ensuring no loops or
		//parallel edges are formed
		A(chosenVertex,neighbor) = 0;
		A(neighbor,chosenVertex) = 0;
		degs[neighbor] = degs[neighbor] - 1;
		newNeighbor = (int) floor(n*(mt1()/normalization));
		while((A(chosenVertex,newNeighbor) != 0) || (newNeighbor == chosenVertex)) {
		    newNeighbor = (int) floor(n*(mt1()/normalization));
		}
		A(chosenVertex,newNeighbor) = 1;
		A(newNeighbor,chosenVertex) = 1;
		degs[newNeighbor] = degs[newNeighbor] + 1;
		if(Opns(chosenVertex,0) != Opns(newNeighbor,0)) {
		    conflictCounter++;
		}
		conflicts += conflictCounter;
	    }
	
	
	    //collect data every collectionInterval steps
	    if(iters % collectionInterval == 0) {
		collectionStep = iters/collectionInterval;
		minorityOpnTimeCourse(collectionStep,0) = opnCounts(0,0)<opnCounts(1,0)?opnCounts(0,0):opnCounts(1,0);
		N10timeCourse(collectionStep,0) = conflicts;
		stepTimeCourse(collectionStep,0) = collectionStep + 1;
	    }
	    iters++;
	}
    }

    /**

*********************************************************************************
Rewire to same is not functional in chooseEdge branch as of 07/04
*********************************************************************************

**/


    else if(rewireTo == "same") {
	while((conflicts > 0) && (iters < maxIter)) {
	    //find vertex with nonzero degree
	    chosenVertex = (int) floor(n*(mt1()/normalization));
	    while(A.block(chosenVertex,0,1,n).sum() == 0) {
		chosenVertex = (int) floor(n*(mt1()/normalization));
	    }
	    //assemble list of possible neighbors, storing their indices in V
	    V = MatrixXi::Zero(n,1);
	    i = 0;
	    for(j = 0; j < n; j++) {
		if((A(chosenVertex,j) != 0) && (Opns(chosenVertex,0) != Opns(j,0))) {
		    V(i,0) = j;
		    i++;
		}
	    }
	    neighborIndex = (int) floor((i)*(mt1()/normalization));
	    neighbor = V(neighborIndex,0);
	    /**	
		V(neighborIndex,0) = V(i,0);
		while((i > 0) && (Opns(chosenVertex,0) == Opns(neighbor,0))) {
		neighborIndex = (int) floor((i--)*(mt1()/normalization));
		neighbor = V(neighborIndex,0);
		V(neighborIndex,0) = V(i,0);
		}
	    **/
	    //check that chosenVertex has a disagreeing neighbor, else do nothing
	    if(i > 0) {
		actionToPerform = mt1()/normalization;
		conflictCounter = 0;
		//case one, change chosenVertex's opinion to match that of neighbor
		if(actionToPerform > a) {
		    //adjust conflicts and opinion-counts
		    opnCounts(Opns(chosenVertex,0) - 1,0) = opnCounts(Opns(chosenVertex,0) - 1,0) - 1;
		    opnCounts(Opns(neighbor,0) - 1,0) = opnCounts(Opns(neighbor,0) - 1,0) + 1;
		    Opns(chosenVertex,0) = Opns(neighbor,0);
		    for(j = 0; j < n; j++) {
			if((A(chosenVertex,j) != 0) && (Opns(j,0) != Opns(chosenVertex,0))) {
			    conflictCounter++;
			}
			else if(A(chosenVertex,j) != 0) {
			    conflictCounter--;
			}
		    }
		    conflicts += conflictCounter;
		}
		//case two, rewire edge
		else {
		    //remove current edge and find new edge from those vertices that share chosenVertex's opini
		    conflictCounter--;
		    A(chosenVertex,neighbor) = 0;
		    A(neighbor,chosenVertex) = 0;
		    newNeighbor = (int) floor(n*(mt1()/normalization));
		    /** notes to self:
		     //in this, the "rewireTo == same" section, the probability of 
		     //finding a neighbor at random decreases as we must ensure not only that
		     //the edge isn't parallel/loop but also that the opinions match
		     //could form list of matching opinion edges, but this would be
		     //computationally expensive for large n, and would probably be slower
		     //than current, random implementation
		     **/
		    while((Opns(chosenVertex,0) != Opns(newNeighbor,0)) || (A(chosenVertex,newNeighbor) != 0) || (newNeighbor == chosenVertex)) {
			newNeighbor = (int) floor(n*(mt1()/normalization));
		    }
		    A(chosenVertex,newNeighbor) = 1;
		    A(newNeighbor,chosenVertex) = 1;
		    conflicts += conflictCounter;
		}
	    }
	    if(iters % collectionInterval == 0) {
		collectionStep = iters/collectionInterval;
		minorityOpnTimeCourse(collectionStep,0) = opnCounts(0,0)<opnCounts(1,0)?opnCounts(0,0):opnCounts(1,0);
		N10timeCourse(collectionStep,0) = conflicts;
		stepTimeCourse(collectionStep,0) = collectionStep + 1;
	    }
	    iters++;
	}
    }
    //output data into csv file, hardcoded for two opinions
    ofstream graphStats;
    graphStats.open(fileName);
    graphStats << setiosflags(ios::left) << setiosflags(ios::fixed);
    for(i = 0; i < collectionStep; i++) {
	graphStats << stepTimeCourse(i,0) << ",";
	graphStats << 1.0*minorityOpnTimeCourse(i,0)/n << ",";
	graphStats << 1.0*N10timeCourse(i,0) << "\n";
    }
    graphStats << "\n\n";
    graphStats.close();
    stringstream ss;
    ss << "bifData_" << rewireTo << "_" << n << "_" << avgDeg << ".csv";
    string bifTitle = ss.str();
    ofstream bifData;
    bifData.open(bifTitle, ios::app);
    bifData << a << ",";
    bifData << (opnCounts(0,0)<opnCounts(1,0)?(1.0*opnCounts(0,0)/n):(1.0*opnCounts(1,0)/n)) << ",";
    //flag to determine whether corresponding data represents frozen state or not
    if(conflicts != 0) {
	bifData << "1" << "\n";
    }
    else {
	bifData << "0" << "\n";
    }
    bifData << initDist[0];
    bifData.close();
    return 0;
}

/**
   The following initializes a simple ER random graph based on the algorithm
   presented in "Efficient generation of large random networks" by V. Batagelj
   and U. Brandes. In parallel, each vertex is randomly assigned an opinion according
   to the initial distribution specified by the user.
       
   The resulting graphs have been tested to ensure the algorithm's performance.
**/
void votingModel::initGraph() {
    double p = avgDeg/(n-1);
    int v = 1;
    int w = -1;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 mt1(seed);
    double normalization = (double) mt1.max();
    double rv1, rv2, partialSum;
    int indexCounter;
    //allocate and init degs array to zero
    degs = new int[n];
    for(int i = 0; i < n; i++) {
	degs[i] = 0;
    }
    while(v <= n) {
	rv1 = mt1()/normalization;
	w = (int) (w + 1 + floor(log(1-rv1)/log(1-p)) + ROUND_CONST);
	while((w >= v) && (v <= n)) {
	    rv2 = mt1()/normalization;
	    partialSum = 0.0;
	    indexCounter = 0;
	    while(partialSum < rv2) {
		partialSum = partialSum + initDist[indexCounter];
		indexCounter++;
	    }
	    Opns(v-1,0) = indexCounter;
	    w = w - v;
	    v = v + 1;
	}
	if(v < n) {
	    A(v,w) = 1;
	    A(w,v) = 1;
	    degs[v] = degs[v] + 1;
	    degs[w] = degs[w] + 1;
	}
    }
}
    
//test function, showed that conflicts were correctly tracked throughout the simulation
int votingModel::consistencyCheck() {
    int sum = 0;
    int currentOpn, i, j;
    int conflicts = 0;
    for(i = 0; i < n; i++) {
	sum = 0;
	for(j = 0; j < n; j++) {
	    sum += A(i,j);
	}
	if(sum != degs[i]) {
	    return -1;
	}
    }

    for(i = 0; i < n; i++) {
	currentOpn = Opns(i,0);
	for(j = i + 1; j < n; j++) {
	    if((A(i,j) != 0) && (Opns(j,0) != currentOpn)) {
		conflicts++;
	    }
	}
    }
    return conflicts;
}

