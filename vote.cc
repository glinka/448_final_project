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
#include <time.h>
#include "vote.h"

using namespace Eigen;
using namespace std;
  
votingModel::votingModel(int n, int k, int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, string rewireTo, string fileName): ROUND_CONST(0.01), n(n), A(MatrixXi::Zero(n,n)), Opns(MatrixXi::Zero(n,1)), k(k), maxIter(maxIter), collectionInterval(collectionInterval), a(a), avgDeg(avgDeg), initDist(initDist), rewireTo(rewireTo), fileName(fileName) {};
  
int votingModel::vote() {
    /* Arguments list:
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
<<<<<<< HEAD
    //constant to allow proper casting/rounding of floating point floor function
    clock_t start, beginWrite, end;
    start = clock();
=======
>>>>>>> a23f10bc8d39f929b8a44c3753c455ddb839f202
    int i, j;
    double sum = 0;
    for(i = 0; i < k; i++) {
	sum += initDist[i];
    }
    if(sum != 1) {
	cerr << "sum(u) != 1, improper initial distribution" << "\n";
	return 1;
    }
    initGraph();
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 mt1(seed);
    double normalization = (double) mt1.max();
    //init matrices
    int totalEdges = 0;
    MatrixXi opnCounts = MatrixXi::Zero(k,1);
    int currentOpn;
    for(int i = 0; i < n; i++) {
	currentOpn = Opns(i,0);
	opnCounts(currentOpn-1,0) = opnCounts(currentOpn-1,0) + 1;
	for(j = i+1; j < n; j++) {
	    totalEdges += A(i,j);
	}
    }
    //count initial number of conflicts
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
    MatrixXi N10timeCourse = MatrixXi::Zero((int) maxIter/collectionInterval, 1);
    MatrixXi minorityOpnTimeCourse = MatrixXi::Zero((int) maxIter/collectionInterval, 1);
    MatrixXi stepTimeCourse = MatrixXi::Zero((int) maxIter/collectionInterval, 1);
    MatrixXi V;
<<<<<<< HEAD
    int step = 0;
    int alphaTracker, actionTracker = 0;
=======
    int step, neighborTracker, chosenVertexTracker = 0;
>>>>>>> a23f10bc8d39f929b8a44c3753c455ddb839f202
    int chosenVertex, neighborIndex, neighbor, conflictCounter, newNeighbor;
    double actionToPerform;
    if(rewireTo == "random") {
      while((conflicts > 0) && (iters < maxIter)) {
	  chosenVertex = (int) floor(n*(mt1()/normalization));
	  while(A.block(chosenVertex,0,1,n).sum() == 0) {
	      chosenVertex = (int) floor(n*(mt1()/normalization));
	      chosenVertexTracker++;
	}
	V = MatrixXi::Zero(n,1);
	i = 0;
	for(j = 0; j < n; j++) {
	  if((A(chosenVertex,j) != 0) && (Opns(chosenVertex,0) != Opns(j,0))) {
	    V(i++,0) = j;
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
<<<<<<< HEAD
	if(Opns(chosenVertex,0) != Opns(neighbor,0)) {
	    actionTracker++;
	  actionToPerform = (mt1()/normalization);
=======
	**/
	if(i > 0) {
	  actionToPerform = mt1()/normalization;
>>>>>>> a23f10bc8d39f929b8a44c3753c455ddb839f202
	  conflictCounter = 0;
	  if(actionToPerform > a) {
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
	  else {
	      alphaTracker++;
	      conflictCounter--;
	      A(chosenVertex,neighbor) = 0;
	      A(neighbor,chosenVertex) = 0;
<<<<<<< HEAD
	      newNeighbor = (int) (floor(n*(mt1()/normalization)) + ROUND_CONST);
	      while((A(chosenVertex,newNeighbor) != 0) || (newNeighbor == chosenVertex)) {
		  newNeighbor = (int) floor(n*(mt1()/normalization)) + ROUND_CONST;
=======
	      newNeighbor = (int) floor(n*(mt1()/normalization));
	      /**
	      nNeighbors = A.block(chosenVertex,0,1,n).sum();
	      newNeighborIndex = (int) (floor(nNeighbors*(mt1()/normalization)) + ROUND_CONST);
	      neighborCount = 0;
	      newNeighbor = 0;
	      while(neighborCount < newNeighborIndex) {
		  neighborCount += A(chosenVertex,j);
		  newNeighbor++;
	      }
	      newNeighbor--;
	      **/
	      //for n = 100, avgDeg = 4, have, on avg, 95% chance of success
	      while((A(chosenVertex,newNeighbor) != 0) || (newNeighbor == chosenVertex)) {
		newNeighbor = (int) floor(n*(mt1()/normalization));
		  neighborTracker++;
>>>>>>> a23f10bc8d39f929b8a44c3753c455ddb839f202
	      }
	      A(chosenVertex,newNeighbor) = 1;
	      A(newNeighbor,chosenVertex) = 1;
	      if(Opns(chosenVertex,0) != Opns(newNeighbor,0)) {
		  conflictCounter++;
	      }
	      conflicts += conflictCounter;
	  }
	}
	if(iters % collectionInterval == 0) {
	    step = iters/collectionInterval;
<<<<<<< HEAD
	    minorityOpnTimeCourse(step,0) = opnCounts(0,0)<opnCounts(1,0)?opnCounts(0,0):opnCounts(1,0);
=======
	    minorityOpnTimeCourse(step,0) = opnCounts(0,0);//<opnCounts(1,0)?opnCounts(0,0):opnCounts(1,0);
>>>>>>> a23f10bc8d39f929b8a44c3753c455ddb839f202
	    N10timeCourse(step,0) = conflicts;
	    stepTimeCourse(step,0) = step + 1;
	}
	iters++;
      }
    }
    beginWrite = clock();
    cout << "alpha: " << a << " simulation: " << (1.0*alphaTracker)/actionTracker << endl;
    cout << "time to complete simulation: " << (double) (beginWrite - start)/CLOCKS_PER_SEC << endl;
    ofstream graphStats;
    graphStats.open(fileName);
    graphStats << setiosflags(ios::left) << setiosflags(ios::fixed);
    /**
    graphStats << setw(8) << "Step";
    graphStats << setw(30) << "Number Opn1 Vertice";
    graphStats << setw(30) << "Number Disagreements" << "\n";
    for(i = 0; i < step; i++) {
      graphStats << setw(8) << stepTimeCourse(i,0);
      graphStats << setw(30) <<  minorityOpnTimeCourse(i,0);
      graphStats << setw(30) << N10timeCourse(i,0) << "\n";
    }
    graphStats << "Total Edges: " << totalEdges << "\n";
    **/
    //the above is human readable, but why?
    //below is csv formatting for easy python importing
    for(i = 0; i < step; i++) {
      graphStats << stepTimeCourse(i,0) << ",";
      graphStats << 1.0*minorityOpnTimeCourse(i,0)/n << ",";
      graphStats << 1.0*N10timeCourse(i,0)/totalEdges << "\n";
    }
    graphStats << "\n\n";
    graphStats.close();
    stringstream ss;
    ss << "bifData_" << n << "_" << avgDeg << ".csv";
    string bifTitle = ss.str();
<<<<<<< HEAD
=======
    cout << "new vertex loop count: " << chosenVertexTracker << endl;
    cout << "conflicts: " << conflicts << endl;
>>>>>>> a23f10bc8d39f929b8a44c3753c455ddb839f202
    ofstream bifData;
    bifData.open(bifTitle, ios::app);
    bifData << a << ",";
    bifData << (opnCounts(0,0)<opnCounts(1,0)?(1.0*opnCounts(0,0)/n):(1.0*opnCounts(1,0)/n)) << ",";
    if(iters == maxIter) {
      bifData << "1" << "\n";
    }
    else {
      bifData << "0" << "\n";
    }
    end = clock();
    cout << "write time: " << (double) (end - beginWrite)/CLOCKS_PER_SEC << endl;
    cout << "percent time writing: " << (double) (end - beginWrite)/(end - start) << endl;
    cout << "total time: " << (double) (end - start)/CLOCKS_PER_SEC << endl;
    return 0;
}

void votingModel::initGraph() {
  double p = avgDeg/(n-1);
  int v = 1;
  int w = -1;
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937 mt1(seed);
  double normalization = (double) mt1.max();
  double rv1, rv2, partialSum;
  int indexCounter;
  //    default_random_engine generator;
  //    uniform_real_distribution<double> distr(0.0, 1.0);
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
    }
  }
}
    

int votingModel::countConflicts() {
    int currentOpn, i, j;
    int conflicts = 0;
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

int votingModel::graphConsistencyCheck() {
    //check graph initialization
    double avgDegCheck = 0;
    double opnCheck = 0;
    int i, j;
    for(i = 0; i < n; i++) {
	if(Opns(i,0) == 2) {
	    opnCheck++;
	}
	if((A(i,i) != 0) || (Opns(i,0) <= 0)) {
	    cout << i << endl;
	    return 1;
	}
	for(j = i; j < n; j++) {
	    if((A(i,j) != A(j,i)) || (A(i,j) > 1)) {
		cout << i << j << endl;
		return 1;
	    }
	    avgDegCheck += A(i,j);
	}
    }
    cout << (2.0*avgDegCheck)/n << endl;
    cout << (1.0*opnCheck)/n << endl;
    return 0;
}
