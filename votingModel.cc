#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <iomanip>
#include "votingModel.h"
#include "votingModelCPI.h"


using namespace std;
  
votingModel::votingModel(int n, int k, int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, string rewireTo, string fileName, bool project, votingModelCPI *CPI):project(project), ROUND_CONST(0.01), n(n), k(k), maxIter(maxIter), collectionInterval(collectionInterval), a(a), avgDeg(avgDeg), initDist(initDist), degs(NULL), Opns(NULL), A(NULL), rewireTo(rewireTo), fileName(fileName), vmCPI(CPI) {
};

/**
   initializes and runs voting model to frozen state or maxIter
**/  
int votingModel::vote() {
  int waitingPeriod = 10000;
  int projectionInterval = 1000;
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
  initGraph(initDist, 0, true);
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937 mt1(seed);
  //adding one to ensure the distribution lies in [0,1)
  double normalization = (double) (mt1.max()+1);
  //init matrices
  int totalEdges = 0;
  int currentOpn;
  int opnCounts[k];
  for(i = 0; i < k; i++) {
    opnCounts[i] = 0;
  }
  //counts initial numbers of each opinion
  for(i = 0; i < n; i++) {
    currentOpn = Opns[i];
    opnCounts[currentOpn-1] = opnCounts[currentOpn-1] + 1;
    for(j = i+1; j < n; j++) {
      totalEdges += A[i][j];
    }
  }
  //count initial number of conflicts
  /**
     Initially the code calculated the number of conflicts from scratch at each step (i.e from
     the adjacency matrix and opinion vector), but this was prohibitively expensive
  **/
  int conflicts = 0;
  nConflicts = new int[n];
  for(i = 0; i < n; i++) {
    nConflicts[i] = 0;
  }
  for(i = 0; i < n; i++) {
    currentOpn = Opns[i];
    for(j = i+1; j < n; j++) {
      if((A[i][j] != 0) && (Opns[j] != currentOpn)) {
	conflicts++;
	nConflicts[i]++;
	nConflicts[j]++;
      }
    }
  }
  sum = 0;
  for(i = 0; i < n; i++) {
    sum += nConflicts[i];
  }
  if(sum != 2*conflicts) {
    cout << "nConflicts not correctly calculated" << endl;
  }
  /**
     iters: number of steps the simulation has executed
     N10timeCourse: stores the number of discordant edges
     minorityOpnTimeCourse: stores the number of vertices holding the minority opinion
     V: during each step, stores the index of potential neighbors to chosenVertex (only used in rewire-to-same)
     ***********
     chosenVertex: in rewire-to-random, randomly chosen vertex of randomly chosen edge, in
     rewire-to-same, randomly chosen vertex from graph
     ***********
     neighborIndex: temporary array index of the chosenVertex's neighbor in V
     neighbor: index of chosenVertex's neighbor
     conflictCounter: tracks changes number of conflicts in system at each step
     newNeighbor: index of new neighbor of chosenVertex (if rewiring)
     collectionStep: number of times data has been collected
     edgeCount: counter used to find randomly chosen edge
     actionToPerform: uniform random number on [0,1), determines what to do with the edge
  **/
  int iters = 0;
  int chosenVertex, chosenEdge, neighborIndex, neighbor, conflictCounter, newNeighbor, collectionStep, edgeCount;
  double actionToPerform;
  int nData = ((int) maxIter/collectionInterval);
  int N10timeCourse[nData];
  int minorityOpnTimeCourse[nData];
  int stepTimeCourse[nData];
  int waitCounter = 0;
  for(i = 0; i < nData; i++) {
    N10timeCourse[i] = 0;
    minorityOpnTimeCourse[i] = 0;
    stepTimeCourse[i] = 0;
  }
  int V[n];
  if(rewireTo == "random") {
    while((conflicts > 0) && (iters < maxIter)) {
      /******************** FROM CONFLICTS ********************/
      //chose edge at random, then select a chosenVertex and a neighbor at random
      //desire chosenEdge in [1, conflicts] instead of [0, conflicts - 1] so add one
      chosenEdge = ((int) floor(2*conflicts*(mt1()/normalization))) + 1;
      i = 0;
      edgeCount = 0;
      while(edgeCount < chosenEdge) {
	edgeCount += nConflicts[i++];
      }
      edgeCount -= nConflicts[--i];
      j = 0;
      while(edgeCount < chosenEdge) {
	if((A[i][j] != 0) && (Opns[i] != Opns[j])) {
	  edgeCount++;
	}
	j++;
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
      // chosenEdge = ((int) floor(2*totalEdges*(mt1()/normalization))) + 1;
      // i = 0;
      // edgeCount = 0;
      // while(edgeCount < chosenEdge) {
      // 	edgeCount += degs[i++];
      // }
      // edgeCount -= degs[--i];
      // j = 0;
      // while(edgeCount < chosenEdge) {
      // 	edgeCount += A[i][j++];
      // }
      // j--;
      // if(mt1()/normalization > 0.5) {
      // 	chosenVertex = i;
      // 	neighbor = j;
      // }
      // else {
      // 	chosenVertex = j;
      // 	neighbor = i;
      // }
      actionToPerform = mt1()/normalization;
      conflictCounter = 0;
      //case one, change opinion of chosenVertex to match neighbor
      if(actionToPerform > a) {
	//do nothing if opinions already match
	if(Opns[chosenVertex] != Opns[neighbor]) {
	  //adjust opinions and opinion counts
	  opnCounts[Opns[chosenVertex] - 1] = opnCounts[Opns[chosenVertex] - 1] - 1;
	  opnCounts[Opns[neighbor] - 1] = opnCounts[Opns[neighbor] - 1] + 1;
	  Opns[chosenVertex] = Opns[neighbor];
	  //adjust conflicts

	  /*************** need to adjust for a more general k-opinion model ***************/

	  for(j = 0; j < n; j++) {
	    if((A[chosenVertex][j] != 0) && (Opns[j] != Opns[chosenVertex])) {
	      conflictCounter++;
	      nConflicts[chosenVertex]++;
	      nConflicts[j]++;
	    }
	    else if(A[chosenVertex][j] != 0) {
	      conflictCounter--;
	      nConflicts[chosenVertex]--;
	      nConflicts[j]--;
	    }
	  }
	  conflicts += conflictCounter;
	}
      }
      //case two, rewire edge
      else {
	if(Opns[chosenVertex] != Opns[neighbor]) {
	  conflictCounter--;
	  nConflicts[chosenVertex]--;
	  nConflicts[neighbor]--;
	}
	//disconnect edge and find new neighbor from all vertices, ensuring no loops or
	//parallel edges are formed
	A[chosenVertex][neighbor] = 0;
	A[neighbor][chosenVertex] = 0;
	degs[neighbor] = degs[neighbor] - 1;
	newNeighbor = (int) floor(n*(mt1()/normalization));
	while((A[chosenVertex][newNeighbor] != 0) || (newNeighbor == chosenVertex)) {
	  newNeighbor = (int) floor(n*(mt1()/normalization));
	}
	A[chosenVertex][newNeighbor] = 1;
	A[newNeighbor][chosenVertex] = 1;
	degs[newNeighbor] = degs[newNeighbor] + 1;
	if(Opns[chosenVertex] != Opns[newNeighbor]) {
	  conflictCounter++;
	  nConflicts[chosenVertex]++;
	  nConflicts[newNeighbor]++;
	}
	conflicts += conflictCounter;
      }
	
	
      //collect data every collectionInterval steps
      if(project && (waitCounter > waitingPeriod) && (conflicts > 0)) {
	if(iters % (projectionInterval/100) == 0) {
	  vector<double> data;
	  data.push_back(opnCounts[0]<opnCounts[1]?1.0*opnCounts[0]/n:1.0*opnCounts[1]/n);
	  data.push_back(conflicts);
	  data.push_back(iters);
	  vmCPI->collectData(data);
	}
	if(iters % projectionInterval == 0 && iters > 0) {
	  waitCounter = 0;
	  vector<double> newPts = vmCPI->project();
	  double newDist[2] = {newPts[0], 1-newPts[0]};
	  initGraph(newDist, newPts[1], false);
	  //reset previously calculated properties
	  totalEdges = 0;
	  for(i = 0; i < k; i++) {
	    opnCounts[i] = 0;
	  }
	  //counts initial numbers of each opinion
	  for(i = 0; i < n; i++) {
	    currentOpn = Opns[i];
	    opnCounts[currentOpn-1] = opnCounts[currentOpn-1] + 1;
	    for(j = i+1; j < n; j++) {
	      totalEdges += A[i][j];
	    }
	  }
	  //count initial number of conflicts
	  /**
	     Initially the code calculated the number of conflicts from scratch at each step (i.e from
	     the adjacency matrix and opinion vector), but this was prohibitively expensive
	  **/
	  conflicts = 0;
	  for(i = 0; i < n; i++) {
	    nConflicts[i] = 0;
	  }
	  for(i = 0; i < n; i++) {
	    currentOpn = Opns[i];
	    for(j = i+1; j < n; j++) {
	      if((A[i][j] != 0) && (Opns[j] != currentOpn)) {
		conflicts++;
		nConflicts[i]++;
		nConflicts[j]++;
	      }
	    }
	  }
	}
      }

      if((iters % collectionInterval == 0) || (conflicts==0)) {
	collectionStep = iters/collectionInterval;
	minorityOpnTimeCourse[collectionStep] = opnCounts[0]<opnCounts[1]?opnCounts[0]:opnCounts[1];
	N10timeCourse[collectionStep] = conflicts;
	stepTimeCourse[collectionStep] = collectionStep + 1;
      }
      iters++;
      waitCounter++;
    }
  }

  /**

*********************************************************************************
Rewire-to-same continues to choose vertex at random instead of edge
*********************************************************************************

**/
  else if(rewireTo == "same") {
    while((conflicts > 0) && (iters < maxIter)) {
      //find vertex with nonzero degree
      chosenVertex = (int) floor(n*(mt1()/normalization));
      int deg = 0;
      for(j = 0; j < n; j++) {
	deg += A[chosenVertex][j];
      }
      while(deg == 0) {
	chosenVertex = (int) floor(n*(mt1()/normalization));
	deg = 0;
	for(j = 0; j < n; j++) {
	  deg += A[chosenVertex][j];
	}
      }
      //assemble list of possible neighbors, storing their indices in V
      for(i = 0; i < n; i++) {
	V[i] = 0;
      }
      i = 0;
      for(j = 0; j < n; j++) {
	if((A[chosenVertex][j] != 0) && (Opns[chosenVertex] != Opns[j])) {
	  V[i] = j;
	  i++;
	}
      }
      neighborIndex = (int) floor((i)*(mt1()/normalization));
      neighbor = V[neighborIndex];
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
	  opnCounts[Opns[chosenVertex] - 1] = opnCounts[Opns[chosenVertex] - 1] - 1;
	  opnCounts[Opns[neighbor] - 1] = opnCounts[Opns[neighbor] - 1] + 1;
	  Opns[chosenVertex] = Opns[neighbor];
	  for(j = 0; j < n; j++) {
	    if((A[chosenVertex][j] != 0) && (Opns[j] != Opns[chosenVertex])) {
	      conflictCounter++;
	    }
	    else if(A[chosenVertex][j] != 0) {
	      conflictCounter--;
	    }
	  }
	  conflicts += conflictCounter;
	}
	//case two, rewire edge
	else {
	  //remove current edge and find new edge from those vertices that share chosenVertex's opini
	  conflictCounter--;
	  A[chosenVertex][neighbor] = 0;
	  A[neighbor][chosenVertex] = 0;
	  newNeighbor = (int) floor(n*(mt1()/normalization));
	  while((Opns[chosenVertex] != Opns[newNeighbor]) || (A[chosenVertex][newNeighbor] != 0) || (newNeighbor == chosenVertex)) {
	    newNeighbor = (int) floor(n*(mt1()/normalization));
	  }
	  A[chosenVertex][newNeighbor] = 1;
	  A[newNeighbor][chosenVertex] = 1;
	  conflicts += conflictCounter;
	}
      }
      if((iters % collectionInterval == 0) || (conflicts==0)) {
	collectionStep = iters/collectionInterval;
	minorityOpnTimeCourse[collectionStep] = opnCounts[0];//<opnCounts[1]?opnCounts[0]:opnCounts[1];
	N10timeCourse[collectionStep] = conflicts;
	stepTimeCourse[collectionStep] = collectionStep + 1;
      }
      iters++;
    }
  }
  //output data into csv file, hardcoded for two opinions
  ofstream graphStats;
  graphStats.open(fileName);
  graphStats << setiosflags(ios::left) << setiosflags(ios::fixed);
  for(i = 0; i < collectionStep; i++) {
    graphStats << stepTimeCourse[i] << ",";
    graphStats << 1.0*minorityOpnTimeCourse[i]/n << ",";
    graphStats << 1.0*N10timeCourse[i] << "\n";
  }
  graphStats << "\n\n";
  graphStats.close();
  stringstream ss;
  ss << "bifData_" << rewireTo << "_" << n << "_" << avgDeg << ".csv";
  string bifTitle = ss.str();
  ofstream bifData;
  bifData.open(bifTitle, ios::app);
  bifData << a << ",";
  bifData << (opnCounts[0]<opnCounts[1]?(1.0*opnCounts[0]/n):(1.0*opnCounts[1]/n)) << ",";
  //flag to determine whether corresponding data represents frozen state or not
  if(conflicts != 0) {
    bifData << "1" << ",";
  }
  else {
    bifData << "0" << ",";
  }
  bifData << initDist[0] << "\n";
  bifData.close();
  ss.str("");
  ss << "convergenceData_" << a << "_" << n << "_" << initDist[0] << ".csv";
  string cnvTitle = ss.str();
  ofstream convData;
  convData.open(cnvTitle, ios::app);
  convData << iters << "," << vmCPI->getProjectionStep();
  convData << "\n";
  convData.close();
  return 0;
}

/**
   The following initializes a simple ER random graph based on the algorithm
   presented in "Efficient generation of large random networks" by V. Batagelj
   and U. Brandes. In parallel, each vertex is randomly assigned an opinion according
   to the initial distribution specified by the user.
       
   The resulting graphs have been tested to ensure the algorithm's performance.
**/
void votingModel::initGraph(double *dist, int conflicts, bool firstInitialization) {
  if(firstInitialization) {
    double p = avgDeg/(n-1);
    int v = 1;
    int w = -1;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 mt1(seed);
    double normalization = (double) mt1.max();
    double rv1, rv2, partialSum;
    int indexCounter, i, j;
    //allocate and init arrays to zero
    degs = new int[n];
    Opns = new int[n];
    A = new int*[n];
    for(i = 0; i < n; i++) {
      degs[i] = 0;
      Opns[i] = 0;
      A[i] = new int[n];
      for(j = 0; j < n; j++) {
	A[i][j] = 0;
      }
    }
    while(v <= n) {
      rv1 = mt1()/normalization;
      w = (int) (w + 1 + floor(log(1-rv1)/log(1-p)) + ROUND_CONST);
      while((w >= v) && (v <= n)) {
	rv2 = mt1()/normalization;
	partialSum = 0.0;
	indexCounter = 0;
	while(partialSum < rv2) {
	  partialSum = partialSum + dist[indexCounter];
	  indexCounter++;
	}
	Opns[v-1] = indexCounter;
	w = w - v;
	v = v + 1;
      }
      if(v < n) {
	A[v][w] = 1;
	A[w][v] = 1;
	degs[v] = degs[v] + 1;
	degs[w] = degs[w] + 1;
      }
    }
  }
  else {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 mt1(seed);
    double normalization = (double) mt1.max();
    int i, j;
    for(i = 0; i < n; i++) {
      degs[i] = 0;
      Opns[i] = 0;
      for(j = 0; j < n; j++) {
	A[i][j] = 0;
      }
    }
    int nEdges = avgDeg*n/2;
    int nOnes = (int) (dist[0]*n);
    int nTwos = n - nOnes;
    for(i = 0; i < n; i++) {
      if(i < nOnes) {
	Opns[i] = 1;
      }
      else {
	Opns[i] = 2;
      }
    }
    double pConflict = 1.0*conflicts/nEdges;
    int edgeCount;
    int chosenVertex, neighbor;
    for(edgeCount = 0; edgeCount < nEdges; edgeCount++) {
      chosenVertex = n*mt1()/normalization;
      if(chosenVertex > nOnes-1) {
	//operate on twos
	if(mt1()/normalization > pConflict) {
	  //no conflict, two to two edge
	  neighbor = nTwos*mt1()/normalization;
	  A[chosenVertex][nOnes-1+neighbor] = 1;
	  A[nOnes-1+neighbor][chosenVertex] = 1;
	}
	else {
	  //conflict, two to one edge
	  neighbor = nOnes*mt1()/normalization;
	  A[chosenVertex][neighbor] = 1;
	  A[neighbor][chosenVertex] = 1;
	}
      }
      else {
	if(mt1()/normalization > pConflict) {
	  //no conflict, one to one edge
	  neighbor = nOnes*mt1()/normalization;
	  A[chosenVertex][neighbor] = 1;
	  A[neighbor][chosenVertex] = 1;
	}
	else {
	  //conflict, one to two edge
	  neighbor = nTwos*mt1()/normalization;
	  A[chosenVertex][nOnes-1+neighbor] = 1;
	  A[nOnes-1+neighbor][chosenVertex] = 1;
	}
      }
    }
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
	degs[i] += A[i][j];
      }
    }
  }
}
    
//test function, showed that conflicts and degrees were correctly tracked throughout the simulation
int votingModel::consistencyCheck() {
  int sum = 0;
  int currentOpn, i, j;
  int conflicts = 0;
  for(i = 0; i < n; i++) {
    sum = 0;
    for(j = 0; j < n; j++) {
      sum += A[i][j];
    }
    if(sum != degs[i]) {
      return -1;
    }
  }

  for(i = 0; i < n; i++) {
    currentOpn = Opns[i];
    for(j = i + 1; j < n; j++) {
      if((A[i][j] != 0) && (Opns[j] != currentOpn)) {
	conflicts++;
      }
    }
  }
  return conflicts;
}
