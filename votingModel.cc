#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "votingModel.h"
#include "votingModelCPI.h"


using namespace std;

int votingModel::id_tracker = 0;
  
votingModel::votingModel(int n, int k, long int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, string rewireTo, string fileName):ROUND_CONST(0.01), n(n), k(k), collectionInterval(collectionInterval), maxIter(maxIter), a(a), avgDeg(avgDeg), initDist(initDist), degs(NULL), Opns(NULL), A(NULL), rewireTo(rewireTo), fileName(fileName), id(id_tracker) {
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  mt = new mt19937(seed);
  rnNormalization = (double) (mt->max()+1);
  degs = new int[n];
  Opns = new int[n];
  nConflicts = new int[n];
  opnCounts = new int[k];
  A = new int*[n];
  for(int i = 0; i < n; i++) {
    A[i] = new int[n];
  }
  initGraph(initDist);
  id_tracker++;
};

/**
   initializes and runs voting model to frozen state or maxIter
**/  

double votingModel::genURN() {
  return (*mt)()/rnNormalization;
}

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
  unsigned int i;
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
  vector<int> n10TimeCourse;
  vector<int> minorityOpnTimeCourse;
  vector<int> stepTimeCourse;
  long int iters = 0;
  while(iters < maxIter && conflicts > 0) {
    step();
    if((iters % collectionInterval == 0) || (conflicts==0)) {
      consistencyCheck();
      int minorityOpn = opnCounts[0]<opnCounts[1]?opnCounts[0]:opnCounts[1];
      minorityOpnTimeCourse.push_back(minorityOpn);
      n10TimeCourse.push_back(conflicts);
      stepTimeCourse.push_back(iters+1);
    }
    iters++;
  }
  //output data into csv file, hardcoded for two opinions
  ofstream graphStats;
  graphStats.open(fileName);
  graphStats << setiosflags(ios::left) << setiosflags(ios::fixed);
  graphStats << "n=" << n << ",";
  graphStats << "avgDeg=" << avgDeg << ",";
  graphStats << "alpha=" << a << "\n";
  vector<vector<int> > data;
  vector<int> v;
  data.push_back(stepTimeCourse);
  data.push_back(minorityOpnTimeCourse);
  data.push_back(n10TimeCourse);
  saveData(data, graphStats);
  graphStats.close();

  stringstream ss;
  ss << "bifData_" << rewireTo << "_" << n << "_" << avgDeg << ".csv";
  string bifTitle = ss.str();
  ofstream bifData;
  bifData.open(bifTitle, ios::app);
  vector<double> rowData;
  rowData.push_back(a);
  rowData.push_back(opnCounts[0]<opnCounts[1]?(1.0*opnCounts[0]/n):(1.0*opnCounts[1]/n));
  if(conflicts != 0) {
    rowData.push_back(1);
  }
  else {
    rowData.push_back(0);
  }
  rowData.push_back(initDist[0]);
  data.clear();
  for(i = 0; i < rowData.size(); i++) {
    data.push_back(v);
    data[i].push_back(rowData[i]);
  }
  saveData(data, bifData);
  bifData.close();

  ss.str("");
  ss << "convergenceData_" << a << "_" << n << "_" << initDist[0] << ".csv";
  string cnvTitle = ss.str();
  ofstream convData;
  convData.open(cnvTitle, ios::app);
  rowData.clear();
  rowData.push_back(iters);
  data.clear();
  for(i = 0; i < rowData.size(); i++) {
    data.push_back(v);
    data[i].push_back(rowData[i]);
  }
  saveData(data, convData);
  convData.close();
  return 0;
}

template <typename dataType>
void votingModel::saveData(const vector<dataType> data, ofstream &fileHandle) {
  for(unsigned int i = 0; i < data.size(); i++) {
    fileHandle << data[i] << "\n";
  }
}
template void votingModel::saveData<int>(const vector<int> data, ofstream &fileHandle);

template <typename dataType>
void votingModel::saveData(const vector<vector<dataType> > data, ofstream &fileHandle) {
  //data should be formatted as column vectors, i.e. the total number of different data types
  //to be saved should be data.size();
  int m = data[0].size();
  int n = data.size();
  unsigned int i, j;
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      fileHandle << data[j][i];
      if(j != n-1) {
	fileHandle << ",";
      }
    }
    fileHandle << endl;
  }
}
template void votingModel::saveData<int>(const vector<vector<int> > data, ofstream &fileHandle);

void votingModel::step() {
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
  if(rewireTo == "random") {

      /******************** FROM CONFLICTS ********************/
      //chose edge at random, then select a chosenVertex and a neighbor at random
      //desire chosenEdge in [1, conflicts] instead of [0, conflicts - 1] so add one
      int chosenEdge = ((int) floor(2*conflicts*(genURN()))) + 1;
      int i = 0;
      int edgeCount = 0;
      while(edgeCount < chosenEdge) {
	edgeCount += nConflicts[i++];
      }
      edgeCount -= nConflicts[--i];
      int j = 0;
      while(edgeCount < chosenEdge) {
	if((A[i][j] != 0) && (Opns[i] != Opns[j])) {
	  edgeCount++;
	}
	j++;
      }
      j--;
      int chosenVertex, neighbor;
      if(genURN() > 0.5) {
      	chosenVertex = i;
      	neighbor = j;
      }
      else {
      	chosenVertex = j;
      	neighbor = i;
      }
      double actionToPerform = genURN();
      int conflictCounter = 0;
      //case one, change opinion of chosenVertex to match neighbor
      if(actionToPerform > a) {
	//do nothing if opinions already match
	if(Opns[chosenVertex] != Opns[neighbor]) {
	  //adjust opinions and opinion counts
	  opnCounts[Opns[chosenVertex] - 1]--;
	  opnCounts[Opns[neighbor] - 1]++;
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
	degs[neighbor]--;
	int newNeighbor = (int) floor(n*(genURN()));
	while((A[chosenVertex][newNeighbor] != 0) || (newNeighbor == chosenVertex)) {
	  newNeighbor = (int) floor(n*(genURN()));
	}
	A[chosenVertex][newNeighbor] = 1;
	A[newNeighbor][chosenVertex] = 1;
	degs[newNeighbor]++;
	if(Opns[chosenVertex] != Opns[newNeighbor]) {
	  conflictCounter++;
	  nConflicts[chosenVertex]++;
	  nConflicts[newNeighbor]++;
	}
	conflicts += conflictCounter;
      }
  }
  /**

*********************************************************************************
Rewire-to-same continues to choose vertex at random instead of edge
*********************************************************************************

**/
  else if(rewireTo == "same") {
      //find vertex with nonzero degree
      int chosenVertex = (int) floor(n*(genURN()));
      int deg = 0;
      int i, j;
      for(j = 0; j < n; j++) {
	deg += A[chosenVertex][j];
      }
      while(deg == 0) {
	chosenVertex = (int) floor(n*(genURN()));
	deg = 0;
	for(j = 0; j < n; j++) {
	  deg += A[chosenVertex][j];
	}
      }
      //assemble list of possible neighbors, storing their indices in V
      int V[n];
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
      int neighborIndex = (int) floor((i)*(genURN()));
      int neighbor = V[neighborIndex];
      //check that chosenVertex has a disagreeing neighbor, else do nothing
      if(i > 0) {
	double actionToPerform = genURN();
	int conflictCounter = 0;
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
	  int newNeighbor = (int) floor(n*(genURN()));
	  while((Opns[chosenVertex] != Opns[newNeighbor]) || (A[chosenVertex][newNeighbor] != 0) || (newNeighbor == chosenVertex)) {
	    newNeighbor = (int) floor(n*(genURN()));
	  }
	  A[chosenVertex][newNeighbor] = 1;
	  A[newNeighbor][chosenVertex] = 1;
	  conflicts += conflictCounter;
	}
      }
  }
}


/**
   The following initializes a simple ER random graph based on the algorithm
   presented in "Efficient generation of large random networks" by V. Batagelj
   and U. Brandes. In parallel, each vertex is randomly assigned an opinion according
   to the initial distribution specified by the user.
       
   The resulting graphs have been tested to ensure the algorithm's performance.
**/
void votingModel::initGraph(double *dist) {
    double p = avgDeg/(n-1);
    int v = 1;
    int w = -1;
    double rv1, rv2, partialSum;
    int indexCounter, i, j;
    //allocate and init arrays to zero
    for(i = 0; i < n; i++) {
      degs[i] = 0;
      Opns[i] = 0;
      nConflicts[i] = 0;
      for(j = 0; j < n; j++) {
	A[i][j] = 0;
      }
    }
    for(i = 0; i < k; i++) {
      opnCounts[i] = 0;
    }
    while(v <= n) {
      rv1 = genURN();
      w = (int) (w + 1 + floor(log(1-rv1)/log(1-p)) + ROUND_CONST);
      while((w >= v) && (v <= n)) {
	rv2 = genURN();
	partialSum = 0.0;
	indexCounter = 0;
	while(partialSum < rv2) {
	  partialSum = partialSum + dist[indexCounter];
	  indexCounter++;
	}
	Opns[v-1] = indexCounter;
	opnCounts[indexCounter - 1]++;
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
    conflicts = 0;
    int currentOpn;
    for(i = 0; i < n; i++) {
      currentOpn = Opns[i];
      for(j = i+1; j < n; j++) {
	if(A[i][j] > 0 && currentOpn != Opns[j]) {
	  conflicts++;
	  nConflicts[i]++;
	  nConflicts[j]++;
	}
      }
    }
}

void votingModel::initGraph(double *dist, int conflictCount) {
    int i, j;
    for(i = 0; i < n; i++) {
      degs[i] = 0;
      Opns[i] = 0;
      nConflicts[i] = 0;
      for(j = 0; j < n; j++) {
	A[i][j] = 0;
      }
    }
    for(i = 0; i < k; i++) {
      opnCounts[i] = 0;
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
    double pConflict = 1.0*conflictCount/nEdges;
    int edgeCount;
    int chosenVertex, neighbor;
    for(edgeCount = 0; edgeCount < nEdges; edgeCount++) {
      chosenVertex = n*genURN();
      if(chosenVertex > nOnes-1) {
	//operate on twos
	if(genURN() >= pConflict) {
	  //no conflict, two to two edge
	  neighbor = nTwos*genURN();
	  while(chosenVertex == nOnes-1+neighbor) {
	    neighbor = nTwos*genURN();
	  }	    
	  A[chosenVertex][nOnes-1+neighbor] = 1;
	  A[nOnes-1+neighbor][chosenVertex] = 1;
	}
	else {
	  //conflict, two to one edge
	  neighbor = nOnes*genURN();
	  while(chosenVertex == neighbor) {
	    neighbor = nOnes*genURN();
	  }	    
	  A[chosenVertex][neighbor] = 1;
	  A[neighbor][chosenVertex] = 1;
	}
      }
      else {
	if(genURN() > pConflict) {
	  //no conflict, one to one edge
	  neighbor = nOnes*genURN();
	  while(chosenVertex == neighbor) {
	    neighbor = nOnes*genURN();
	  }	    
	  A[chosenVertex][neighbor] = 1;
	  A[neighbor][chosenVertex] = 1;
	}
	else {
	  //conflict, one to two edge
	  neighbor = nTwos*genURN();
	  while(chosenVertex == nOnes-1+neighbor) {
	    neighbor = nTwos*genURN();
	  }	    
	  A[chosenVertex][nOnes-1+neighbor] = 1;
	  A[nOnes-1+neighbor][chosenVertex] = 1;
	}
      }
    }
    conflicts = 0;
    for(i = 0; i < n; i++) {
      opnCounts[Opns[i] - 1]++;
      for(j = 0; j < n; j++) {
	degs[i] += A[i][j];
	if(A[i][j] > 0 && Opns[i] != Opns[j]) {
	  nConflicts[i]++;
	  conflicts++;
	}
      }
    }
    conflicts /= 2;
}
    
//test function, showed that conflicts and degrees were correctly tracked throughout the simulation
int votingModel::consistencyCheck() {
  int sum = 0;
  int currentOpn, i, j;
  int conflictCount = 0;
  for(i = 0; i < n; i++) {
    if(A[i][i] != 0) {
      return -2;
    }
    sum = 0;
    for(j = 0; j < n; j++) {
      sum += A[i][j];
    }
    if(sum != degs[i]) {
      return -1;
    }
  }
  for(i = 0; i < n; i++) {
    int nconflict_check = 0;
    currentOpn = Opns[i];
    for(j = 0; j < n; j++) {
      if((A[i][j] != 0) && (Opns[j] != currentOpn)) {
	nconflict_check++;
	conflictCount++;
      }
    }
    if(nconflict_check != nConflicts[i]) {
      return -3;
    }
  }
  return conflictCount/2;
}
