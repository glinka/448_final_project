#include <stdlib.h>
#include <string>
#include <sstream>
#include "votingModel.h"
#include "votingModelCPI.h"

using namespace std;

int main(int argc, char *argv[]) {
  double avgDeg = 4;
  int i, j;
  bool project = false;
  string rewireTo = "random";
  int n = 200;
  long int maxIter = 5*n*n;
  int collectionInterval;
  double a = 0.5;
  int k = 2;
  double *initDist = new double[2];
  int nVMS = 4;
  double projectionStep;
  initDist[0] = 0.5;
  initDist[1] = 0.5;
  //loop through all arguments, assign variables as needed
  for(i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {
      string currentLabel = argv[i];
      char *currentArg = argv[i+1];
      if(currentLabel == "-n" || currentLabel == "-nNodes" || currentLabel == "-nodes") {
	n = atoi(currentArg);
      }
      else if(currentLabel == "-maxIter" || currentLabel == "-mi" || currentLabel == "-iter") {
	maxIter = atoi(currentArg);
      }
      else if(currentLabel == "-collectionInterval" || currentLabel == "-ci" || currentLabel == "-colliter") {
	collectionInterval = atoi(currentArg);
      }
      else if(currentLabel == "-a" || currentLabel == "-alpha") {
	a = atof(currentArg);
      }
      else if(currentLabel == "-k" || currentLabel == "-nOpinions" || currentLabel == "-opinions") {
	k = atoi(currentArg);
	delete[] initDist;
	initDist = new double[k];
	int distOffset = 2;
	for(j = distOffset; j < k+distOffset; j++) {
	  initDist[j-distOffset] = atof(argv[i+j]);
	}
      }
      else if(currentLabel == "-avgDeg" || currentLabel == "-deg" || currentLabel == "-ad") {
	avgDeg = atof(currentArg);
      }
      else if(currentLabel == "-project" || currentLabel == "-toProject") {
	//projection false by default
	if(currentArg[0] == 't' || currentArg[0] == 'T' || currentArg[0] == '1') {
	  project = true;
	}
      }
      else if(currentLabel == "-projectionStep" || currentLabel == "-ps" || currentLabel == "-pi") {
	project = true;
	projectionStep = atof(currentArg);
      }
      else if(currentLabel == "-numberVotingModels" || currentLabel == "-numberVMS" || currentLabel == "-nVMS" || currentLabel == "-nvms") {
	project = true;
	nVMS = atoi(currentArg);
      }
    }
  }
  int nMS = n;
  int waitingPeriod = 100;
  projectionStep = n;
  collectionInterval = nMS/10;
  stringstream ss;
  ss << "graphStats_" << "n_" << n << "_avgDeg_" << avgDeg << "_alpha_" << a << "_projectionStep_" << projectionStep << "_initDist";
  for(i = 0; i < k; i++) {
    ss << "_" << initDist[i];
  }
  ss << ".csv";
  string fileName = ss.str();
  vector<votingModel> vmV;
  for(i = 0; i < nVMS; i++) {
    vmV.push_back(votingModel(n, k, maxIter, collectionInterval, a, avgDeg, initDist, rewireTo, fileName, project));
  }
  votingModelCPI *cpi = new votingModelCPI(vmV, projectionStep, waitingPeriod, collectionInterval, nMS);
  cpi->run(maxIter);
  delete cpi;
  delete[] initDist;
  return 0;
}

