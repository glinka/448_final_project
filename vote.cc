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
  double a = 0.5;
  int k = 2;
  double *initDist = new double[k];
  int nVMS = 4;
  double projectionStep = n;
  int save_interval = 1000;
  int nMS = n;
  int collectionInterval = nMS/10;
  initDist[0] = 0.5;
  initDist[1] = 0.5;
  int nruns = 64;
  int waitingPeriod = 100;
  bool alpha_range = false;
  bool init_range = false;
  bool alter = false;
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
	cout << "collection interval: " << collectionInterval << endl;
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
      else if(currentLabel == "-saveInterval" || currentLabel == "-saveinterval" || currentLabel == "-si") {
	project = true;
	save_interval = atoi(currentArg);
      }
      else if(currentLabel == "-nMicroSteps" || currentLabel == "-microSteps" || currentLabel == "-nms" || currentLabel == "-ms") {
	project = true;
	nMS = atoi(currentArg);
      }
      else if(currentLabel == "-waitingperiod" || currentLabel == "-wp") {
	project = true;
	waitingPeriod = atoi(currentArg);
      }
      else if(currentLabel == "-nruns") {
	project = false;
	nruns = atoi(currentArg);
      }
      else if(currentLabel == "-arange" || currentLabel == "-alpharange") {
	alpha_range = true;
      }
      else if(currentLabel == "-irange" || currentLabel == "-initrange") {
	init_range = true;
      }
      else if(currentLabel == "-alter") {
	alter = true;
      }
      else {
	cout << currentLabel << " -- flag not recognized" << endl;
      }
    }
  }
  stringstream ss;
  ss << "n=" << n;
  ss << ",nVms=" << nVMS;
  ss << ",nsteps=" << maxIter;
  ss << ",avg_deg=" << avgDeg;
  ss << ",alpha=" << a;
  /**
     
     left out for ease of parsing arguments in python

  ss << ",init_ndist=";
  for(i = 0; i < k; i++) {
    ss << initDist[i] << ",";
  }
  //leave out last comma due to previous loop
  ss << "rewire_to=" << rewireTo;
  **/
  string file_header = ss.str();
  ss.str("");
  ss << "_nvms_" << nVMS;
  ss << "n_" << n;
  ss << "_nsteps_" << maxIter;
  ss << "_nms_" << nMS;
  ss << "_projstep_" << projectionStep;
  string file_name = ss.str();
  if(project) {
    vector<votingModel> vmV;
    for(i = 0; i < nVMS; i++) {
      vmV.push_back(votingModel(n, k, maxIter, collectionInterval, a, avgDeg, initDist, rewireTo, ""));
    }
    votingModelCPI *cpi = new votingModelCPI(vmV, waitingPeriod, collectionInterval, nMS, file_header, file_name);
    //TEST


    vector< vector<double> > testdata(1);
    cpi->easy_average<double>(testdata);


    //ENDTEST
    cpi->run(maxIter, projectionStep, save_interval);
    delete cpi;
  }
  else {
    for(i = 0; i < nruns; i++) {
      if(alpha_range) {
	  if(init_range) {
	    // range through both alpha and init fracs
	    for(double alpha = 0.05; alpha < 1; alpha+=0.1) {
	      for(double initfrac = 0.05; initfrac < 1; initfrac+=0.1) {
		ss.str("");
		ss << "single_runs/graphstats_n_" << n;
		ss << "_alpha_" << alpha;
		ss << "_nsteps_" << maxIter;
		ss << "_initfrac_" << initfrac;
		ss << "_rewireto_" << rewireTo;
		ss << "_" << i;
		file_name = ss.str();
		double initdist[2] = {initfrac, 1-initfrac};
		votingModel vm(n, k, maxIter, collectionInterval, alpha, avgDeg, initdist, rewireTo, file_name);
		vm.vote(alter, "");
	      }		
	    }
	  }
	  else {
	    // range through alpha
	    for(double alpha = 0.05; alpha < 1; alpha+=0.1) {
		ss.str("");
		ss << "single_runs/graphstats_n_" << n;
		ss << "_alpha_" << alpha;
		ss << "_nsteps_" << maxIter;
		ss << "_initfrac_" << initDist[0];
		ss << "_rewireto_" << rewireTo;
		ss << "_" << i;
		file_name = ss.str();
		votingModel vm(n, k, maxIter, collectionInterval, alpha, avgDeg, initDist, rewireTo, file_name);
		vm.vote(alter, "");
	    }
	  }
      }
      else {
	// range through init fracs
	if(init_range) {
	      for(double initfrac = 0.05; initfrac < 1; initfrac+=0.1) {
		ss.str("");
		ss << "single_runs/graphstats_n_" << n;
		ss << "_alpha_" << a;
		ss << "_nsteps_" << maxIter;
		ss << "_initfrac_" << initfrac;
		ss << "_rewireto_" << rewireTo;
		ss << "_" << i;
		file_name = ss.str();
		double initdist[2] = {initfrac, 1-initfrac};
		votingModel vm(n, k, maxIter, collectionInterval, a, avgDeg, initdist, rewireTo, file_name);
		vm.vote(alter, "");
	      }		
	}
	else if(alter) {
	  vector<bool> alters(2);
	  alters[0] = true;
	  alters[1] = false;
	  for(int ii = 4; ii < 5; ii++) {
	    for(int j = 0; j < 1; j++) {
	      initDist[0] = 0.1*ii;
	      initDist[1] = 1 - initDist[0];
	      votingModel vm(n, k, maxIter, collectionInterval, a, avgDeg, initDist, rewireTo, file_name);
	      ss.str("");
	      ss << ii;
	      // cout << id << endl;
	      vm.vote(alters[j], ss.str());
	    }
	  }
	}
	// don't range through any vars
	else {
	  ss.str("");
	  ss << "single_runs/graphstats_n_" << n;
	  ss << "_alpha_" << a;
	  ss << "_nsteps_" << maxIter;
	  ss << "_initfrac_" << initDist[0];
	  ss << "_rewireto_" << rewireTo;
	  ss << "_" << i;
	  file_name = ss.str();
	  votingModel vm(n, k, maxIter, collectionInterval, a, avgDeg, initDist, rewireTo, file_name);
	  vm.vote(alter, "");
	}
      }
    }
  }
  delete[] initDist;
  return 0;
}
