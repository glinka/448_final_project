#include <string>
#include <sstream>
#include "fitCurves.h"
#include "votingModel.h"
#include "votingModelCPI.h"

using namespace std;

typedef vector<votingModel>::iterator vmIt

votingModelCPI::votingModelCPI(vector<votingModel> vms, const int projStep, int waitingPeriod, int collectionInterval): vms(vms), projStep(projStep), waitingPeriod(waitingPeriod), collectionInterval(collectionInterval), nMicroSteps(nMicroSteps) {
};

votingModelCPI::run(int nSteps) {
  int step = 0;
  int microStepCount = 0;
  int i, j;
  int nVms = vms.size();
  //initialize space for each vm's time-courses, all will share the same time vector
  stringstream ss;
  ss << "CPIAdj_nSteps_" << nSteps << "_projStep_" << projStep << ".csv";
  string adjFilename = ss.str();
  ofstream cpiAdjData;
  cpiAdjData.open(adjFilename);
  ss.str("");
  ss << "CPIOpns_nSteps_" << nSteps << "_projStep_" << projStep << ".csv";
  string opnsFilename = ss.str();
  ofstream cpiOpnsData;
  cpiOpnsData.open(opnsFilename);
  bool addHeader = true;
  vector<matrix> tempM;
  vector<vect> tempV;
  while(step < nSteps) {
    for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
      vm->step();
    }
    step++;
    microStepCount++;
    if(microStepCount > waitingPeriod) {
      int onManifoldStep = microStepCount - waitingPeriod;
      if(onManifoldStep % collectionInterval == 0) {
	adjMatrices.push_back(tempM);
	opns.push_back(tempV);
	for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
	  adjMatrices.back().push_back(vm->getAdjMatrix());
	  opns.back().push_back(vm->getOpns());
	}
	times.push_back(onManifoldStep);
	if(addHeader) {
	  int n = adjMatrices[0][0].size();
	  cpiAdjData << "n=" << n << ",nVms=" << nVms << endl;
	  cpiOpnsData << "n=" << n << ",nVms=" << nVms << endl;
	  addHeader = false;
	}
      }
      if(microStepCount == nMicroSteps) {
	//just project the minority fraction you fool
	vmVects opnsTC = average(opns);
	vmMatrices adjTC = average(adjMatrices);
	saveData(opnsTC, cpiOpnsData);
	saveData(adjTC, cpiAdjData);
	vector<vect> minorityFracsTC = findMinorityFractions(opns);
	vector<vect> conflictsTC = findConflicts(adjMatrices, opns);
	double newOpns = project(times, average(minorityFractsTC));
	int newConlicts = (int) (project(times, average(conflictsTC)) + 0.5);
	for(i = 0; i < nVms; i++) {
	  vms[i].initGraph(newDist, newConflicts);
	  adjMatrices[i].clear();
	  opns[i].clear();
	  times[i].clear();
	}
	microStepCount = 0;
	step += projStep;
      }
    }
  }
  cpiAdjData.close();
  cpiOpnsData.close();
}

template <typename T>
double votingModel::average(const vector<T> &data) {
  T sum = 0;
  for(vector<T>::iterator v = data.begin(); v != data.end(); v++) {
    sum += (*v);
  }
  return sum/data.size();
}
template double votingModelCPI::average<int>(const vector<int> &data);

vect votingModel::average(const vector<vect> &data) {
  vect avgdData;
  vector<int> temp;
  int nVs = data.size();
  int nElements = data[0].size();
  int i;
  for(i = 0; i < nElements; i++) {
    for(vector<vect>::iterator v = data.begin(); v != data.end(); v++) {
      temp.push_back(v[i]);
    }
    avgdData.push_back(average<int>(temp) + 0.5);
    temp.clear();
  }
}

vmVects votingModel::average(const vector<vmVects> &data) {
  vmVects avgdData;
  for(vector<vmVects>::iterator vmVs = data.begin(); vmVs != data.end(); vmVs++) {
    avgdData.push_back(average(*vmVs));
  }
  return avgdData();
}

vmMatrices votingModel::average(const std::vector<vmMatrices> &data) {
  vmMatrics avgdData;
  matrix M;
  vect v;
  int n = data[0][0].size();
  int i, j;
  for(vector<vmMatrices>::iterator vmMs = data.begin(); vmMs != data.end(); vmMs++) {
    avgdData.push_back(M);
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
	for(vector<matrix>::iterator m = vmMs->begin(); m != v->end(); m++) {
	  v.push_back((*m)[i][j]);
	}
	avgdData.back().push_back
  
  

double votingModelCPI::project(vect data, vect times) {
  
}
