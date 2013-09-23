#include <fstream>
#include <string>
#include <sstream>
#include "fitCurves.h"
#include "votingModel.h"
#include "votingModelCPI.h"

using namespace std;

typedef vector<votingModel>::iterator vmIt;

template<typename T>
bool inArray(T val, vector<T> v) {
  typedef typename vector<T>::iterator myIt;
  for(myIt vi = v.begin(); vi != v.end(); vi++) {
    if(val == *vi) {
      return true;
    }
  }
  return false;
}
template bool inArray<int>(int val, vector<int> v);

template<typename T>
int whereIs(T val, vector<T> v) {
  typedef typename vector<T>::iterator myIt;
  int count = 0;
  for(myIt vi = v.begin(); vi != v.end(); vi++) {
    if(*vi == val) {
      return count;
    }
    count++;
  }
  return -1;
}
template int whereIs<int>(int val, vector<int> v);


votingModelCPI::votingModelCPI(vector<votingModel> vms, const int projStep, int waitingPeriod, int collectionInterval, int nMicroSteps): vms(vms), projStep(projStep), waitingPeriod(waitingPeriod), collectionInterval(collectionInterval), nMicroSteps(nMicroSteps) {
};

void votingModelCPI::run(int nSteps) {
  int step = 0;
  int microStepCount = 0;
  int i;
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
	vms[0].step();
	(*vm).step();
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
	this->saveData(opnsTC, cpiOpnsData);
	this->saveData(adjTC, cpiAdjData);
	vector<double> minorityFracsTC = findAvgdMinorityFractions(opns);
	vect conflictsTC = findAvgdConflicts(adjMatrices, opns);
	double newOpns = project<double>(times, minorityFracsTC);
	int newConflicts = (int) (project<int>(times, conflictsTC) + 0.5);
	for(i = 0; i < nVms; i++) {
	  vms[i].initGraph(&newOpns, newConflicts);
	  adjMatrices[i].clear();
	  opns[i].clear();
	}
	times.clear();
	microStepCount = 0;
	step += projStep;
      }
    }
  }
  cpiAdjData.close();
  cpiOpnsData.close();
}

template <typename T>
double votingModelCPI::average(const vector<T> &data) {
  T sum = 0;
  typedef typename vector<T>::const_iterator myIt;
  for(myIt v = data.begin(); v != data.end(); v++) {
    sum += (*v);
  }
  return sum/data.size();
}
template double votingModelCPI::average<int>(const vector<int> &data);

vect votingModelCPI::average(const vector<vect> &data) {
  vect avgdData;
  vector<int> temp;
  int nElements = data[0].size();
  int i;
  for(i = 0; i < nElements; i++) {
    for(vector<vect>::const_iterator v = data.begin(); v != data.end(); v++) {
      temp.push_back((*v)[i]);
    }
    avgdData.push_back(average<int>(temp) + 0.5);
    temp.clear();
  }
  return avgdData;
}

vmVects votingModelCPI::average(const vector<vmVects> &data) {
  vmVects avgdData;
  for(vector<vmVects>::const_iterator vmVs = data.begin(); vmVs != data.end(); vmVs++) {
    avgdData.push_back(average(*vmVs));
  }
  return avgdData;
}

vmMatrices votingModelCPI::average(const std::vector<vmMatrices> &data) {
  vmMatrices avgdData;
  matrix M;
  vect v;
  int n = data[0][0].size();
  int i, j;
  for(vector<vmMatrices>::const_iterator vmMs = data.begin(); vmMs != data.end(); vmMs++) {
    avgdData.push_back(M);
    for(i = 0; i < n; i++) {
      avgdData.back().push_back(v);
      for(j = 0; j < n; j++) {
	for(vector<matrix>::const_iterator m = (*vmMs).begin(); m != (*vmMs).end(); m++) {
	  v.push_back((*m)[i][j]);
	}
	avgdData.back().back().push_back(average<int>(v) + 0.5);
	v.clear();
      }
    }
  }
  return avgdData;
}

void votingModelCPI::saveData(const vector<vect> &data, ofstream &fileHandle) {
  for(vector<vect>::const_iterator v = data.begin(); v != data.end(); v++) {
    for(vect::const_iterator i = (*v).begin(); i != (*v).end(); i++) {
      fileHandle << *i << endl;
    }
  }
}

  
void votingModelCPI::saveData(const vector<matrix> &data, ofstream &fileHandle) {
  for(vector<matrix>::const_iterator m = data.begin(); m != data.end(); m++) {
    for(matrix::const_iterator row = (*m).begin(); row != (*m).end(); row++) {
      //row may not actually point to rows in the memory layout, but as our data is symmetric
      //it doesn't matter if it's looping over rows or columns
      for(vect::const_iterator i = (*row).begin(); i != (*row).end(); i++) {
	fileHandle << *i;
	if(i != --((*row).end())) {
	  fileHandle << ",";
	}
      }
      fileHandle << endl;
    }
  }
}

vector<double> votingModelCPI::findAvgdMinorityFractions(const vector<vmVects> &opns) {
  vector<vect> avgdOpnsTC = average(opns);
  vector<double> avgdMinorityFrac;
  int n = avgdOpnsTC[0].size();
  int index;
  vect ids;
  vect counts;
  for(vector<vect>::const_iterator v = avgdOpnsTC.begin(); v != avgdOpnsTC.end(); v++) {
    for(vect::const_iterator i = (*v).begin(); i != (*v).end(); i++) {
      if(!inArray<int>(*i, ids)) {
	ids.push_back(*i);
	index = ids.size() - 1;
      }
      else {
	index = whereIs<int>(*i, ids);
      }
      counts[index]++;
    }
    int max = 0;
    for(vect::iterator c = counts.begin(); c != counts.end(); c++) {
      if(max > *c) {
	max = *c;
      }
    }
    avgdMinorityFrac.push_back(1.0*max/n);
  }
  return avgdMinorityFrac;
}

vect votingModelCPI::findAvgdConflicts(const vector<vmMatrices> &As, const vector<vmVects> &opns) {
  vector<matrix> avgdAs = average(As);
  vector<vect> avgdOpns = average(opns);
  vector<vect>::iterator opn = avgdOpns.begin();
  vect conflictsV;
  int n = opns[0][0].size();
  for(vector<matrix>::iterator M = avgdAs.begin(); M != avgdAs.end(); M++) {
    vect currentOpns = *opn;
    opn++;
    int conflicts = 0;
    int i, j;
    for(i = 0; i < n; i++) {
      int currentOpn = currentOpns[i];
      for(j = i; j < n; j++) {
	if(currentOpn != currentOpns[j]) {
	  conflicts += (*M)[i][j];
	}
      }
    }
    conflictsV.push_back(conflicts);
  }
  return conflictsV;
}    

template <typename T>
double votingModelCPI::project(const vect &times, const vector<T> &data) {
  fxs linearFit;
  linearFit.push_back([] (double x) {return 1.0;});
  linearFit.push_back([] (double x) {return x;});
  vector<double> doubleTimes(times.begin(), times.end());
  vector<double> doubleData(data.begin(), data.end());
  vector<double> coeffs = fitCurves::fitFx(doubleTimes, doubleData, linearFit);
  int newVal = 0;
  vector<double>::iterator coeff = coeffs.begin();
  for(fxs::iterator f = linearFit.begin(); f != linearFit.end(); f++) {
    newVal += (*coeff)*((*f)(times.back() + projStep));
    coeff++;
  }
  return newVal;
}
template double votingModelCPI::project<int>(const vect &times, const vector<int> &data);
template double votingModelCPI::project<double>(const vect &times, const vector<double> &data);
