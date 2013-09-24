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

int votingModelCPI::run(long int nSteps) {
  long int step = 0;
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
  ss.str("");
  ss << "CPITimes_nSteps_" << nSteps << "_projStep_" << projStep << ".csv";
  string timesFilename = ss.str();
  ofstream cpiTimesData;
  cpiTimesData.open(timesFilename);
  bool addHeader = true;
  vector<matrix> tempM;
  vector<vect> tempV;
  vect real_times;
  while(step < nSteps) {
    unsigned int nCompletedVMs = 0;
    for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
      if(vm->getConflicts() > 0) {
	vm->step();
      }
      else {
	nCompletedVMs++;
      }
    }
    if(nCompletedVMs == vms.size()) {
      return 0;
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
	real_times.push_back(step);
	if(addHeader) {
	  int n = adjMatrices[0][0].size();
	  cpiAdjData << "n=" << n << ",nVms=" << nVms << endl;
	  cpiOpnsData << "n=" << n << ",nVms=" << nVms << endl;
	  cpiTimesData << "n=" << n << ",nVms=" << nVms << endl;
	  addHeader = false;
	}
      }
      if(microStepCount == nMicroSteps) {
	//just project the minority fraction you fool
	this->saveData(opns, cpiOpnsData);
	this->saveData(adjMatrices, cpiAdjData);
	this->saveData(real_times, cpiTimesData);
	vector<double> minorityFracsTC = findAvgdMinorityFractions(opns);
	vect conflictsTC = findAvgdConflicts(adjMatrices, opns);
	double newMinorityFrac = project<double>(times, minorityFracsTC);
	double newDist[2] = {newMinorityFrac, 1 - newMinorityFrac};
	int newConflicts = (int) (project<int>(times, conflictsTC) + 0.5);
	if(newConflicts < 0) {
	  return 2;
	}
	for(i = 0; i < nVms; i++) {
	  vms[i].initGraph(newDist, newConflicts);
	}
	adjMatrices.clear();
	opns.clear();
	times.clear();
	real_times.clear();
	microStepCount = 0;
	step += projStep;
      }
    }
  }
  cpiAdjData.close();
  cpiOpnsData.close();
  return 0;
}

template <typename T>
double votingModelCPI::average(const vector<T> &data) {
  T sum = 0;
  typedef typename vector<T>::const_iterator myIt;
  for(myIt v = data.begin(); v != data.end(); v++) {
    sum += (*v);
  }
  return 1.0*sum/data.size();
}
template double votingModelCPI::average<int>(const vector<int> &data);
template double votingModelCPI::average<double>(const vector<double> &data);

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

vmMatrices votingModelCPI::average(const vector<vmMatrices> &data) {
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

void votingModelCPI::saveData(const vect &data, ofstream &fileHandle) {
    for(vect::const_iterator i = data.begin(); i != data.end(); i++) {
	fileHandle << *i << endl;
    }
}

void votingModelCPI::saveData(const vector<vmVects>  &data, ofstream &fileHandle) {
  for(vector<vmVects>::const_iterator vmV = data.begin(); vmV != data.end(); vmV++) {
    for(vector<vect>::const_iterator v = (*vmV).begin(); v != (*vmV).end(); v++) {
      for(vect::const_iterator i = (*v).begin(); i != (*v).end(); i++) {
	fileHandle << *i << endl;
      }
    }
  }
}

  
void votingModelCPI::saveData(const vector<vmMatrices> &data, ofstream &fileHandle) {
  for(vector<vmMatrices>::const_iterator vmM = data.begin(); vmM != data.end(); vmM++) {
    for(vector<matrix>::const_iterator m = (*vmM).begin(); m != (*vmM).end(); m++) {
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
}

vector<double> votingModelCPI::findAvgdMinorityFractions(const vector<vmVects> &opns) {
  vector<vector<double> > minorityFracs;
  vector<double> v;
  for(vector<vmVects>::const_iterator vmV = opns.begin(); vmV != opns.end(); vmV++) {
    minorityFracs.push_back(v);
    for(vmVects::const_iterator opnV = (*vmV).begin(); opnV != (*vmV).end(); opnV++) {
      minorityFracs.back().push_back(getMinorityFrac(*opnV));
    }
  }
  vector<double> avgdMinorityFracs;
  for(vector<vector<double> >::iterator mfV = minorityFracs.begin(); mfV != minorityFracs.end(); mfV++) {
    avgdMinorityFracs.push_back(average<double>(*mfV));
  }
  return avgdMinorityFracs;
}

double votingModelCPI::getMinorityFrac(const vect &opns) {
  vect counts;
  vect ids;
  for(vect::const_iterator i = opns.begin(); i != opns.end(); i++) {
    if(inArray<int>(*i, ids)) {
      counts[whereIs<int>(*i, ids)]++;
    }
    else {
      ids.push_back(*i);
      counts.push_back(1);
    }
  }
  int min = opns.size();
  for(vect::iterator i = counts.begin(); i != counts.end(); i++) {
    if(*i < min) {
      min = *i;
    }
  }
  return 1.0*min/opns.size();
}

vect votingModelCPI::findAvgdConflicts(const vector<vmMatrices> &As, const vector<vmVects> &opns) {
  vector<vect> vmConflicts;
  vector<vmVects>::const_iterator vmVs = opns.begin();
  for(vector<vmMatrices>::const_iterator vmMs = As.begin(); vmMs != As.end(); vmMs++) {
    vmConflicts.push_back(vect());
    vmVects::const_iterator V = (*vmVs).begin();
    for(vmMatrices::const_iterator M = (*vmMs).begin(); M != (*vmMs).end(); M++) {
      vmConflicts.back().push_back(getConflicts(*M, *V));
      V++;
    }
    vmVs++;
  }
  vect avgdConflicts;
  for(vector<vect>::iterator V = vmConflicts.begin(); V != vmConflicts.end(); V++) {
    avgdConflicts.push_back((int) average<int>(*V));
  }
  return avgdConflicts;
}

int votingModelCPI::getConflicts(const matrix &A, const vect &opns) {
  int conflicts = 0;
  int n = opns.size();
  for(int i = 0; i < n; i++) {
    int currentOpn = opns[i];
    for(int j = i+1; j < n; j++) {
      if(A[i][j] > 0 && opns[j] != currentOpn) {
	conflicts++;
      }
    }
  }
  return conflicts;
}

template <typename T>
double votingModelCPI::project(const vect &times, const vector<T> &data) {
  fxs linearFit;
  linearFit.push_back([] (double x) {return 1.0;});
  linearFit.push_back([] (double x) {return x;});
  vector<double> doubleTimes(times.begin(), times.end());
  vector<double> doubleData(data.begin(), data.end());
  vector<double> coeffs = fitCurves::fitFx(doubleTimes, doubleData, linearFit);
  double newVal = 0;
  vector<double>::iterator coeff = coeffs.begin();
  for(fxs::iterator f = linearFit.begin(); f != linearFit.end(); f++) {
    newVal += (*coeff)*((*f)(times.back() + projStep));
    coeff++;
  }
  return newVal;
}
template double votingModelCPI::project<int>(const vect &times, const vector<int> &data);
template double votingModelCPI::project<double>(const vect &times, const vector<double> &data);
