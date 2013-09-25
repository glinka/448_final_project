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


votingModelCPI::votingModelCPI(vector<votingModel> vms, int waitingPeriod, int collectionInterval, int nMicroSteps, string file_header, string file_name): vms(vms), waitingPeriod(waitingPeriod), collectionInterval(collectionInterval), nMicroSteps(nMicroSteps), file_header(file_header), file_name(file_name) {
};

int votingModelCPI::run(long int nSteps, int proj_step, int save_data_interval) {
  long int step = 0;
  int microStepCount = 0;
  stringstream ss;
  string folder = "./csv_data/";
  ss << folder << "CPIAdj_" << file_name << ".csv";
  string adjFilename = ss.str();
  ofstream cpiAdjData;
  cpiAdjData.open(adjFilename);
  ss.str("");
  ss << folder << "CPIOpns_" << file_name << ".csv";
  string opnsFilename = ss.str();
  ofstream cpiOpnsData;
  cpiOpnsData.open(opnsFilename);
  ss.str("");
  ss << folder << "CPITimes_" << file_name << ".csv";
  string timesFilename = ss.str();
  ofstream cpiTimesData;
  cpiTimesData.open(timesFilename);
  ss.str("");
  ss << folder << "CPINvms_" << file_name << ".csv";
  string nvmsFilename = ss.str();
  ofstream cpiNvmsData;
  cpiNvmsData.open(nvmsFilename);
  cpiAdjData << file_header << endl;
  cpiOpnsData << file_header << endl;
  cpiTimesData << file_header << endl;
  cpiNvmsData << file_header << endl;
  //initialize space for each vm's time-courses, all will share the same time vector
  vect times_to_save;
  vect nvms_to_save;
  vector<vmMatrices> adj_to_save;
  vector<vmVects> opns_to_save;
  unsigned int nCompletedVMs = 0;
  while(step < nSteps) {
    for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
      if(vm->getConflicts() > 0) {
	vm->step();
      }
      else {
	//remove finished runs from simulation
	vms.erase(vm);
	nCompletedVMs++;
      }
    }
    step++;
    microStepCount++;
    //collect data to save every save_data_interval steps, and also at the 
    //beginning of each projection iteration and end of simulation
    if(microStepCount % save_data_interval == 0 || microStepCount == 1 || nCompletedVMs == vms.size()) {
      adj_to_save.push_back(vector<matrix>());
      opns_to_save.push_back(vector<vect>());
      for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
	adj_to_save.back().push_back(vm->getAdjMatrix());
	opns_to_save.back().push_back(vm->getOpns());
      }
      times_to_save.push_back(step);
      nvms_to_save.push_back(vms.size());
    }
    if(nCompletedVMs == vms.size()) {
      this->saveData(adj_to_save, cpiAdjData);
      this->saveData(opns_to_save, cpiOpnsData);
      this->saveData(times_to_save, cpiTimesData);
      this->saveData(nvms_to_save, cpiNvmsData);
      cpiAdjData.close();
      cpiOpnsData.close();
      cpiTimesData.close();
      cpiNvmsData.close();
      return 0;
    }
    if(microStepCount > waitingPeriod) {
      int onManifoldStep = microStepCount - waitingPeriod;
      if(onManifoldStep % collectionInterval == 0) {
	adjMatrices.push_back(vector<matrix>());
	opns.push_back(vector<vect>());
	for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
	  adjMatrices.back().push_back(vm->getAdjMatrix());
	  opns.back().push_back(vm->getOpns());
	}
	times.push_back(onManifoldStep);
      }
      if(microStepCount == nMicroSteps) {
	//just project the minority fraction you fool
	vector<double> minorityFracsTC = findAvgdMinorityFractions(opns);
	vect conflictsTC = findAvgdConflicts(adjMatrices, opns);
	double newMinorityFrac = project<double>(times, minorityFracsTC, proj_step);
	//hard coded for two opinions
	int newConflicts = (int) (project<int>(times, conflictsTC, proj_step) + 0.5);
	int temp_proj_step = proj_step;
	while(newConflicts < 0 || newMinorityFrac < 0) {
	  temp_proj_step /= 2;
	  minorityFracsTC = findAvgdMinorityFractions(opns);
	  conflictsTC = findAvgdConflicts(adjMatrices, opns);
	  newMinorityFrac = project<double>(times, minorityFracsTC, temp_proj_step);
	  newConflicts = (int) (project<int>(times, conflictsTC, temp_proj_step) + 0.5);
	}
	double newDist[2] = {newMinorityFrac, 1 - newMinorityFrac};
	for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
	  vm->initGraph(newDist, newConflicts);
	}
	this->saveData(adj_to_save, cpiAdjData);
	this->saveData(opns_to_save, cpiOpnsData);
	this->saveData(times_to_save, cpiTimesData);
	this->saveData(nvms_to_save, cpiNvmsData);
	adj_to_save.clear();
	opns_to_save.clear();
	times_to_save.clear();
	nvms_to_save.clear();
	adjMatrices.clear();
	opns.clear();
	times.clear();
	microStepCount = 0;
	step += proj_step;
      }
    }
  }
  this->saveData(adj_to_save, cpiAdjData);
  this->saveData(opns_to_save, cpiOpnsData);
  this->saveData(times_to_save, cpiTimesData);
  this->saveData(nvms_to_save, cpiNvmsData);
  cpiAdjData.close();
  cpiOpnsData.close();
  cpiTimesData.close();
  cpiNvmsData.close();
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
double votingModelCPI::project(const vect &times, const vector<T> &data, const int proj_step) {
  fxs linearFit;
  linearFit.push_back([] (double x) {return 1.0;});
  linearFit.push_back([] (double x) {return x;});
  vector<double> doubleTimes(times.begin(), times.end());
  vector<double> doubleData(data.begin(), data.end());
  vector<double> coeffs = fitCurves::fitFx(doubleTimes, doubleData, linearFit);
  double newVal = 0;
  vector<double>::iterator coeff = coeffs.begin();
  for(fxs::iterator f = linearFit.begin(); f != linearFit.end(); f++) {
    newVal += (*coeff)*((*f)(times.back() + proj_step));
    coeff++;
  }
  return newVal;
}
template double votingModelCPI::project<int>(const vect &times, const vector<int> &data, const int proj_step);
template double votingModelCPI::project<double>(const vect &times, const vector<double> &data, const int proj_step);
