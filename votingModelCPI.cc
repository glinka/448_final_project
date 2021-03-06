#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
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
  const int MAX_PROJECTION_ATTEMPTS = 10;
  int nvms = vms.size();
  long int step = 0;
  int microStepCount = 0;
  //open all the files
  stringstream ss;
  int folder_counter = 0;
  string dir_base = "./csv_data_fast";
  string dir;
  bool isdir = true;
  struct stat stat_dir;
  do {
    ss.str("");
    ss << dir_base << folder_counter << "/";
    folder_counter++;
    dir = ss.str();
    int check = stat(dir.c_str(), &stat_dir);
    if(check == -1) {
      mkdir(dir.c_str(), 0700);
      isdir = false;
    }
  } while (isdir);
  ss.str("");
  ss << dir << "CPIAdj_" << file_name << ".csv";
  string adjFilename = ss.str();
  ofstream cpiAdjData;
  cpiAdjData.open(adjFilename);
  ss.str("");
  ss << dir << "CPIOpns_" << file_name << ".csv";
  string opnsFilename = ss.str();
  ofstream cpiOpnsData;
  cpiOpnsData.open(opnsFilename);
  ss.str("");
  ss << dir << "CPIConflicts_" << file_name << ".csv";
  string conflictsFilename = ss.str();
  ofstream cpiConflictsData;
  cpiConflictsData.open(conflictsFilename);
  ss.str("");
  ss << dir << "CPITimes_" << file_name << ".csv";
  string timesFilename = ss.str();
  ofstream cpiTimesData;
  cpiTimesData.open(timesFilename);
  ss.str("");
  ss << dir << "CPIids_" << file_name << ".csv";
  string idsFilename = ss.str();
  ofstream cpiIdsData;
  cpiIdsData.open(idsFilename);
  ss.str("");
  ss << dir << "CPIMinorities_" << file_name << ".csv";
  string minsFilename = ss.str();
  ofstream cpiMinsData;
  cpiMinsData.open(minsFilename);
  cpiAdjData << file_header << endl;
  cpiOpnsData << file_header << endl;
  cpiTimesData << file_header << endl;
  cpiIdsData << file_header << endl;
  cpiMinsData << file_header << endl;
  //initialize space for each vm's time-courses, all will share the same time vector
  vect times_to_save;
  vector< vect > ids_to_save;
  vector< vector<double> > mins_to_save;
  //far too large a file, don't save for now. would be nice to have a non-commented out way of adjusting what was saved
  //  vector<vmMatrices> // adj_to_save;
  vector<vmVects> opns_to_save;
  vector< vect > conflicts_to_save;
  int nCompletedVMs = 0;
  int init_waiting_period = 25000;
  bool first_projection = true;
  for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
    vm->initGraph();
  }
  vector< vmIt > to_delete;
  while(step < nSteps) {
    int finishedVMs = 0;
    for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
      if(vm->getConflicts() > 0) {
	vm->step();
      }
      else {
	//remove finished runs from simulation
	to_delete.push_back(vm);
	finishedVMs++;
      }
    }
    //    nCompletedVMs += finishedVMs;
    step++;
    microStepCount++;
    //collect data to save every save_data_interval steps, and also at the 
    //beginning of each projection iteration and end of simulation
    if(microStepCount % save_data_interval == 0 || microStepCount == 1 || finishedVMs == nvms) {
      // adj_to_save.push_back(vector<matrix>());
      opns_to_save.push_back(vector<vect>());
      ids_to_save.push_back(vector<int>());
      mins_to_save.push_back(vector<double>());
      conflicts_to_save.push_back(vector<int>());
      for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
	// adj_to_save.back().push_back(vm->getAdjMatrix());
	opns_to_save.back().push_back(vm->getOpns());
	ids_to_save.back().push_back(vm->get_id());
	mins_to_save.back().push_back(vm->getMinorityFraction());
	conflicts_to_save.back().push_back(vm->getConflicts());
      }
      times_to_save.push_back(step);
      //      if(nCompletedVMs == nvms) {
      if(finishedVMs == nvms) {
	//	this->saveData(// adj_to_save, cpiAdjData);
	this->saveData(opns_to_save, cpiOpnsData);
	this->saveData(times_to_save, cpiTimesData);
	this->saveData(ids_to_save, cpiIdsData);
	this->saveData(mins_to_save, cpiMinsData);
	this->saveData(conflicts_to_save, cpiConflictsData);
	cpiAdjData.close();
	cpiOpnsData.close();
	cpiTimesData.close();
	cpiIdsData.close();
	cpiMinsData.close();
	cpiConflictsData.close();
	return 0;
      }
    }
    //******************************
    //faster to check if size > 0 ?
    //******************************
    /**
    if(to_delete.size() > 0) {
      for(vector< vmIt >::const_iterator vm = to_delete.begin(); vm != to_delete.end(); vm++) {
	vms.erase(*vm);
      }
      to_delete.clear();
    }
    **/
    if((microStepCount > waitingPeriod) && (microStepCount > init_waiting_period)) {
      int onManifoldStep = 1;
      if(first_projection) {
	microStepCount = 0;
	init_waiting_period = 0;
	first_projection = false;
      }
      else {
	onManifoldStep = microStepCount - waitingPeriod;
      }
      if(onManifoldStep % collectionInterval == 0) {
	//	adjMatrices.push_back(vector<matrix>());
	//	cpi_opns.push_back(vector<vect>());
	cpi_conflicts.push_back(vector<int>());
	cpi_minorities.push_back(vector<double>());
	for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
	  //	  adjMatrices.back().push_back(vm->getAdjMatrix());
	  //	  cpi_opns.back().push_back(vm->getOpns());
	  cpi_conflicts.back().push_back(vm->getConflicts());
	  cpi_minorities.back().push_back(vm->getMinorityFraction());
	}
	times.push_back(step);
      }
      if(microStepCount == nMicroSteps) {
	//just project the minority fraction you fool
	//	vector<double> minorityFracsTC = findAvgdMinorityFractions(cpi_opns);
	//	vect conflictsTC = findAvgdConflicts(adjMatrices, cpi_opns);
	vector<double> minorityFracsTC = average_mf(cpi_minorities);
	vector< int > conflictsTC = easy_average<int>(cpi_conflicts);
	double newMinorityFrac = project<double>(times, minorityFracsTC, proj_step);
	//hard coded for two opinions
	int newConflicts = (int) (project<int>(times, conflictsTC, proj_step) + 0.5);
	//	int newConflicts = conflictsTC.back();
	int temp_proj_step = proj_step;
	int nattempts = 0;
	while((newConflicts <= 0 || newMinorityFrac <= 0) && (nattempts < MAX_PROJECTION_ATTEMPTS)) {
	  temp_proj_step /= 2;
	  //	  minorityFracsTC = findAvgdMinorityFractions(cpi_opns);
	  //	  conflictsTC = findAvgdConflicts(adjMatrices, cpi_opns);
	  newMinorityFrac = project<double>(times, minorityFracsTC, temp_proj_step);
	  newConflicts = (int) (project<int>(times, conflictsTC, temp_proj_step) + 0.5);
	  nattempts++;
	}
	cout << nattempts << endl;
	if(nattempts == MAX_PROJECTION_ATTEMPTS) {
	  cout << "failed projection" << endl;
	  return -1;
	}
	double newDist[2] = {newMinorityFrac, 1 - newMinorityFrac};
	for(vmIt vm = vms.begin(); vm != vms.end(); vm++) {
	  //	  if(vm->getConflicts() > 0) {
	    vm->initGraph(newDist, newConflicts);
	    //	  }
	}
//	this->saveData(// adj_to_save, cpiAdjData);
	this->saveData(opns_to_save, cpiOpnsData);
	this->saveData(times_to_save, cpiTimesData);
	this->saveData(ids_to_save, cpiIdsData);
	this->saveData(mins_to_save, cpiMinsData);
	this->saveData(conflicts_to_save, cpiConflictsData);
	// adj_to_save.clear();
	opns_to_save.clear();
	times_to_save.clear();
	ids_to_save.clear();
	mins_to_save.clear();
	conflicts_to_save.clear();
	adjMatrices.clear();
	cpi_opns.clear();
	cpi_minorities.clear();
	cpi_conflicts.clear();
	times.clear();
	microStepCount = 0;
	step += proj_step;
      }
    }
  }
    //  this->saveData(// adj_to_save, cpiAdjData);
  this->saveData(opns_to_save, cpiOpnsData);
  this->saveData(times_to_save, cpiTimesData);
  this->saveData(ids_to_save, cpiIdsData);
  this->saveData(mins_to_save, cpiMinsData);
  this->saveData(conflicts_to_save, cpiConflictsData);
  cpiAdjData.close();
  cpiOpnsData.close();
  cpiTimesData.close();
  cpiIdsData.close();
  cpiConflictsData.close();
  cpiMinsData.close();
  return 0;
}


vector<double> votingModelCPI::average_mf(const vector< vector<double> > &data) {
  vector<double> avgs;
  for(vector< vector<double> >::const_iterator v = data.begin(); v != data.end(); v++) {
    double n = 0.;
    double sum = 0.;
    for(vector<double>::const_iterator val = (*v).begin(); val != (*v).end(); val++) {
      n++;
      sum += (*val);
    }
    avgs.push_back(sum/n);
  }
  return avgs;
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

template <typename T>
vector< T > votingModelCPI::easy_average(const vector< vector< T > > &data) {
  vector < T > avgs;
  typedef typename vector < vector< T > >::const_iterator it;
  for(it v = data.begin(); v != data.end(); v++) {
    avgs.push_back(average<T>(*v));
  }
  return avgs;
}
template vector< int > votingModelCPI::easy_average<int>(const vector< vector< int > > &data);
template vector< double > votingModelCPI::easy_average<double>(const vector< vector< double > > &data);

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

void votingModelCPI::saveData(const vector< vect > &data, ofstream &fileHandle) {
  for(vector< vect >::const_iterator v = data.begin(); v != data.end(); v++) {
    for(vect::const_iterator i = (*v).begin(); i != (*v).end(); i++) {
      fileHandle << *i;
      //add comma if not last element
      if(i != ((*v).end() - 1)) {
	fileHandle << ",";
      }
    }
    //don't add newline if content wasn't added
    if(!(*v).empty()) {
      fileHandle << endl;
    }
  }
}

void votingModelCPI::saveData(const vector< vector<double> > &data, ofstream &fileHandle) {
  for(vector< vector<double> >::const_iterator v = data.begin(); v != data.end(); v++) {
    for(vector< double >::const_iterator i = (*v).begin(); i != (*v).end(); i++) {
      fileHandle << *i;
      //add comma if not last element
      if(i != ((*v).end() - 1)) {
	fileHandle << ",";
      }
    }
    //don't add newline if content wasn't added
    if(!(*v).empty()) {
      fileHandle << endl;
    }
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
  linearFit.push_back([] (double x) {return x;});
  linearFit.push_back([] (double x) {return 1.0;});
  vector<double> doubleTimes(times.begin(), times.end());
  //shift to 0
  double max_time = doubleTimes.back();
  /*
  for(vector<double>::iterator t = doubleTimes.begin(); t != doubleTimes.end(); t++) {
    (*t) -= max_time;
  } 
  */
  vector<double> doubleData(data.begin(), data.end());
  vector<double> coeffs = fitCurves::fitFx(doubleTimes, doubleData, linearFit);
  double newVal = 0;
  vector<double>::iterator coeff = coeffs.begin();
  for(fxs::iterator f = linearFit.begin(); f != linearFit.end(); f++) {
    cout << (*coeff) << endl;
    newVal += (*coeff)*((*f)(times.back() + proj_step));
    coeff++;
  }
  cout << newVal << endl;
  return newVal;
}
template double votingModelCPI::project<int>(const vect &times, const vector<int> &data, const int proj_step);
template double votingModelCPI::project<double>(const vect &times, const vector<double> &data, const int proj_step);
