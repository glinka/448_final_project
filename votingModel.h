#ifndef VOTINGMODEL_H
#define VOTINGMODEL_H
#include <string>
#include <vector>
#include <chrono>
#include <random>

class votingModelCPI;

class votingModel {
 private:
  const double ROUND_CONST;
  const int n, k, collectionInterval;
  const long int maxIter;
  const double a, avgDeg;
  double *initDist;
  int *degs, *Opns, *opnCounts, *nConflicts;
  int **A;
  double rnNormalization;
  int conflicts;
  //tracks the number of conflicting edges at each vertex
  std::string rewireTo, fileName;
  std::mt19937 *mt;
  void initGraph(double *dist);
  double genURN();
 public:
  int consistencyCheck();
  int vote();
  void step();
  void initGraph(double *dist, int conflicts);
  template <typename dataType>
    void saveData(const std::vector<dataType> data, std::ofstream &fileHandle);
  template <typename dataType>
    void saveData(const std::vector<std::vector<dataType> > data, std::ofstream &fileHandle);
  votingModel(int n, int k, long int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, std::string rewireTo, std::string fileName="");

  int getConflicts() {
    return conflicts;
  };

  std::vector<int> getOpns() {
    std::vector<int> newOpns;
    for(int i = 0; i < n; i++) {
      newOpns.push_back(Opns[i]);
    }
    return newOpns;
  };

  std::vector<std::vector<int> > getAdjMatrix() {
    std::vector<std::vector<int> > newA;
    std::vector<int> v;
    for(int i = 0; i < n; i++) {
      newA.push_back(v);
      for(int j = 0; j < n; j++) {
	newA.back().push_back(A[i][j]);
      }
    }
    return newA;
  };

  ~votingModel() {
    for(int i = 0; i < n; i++) {
      delete[] A[i];
    }
    delete[] A;
    delete[] degs;
    delete[] Opns;
    delete[] opnCounts;
    delete[] nConflicts;
    delete mt;
  };

 votingModel(const votingModel &toCopy): ROUND_CONST(toCopy.ROUND_CONST), n(toCopy.n), k(toCopy.k), collectionInterval(toCopy.collectionInterval), maxIter(toCopy.maxIter), a(toCopy.a), avgDeg(toCopy.avgDeg) {
      A = new int*[n];
      for(int i = 0; i < n; i++) {
	  A[i] = new int[n];
      }
      degs = new int[n];
      Opns = new int[n];
      nConflicts = new int[n];
      opnCounts = new int[k];
      for(int i = 0; i < n; i++) {
	  degs[i] = toCopy.degs[i];
	  Opns[i] = toCopy.Opns[i];
	  nConflicts[i] = toCopy.nConflicts[i];
	  for(int j = 0; j < n; j++) {
	      A[i][j] = toCopy.A[i][j];
	  }
      }
      for(int i = 0; i < k; i++) {
	  opnCounts[i] = toCopy.opnCounts[i];
      }
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      mt = new std::mt19937(seed);
      rnNormalization = (double) (mt->max()+1);
      conflicts = toCopy.conflicts;
      rewireTo = toCopy.rewireTo;
      fileName = toCopy.fileName;
      initDist = toCopy.initDist;
  };

  votingModel &operator=(const votingModel &rhs) {
      for(int i = 0; i < n; i++) {
	  degs[i] = rhs.degs[i];
	  Opns[i] = rhs.Opns[i];
	  nConflicts[i] = rhs.nConflicts[i];
	  for(int j = 0; j < n; j++) {
	      A[i][j] = rhs.A[i][j];
	  }
      }
      for(int i = 0; i < k; i++) {
	  opnCounts[i] = rhs.opnCounts[i];
      }
      conflicts = rhs.conflicts;
      rewireTo.assign(rhs.rewireTo);
      fileName.assign(rhs.fileName);
      initDist = rhs.initDist;
      return *this;
  };
};

#endif
