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
  int id;
  static int id_tracker;
  const double *initDist;
  int *degs, *Opns, *opnCounts, *nConflicts;
  int **A;
  double rnNormalization;
  int conflicts;
  //tracks the number of conflicting edges at each vertex
  std::string rewireTo, fileName;
  std::mt19937 *mt;
  double genURN();
  void alter_opinions();
 public:
  int consistencyCheck();
  int vote(bool alter, const std::string run_id);
  std::vector< std::vector<int> > run_nsteps(const int nsteps, bool& converged);
  void step();
  void initGraph();
  std::vector<int> get_opns();
  void initGraph(double *dist, int conflicts);
  template <typename dataType>
    void saveData(const std::vector<dataType> data, std::ofstream &fileHandle);
  template <typename dataType>
    void saveData(const std::vector<std::vector<dataType> > data, std::ofstream &fileHandle);
  votingModel(const int n, const double avgDeg, const int k, const double *initDist, const double a, const std::string rewireTo);
  votingModel(int n, int k, long int maxIter, int collectionInterval, double a, double avgDeg, const double *initDist, std::string rewireTo, std::string fileName);

  int get_id() {
    return id;
  };

  int getConflicts() {
    return conflicts;
  };

  double getMinorityFraction() {
    int min = n;
    for(int i = 0; i < k; i++) {
      if(opnCounts[i] < min) {
	min = opnCounts[i];
      }
    }
    return (1.0*min)/n;
  }

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
    A = 0;
    degs = 0;
    Opns = 0;
    opnCounts = 0;
    nConflicts = 0;
    mt = 0;
  };

 votingModel(const votingModel &toCopy): ROUND_CONST(toCopy.ROUND_CONST), n(toCopy.n), k(toCopy.k), collectionInterval(toCopy.collectionInterval), maxIter(toCopy.maxIter), a(toCopy.a), avgDeg(toCopy.avgDeg), id(toCopy.id), initDist(toCopy.initDist), conflicts(toCopy.conflicts), rewireTo(toCopy.rewireTo), fileName(toCopy.fileName) {
      A = new int*[toCopy.n];
      for(int i = 0; i < toCopy.n; i++) {
	  A[i] = new int[toCopy.n];
      }
      degs = new int[toCopy.n];
      Opns = new int[toCopy.n];
      nConflicts = new int[toCopy.n];
      opnCounts = new int[toCopy.k];
      for(int i = 0; i < toCopy.n; i++) {
	  degs[i] = toCopy.degs[i];
	  Opns[i] = toCopy.Opns[i];
	  nConflicts[i] = toCopy.nConflicts[i];
	  for(int j = 0; j < toCopy.n; j++) {
	      A[i][j] = toCopy.A[i][j];
	  }
      }
      for(int i = 0; i < toCopy.k; i++) {
	  opnCounts[i] = toCopy.opnCounts[i];
      }
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      mt = new std::mt19937(seed);
      rnNormalization = (double) (mt->max()+1);
  };

  votingModel &operator=(const votingModel &rhs) {
    if(this == &rhs) {
      return *this;
    }
    else {
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
      initDist = rhs.initDist;
      rewireTo.assign(rhs.rewireTo);
      fileName.assign(rhs.fileName);
      id = rhs.id;
      return *this;
    }
  };
};

#endif
