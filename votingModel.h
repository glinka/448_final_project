#ifndef VOTINGMODEL_H
#define VOTINGMODEL_H
#include <string>
#include <vector>
#include <random>

class votingModelCPI;

class votingModel {
 private:
  const bool project;
  const double ROUND_CONST;
  const int n, k, maxIter, collectionInterval;
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
  int consistencyCheck();
  double genURN();
 public:
  int vote();
  void step();
  void initGraph(double *dist, int conflicts);
  template <typename dataType>
    void saveData(const std::vector<dataType> data, std::ofstream &fileHandle);
  template <typename dataType>
    void saveData(const std::vector<std::vector<dataType> > data, std::ofstream &fileHandle);
  votingModel(int n, int k, int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, std::string rewireTo, std::string fileName, bool project);

  std::vector<int> getOpns() {
    std::vector<int> newOpns;
    for(int i = 0; i < n; i++) {
      newOpns.push_back(Opns[i]);
    }
    return newOpns;
  }
  std::vector<std::vector<int> > getAdjMatrix() {
    std::vector<std::vector<int> > newA;
    std::vector<int> v;
    for(int i = 0; i < n; i++) {
      newA.push_back(v);
      for(int j = 0; i < n; j++) {
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
  };
};

#endif
