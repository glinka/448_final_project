#ifndef VOTINGMODEL_H
#define VOTINGMODEL_H
#include <string>

class votingModelCPI;

class votingModel {
 private:
  const bool project;
  const double ROUND_CONST;
  const int n, k, maxIter, collectionInterval;
  const double a, avgDeg;
  double *initDist;
  int *degs;
  int *Opns;
  int **A;
  //tracks the number of conflicting edges at each vertex
  int *nConflicts;
  std::string rewireTo, fileName;
  votingModelCPI *vmCPI;
  void initGraph(double *dist, int conflicts, bool firstInitialization);
  int consistencyCheck();
 public:
  int vote();
  votingModel(int n, int k, int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, std::string rewireTo, std::string fileName, bool project, votingModelCPI *CPI);
  ~votingModel() {
    for(int i = 0; i < n; i++) {
      delete[] A[i];
    }
    delete[] A;
    delete[] degs;
    delete[] Opns;
    delete[] nConflicts;
  };
};

#endif
