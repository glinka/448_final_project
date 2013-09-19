#ifndef VOTINGMODELCPI_H
#define VOTINGMODELCPI_H
#include <vector>
#include <iostream>

class votingModel;

typedef std::vector<int> vect;
typedef std::vector<vect> vmVects;
typedef std::vector<std::vector<int> > matrix;
typedef std::vector<matrix> vmMatrices;

class votingModelCPI {
 private:
  std::vector<votingModel> vms;
  const int projStep, waitingPeriod, collectionInterval, nMicroSteps;
  vect times;
  std::vector<vmVects > opns;
  std::vector<vmMatrices > adjMatrices;
  template <typename T>
    double project(const vect &times, const std::vector<T> &data);
  template <typename T>
    double average(const std::vector<T> &data);
  vect average(const std::vector<vect> &data);
  vmVects average(const std::vector<vmVects> &data);
  vmMatrices average(const std::vector<vmMatrices> &data);
  void saveData(const vmMatrices &data, std::ofstream &fileHandle);
  void saveData(const vmVects &data, std::ofstream &fileHandle);
  std::vector<double> findAvgdMinorityFractions(const std::vector<vmVects> &opns);
  vect findAvgdConflicts(const std::vector<vmMatrices> &As, const std::vector<vmVects> &opns);
 public:
  votingModelCPI(std::vector<votingModel> vms, const int projStep, int waitingPeriod, int collectionInterval, int nMicroSteps);
  ~votingModelCPI() {};
  void run(int nSteps);
};

#endif
