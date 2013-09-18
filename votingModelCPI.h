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
  double project(const vect &data, const vect &times);
  vect average(const std::vector<vect> &data);
  vmVects average(const std::vector<vmVects> &data);
  vmMatrices average(const std::vector<vmMatrices> &data);
  void saveData(const std::vector<matrix> &data, std::ofstream &fileHandle);
  void saveData(const std::vector<vec> &data, std::ofstream &fileHandle);
 public:
  votingModelCPI(std::vector<votingModel> vms, const int projStep, int waitingPeriod, int collectionInterval, int nMicroSteps);
  ~votingModelCPI() {};
  void run(int nSteps);
};

#endif
