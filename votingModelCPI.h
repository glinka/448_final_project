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
  const int waitingPeriod, collectionInterval, nMicroSteps;
  vect times;
  std::vector<vmVects > opns;
  std::vector<vmMatrices > adjMatrices;
  std::string file_header, file_name;
  template <typename T>
    double project(const vect &times, const std::vector<T> &data, const int proj_step);
  template <typename T>
    double average(const std::vector<T> &data);
  vect average(const std::vector<vect> &data);
  vmVects average(const std::vector<vmVects> &data);
  vmMatrices average(const std::vector<vmMatrices> &data);
  void saveData(const vect &data, std::ofstream &fileHandle);
  void saveData(const std::vector<vmMatrices> &data, std::ofstream &fileHandle);
  void saveData(const std::vector<vmVects> &data, std::ofstream &fileHandle);
  std::vector<double> findAvgdMinorityFractions(const std::vector<vmVects> &opns);
  vect findAvgdConflicts(const std::vector<vmMatrices> &As, const std::vector<vmVects> &opns);
  double getMinorityFrac(const vect &opns);
  int getConflicts(const matrix &A, const vect &opns);
 public:
  votingModelCPI(std::vector<votingModel> vms, int waitingPeriod, int collectionInterval, int nMicroSteps, std::string file_header, std::string file_name);
  ~votingModelCPI() {};
  int run(long int nSteps, int proj_step, int save_data_interval=100);
};

#endif
