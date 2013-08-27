#ifndef VOTINGMODELCPI_H
#define VOTINGMODELCPI_H
#include <vector>

class votingModelCPI {
 private:
  std::vector<double> dataPts;
  std::vector<double> intervals;
  std::vector<std::vector<double>> _data;
  const double projectionInterval;
 public:
  std::vector<double> *project();
  void collectData(const std::vector<double> data);
  double getProjectionstep() {
    return projectionInterval;
  };
  votingModelCPI(double projectionInterval);
  ~votingModelCPI() {};
};

#endif
