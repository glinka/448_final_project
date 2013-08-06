#ifndef VOTINGMODELCPI_H
#define VOTINGMODELCPI_H
#include <vector>

class votingModelCPI {
 private:
  std::vector<double> dataPts;
  std::vector<double> intervals;
  const double projectionInterval;
 public:
  std::vector<double> *project();
  void collectData(const double data, const double interval);
  double getProjectionstep() {
    return projectionInterval;
  };
  votingModelCPI(double projectionInterval);
  ~votingModelCPI() {};
};

#endif
