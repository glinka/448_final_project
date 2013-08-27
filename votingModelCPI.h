#ifndef VOTINGMODELCPI_H
#define VOTINGMODELCPI_H
#include <vector>

class votingModelCPI {
 private:
  /** 
      currently contains:
      1. minority opn fraction
      2. number of conflicts
      3. number of iterations
  **/
  std::vector<std::vector<double>> _data;
  const double projectionInterval;
  std::vector<double> fitLine(const std::vector<double> x, const std::vector<double> y);
 public:
  std::vector<double> project();
  double getProjectionStep() {
    return projectionInterval;
  }
  void collectData(const std::vector<double> data);
  votingModelCPI(double projectionInterval);
  ~votingModelCPI() {};
};

#endif
