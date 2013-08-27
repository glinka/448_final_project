#ifndef VOTINGMODELCPI_H
#define VOTINGMODELCPI_H
#include <vector>

class votingModelCPI {
 private:
  std::vector<double> dataPts;
  std::vector<double> intervals;
  /** 
      currently contains:
      1. minority opn fraction
      2. number of conflicts
      3. number of iterations
  **/
  std::vector<std::vector<double>> _data;
  const double projectionInterval;
  double getLeastSquaresSlope(const vector<double> x, const vector<double> y);
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
