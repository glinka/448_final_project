#ifndef VOTE_H
#define VOTE_H
#include <string>

class votingModel {
 private:
  const int n, k, maxIter, collectionInterval;
  const double a, avgDeg;
  double *initDist;
  std::string rewireTo, fileName;
 public:
  int vote();
  votingModel(int n, int k, int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, std::string rewireTo, std::string fileName);
  ~votingModel(){};
};

#endif
