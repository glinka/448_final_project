#include <stdlib.h>
#include <string>
#include <sstream>
#include "votingModel.h"

using namespace std;

int main(int argc, char *argv[]) {
  double avgDeg = 4;
  string rewireTo = argv[1];
  int n = atoi(argv[2]);
  int maxIter = atoi(argv[3]);
  int collectionInterval = atoi(argv[4]);
  double a = atof(argv[5]);
  int k = atoi(argv[6]);
  double initDist[k];
  stringstream ss;
  ss << "graphStats_" << n << "_" << avgDeg << "_" << a;
  for(int i = 0; i < k; i++) {
    initDist[i] = atof(argv[argc+i-k]);
    ss << "_" << initDist[i];
  }
  ss << ".csv";
  string fileName = ss.str();
  votingModel model(n, k, maxIter, collectionInterval, a, avgDeg, initDist, rewireTo, fileName);
  model.vote();
  return 0;
}
