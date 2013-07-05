#include <stdlib.h>
#include <string>
#include <sstream>
#include "vote.h"

using namespace std;

int main(int argc, char *argv[]) {
  double avgDeg = 4;
  int collectionInterval = 1000;
  int n = atoi(argv[argc-3]);
  double a = atof(argv[1]);
  int k = atoi(argv[2]);
  double initDist[k];
  stringstream ss;
  ss << "graphStats_" << n << "_" << avgDeg << "_" << a;
  for(int i = 0; i < k; i++) {
    initDist[i] = atof(argv[i+3]);
    ss << "_" << initDist[i];
  }
  int maxIter = atoi(argv[argc-2]);
  string rewireTo = argv[argc-1];
  ss << ".csv";
  string fileName = ss.str();
  votingModel model(n, k, maxIter, collectionInterval, a, avgDeg, initDist, rewireTo, fileName);
  model.vote();
  return 0;
}
