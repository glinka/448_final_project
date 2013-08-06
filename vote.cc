#include <stdlib.h>
#include <string>
#include <sstream>
#include "votingModel.h"

using namespace std;

int main(int argc, char *argv[]) {
  double avgDeg = 4;
  int offset = 1;
  bool project = atoi(argv[1]);
  string rewireTo = argv[1+offset];
  int n = atoi(argv[2+offset]);
  int maxIter = atoi(argv[3+offset]);
  int collectionInterval = atoi(argv[4+offset]);
  double a = atof(argv[5+offset]);
  int k = atoi(argv[6+offset]);
  double initDist[k];
  stringstream ss;
  ss << "graphStats_" << n << "_" << avgDeg << "_" << a;
  for(int i = 0; i < k; i++) {
    initDist[i] = atof(argv[argc+i-k+offset]);
    ss << "_" << initDist[i];
  }
  for(int i = 0; i < argc; i++) {
    cout << argv[i] << endl;
  }
  ss << ".csv";
  string fileName = ss.str();
  votingModel model(n, k, maxIter, collectionInterval, a, avgDeg, initDist, rewireTo, fileName, project);
  model.vote();
  return 0;
}
