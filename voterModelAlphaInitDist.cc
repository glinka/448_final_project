#include <stdlib.h>
#include <string>
#include <sstream>
#include "vote.h"

using namespace std;

int main(int argc, char *argv[]) {
  /** 
      argc[1]: a/alpha
      argc[2]: k, number of opinions
      argc[3 to k+3]: initial distributions of 'k' opinions
      
      Constructs voter-model with variable alpha/init distr.
      with the following hard-coded parameters:
      n = 1000
      avgDeg = 4
      rewireTo = random
      maxIter = 100000
      collectionInterval = 50
      
      fileName is constructed as: fileName_n_avgDeg_a_initD[0]_initD[1]_...
  **/
  int n = atoi(argv[argc-2]);
  int maxIter = atoi(argv[argc-1]);
  double avgDeg = 4;
  string rewireTo = "same";
  int collectionInterval = 1000;
  double a = atof(argv[1]);
  int k = atoi(argv[2]);
  double initDist[k];
  stringstream ss;
  ss << "graphStats_" << n << "_" << avgDeg << "_" << a;
  for(int i = 0; i < k; i++) {
    initDist[i] = atof(argv[i+3]);
    ss << "_" << initDist[i];
  }
  ss << ".csv";
  string fileName = ss.str();
  votingModel model(n, k, maxIter, collectionInterval, a, avgDeg, initDist, rewireTo, fileName);
  model.vote();
  return 0;
}
