#include "vote.h"

using namespace std;

int main(int argc, char *argc[]) {
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
  int n = 1000;
  double avgDeg = 4;
  string rewireTo = "random";
  int maxIter = 100000;
  int collectionInterval = 50;
  double a = atof(argc[1]);
  int k = atoi(argc[2]);
  double initDist[k];
  stringstream ss;
  ss << "graphStats_" << n << "_" << avgDeg << "_" << a << "_";
  for(int i = 0; i < k; i++) {
    initDist[i] = atof(argc[i+3]);
    ss << "_" << initDist[i];
  }
  string fileName = ss.str();
  vote(n, k, maxIter, collectionInterval, initDist, a, avgDeg, rewireTo, fileName);
  return 0;
}
