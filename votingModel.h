#ifndef VOTINGMODEL_H
#define VOTINGMODEL_H
#include <string>

class votingModel {
 private:
    const double ROUND_CONST;
    const int n, k, maxIter, collectionInterval;
    const double a, avgDeg;
    double *initDist;
    int *degs;
    int *Opns;
    int **A;
    std::string rewireTo, fileName;
    void initGraph();
    int consistencyCheck();
 public:
    int vote();
    votingModel(int n, int k, int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, std::string rewireTo, std::string fileName);
    ~votingModel() {
	for(int i = 0; i < n; i++) {
	  delete[] A[i];
	}
	delete[] A;
	delete[] degs;
	delete[] Opns;
    };
};

#endif
