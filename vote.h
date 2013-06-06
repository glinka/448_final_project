#ifndef VOTE_H
#define VOTE_H
#include <string>
#include <Eigen/Sparse>

class votingModel {
 private:
    const double ROUND_CONST;
    const int n, k, maxIter, collectionInterval;
    const double a, avgDeg;
    double *initDist;
    Eigen::MatrixXi A, Opns;
    std::string rewireTo, fileName;
    void initGraph();
    int countConflicts();
    int graphConsistencyCheck();
 public:
    int vote();
    votingModel(int n, int k, int maxIter, int collectionInterval, double a, double avgDeg, double *initDist, std::string rewireTo, std::string fileName);
    ~votingModel(){};
};

#endif
