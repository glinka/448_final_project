#ifndef VOTE_H
#define VOTE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <Eigen/Sparse>
#include <chrono>
#include <iomanip>

int vote(const int n, const int k, const int maxIter, const int collectionInterval, double *initDist, const double a, const double avgDeg, const string rewireTo, const string fileName);

#endif
