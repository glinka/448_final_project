#include <vector>
#include <algorithm>
#include <iostream>
#include "votingModelCPI.h"

using namespace std;

votingModelCPI::votingModelCPI(double projectionInterval): projectionInterval(projectionInterval) {
  int nData = 3;
  int i;
  vector<double> x;
  for(i = 0; i < nData; i++) {
    _data.push_back(x);
  }
};
void votingModelCPI::collectData(const vector<double> data) {
  for(int i = 0; i < data.size(); i++) {
    _data[i].push_back(data[i]);
  }
}

vector<double> votingModelCPI::project() {
  vector<double> minorityVsTime = fitLine(_data[2], _data[0]);
  vector<double> conflictsVsTime = fitLine(_data[2], _data[1]);
  double newTime = (_data[2]).back() + projectionInterval;
  double newMinorityFrac = minorityVsTime[0]*newTime + minorityVsTime[1];
  double newConflicts = conflictsVsTime[0]*newTime + conflictsVsTime[1];
  vector<double> newData;
  newData.push_back(newMinorityFrac);
  newData.push_back(newConflicts);
  return newData;
}

vector<double> votingModelCPI::fitLine(const vector<double> x, const vector<double> y) {
  //assume x.size() == y.size()
  int i;
  int n = x.size();
  double xSum = 0, ySum = 0, xSqSum = 0, xySum = 0;
  for(i = 0; i < n; i++) {
    xSum += x[i];
    ySum += y[i];
    xSqSum += (x[i]*x[i]);
    xySum += (x[i]*y[i]);
  }
  double a = (ySum*xSqSum - xSum*xySum)/(n*xSqSum - xSum*xSum);
  double b = (n*xySum - xSum*ySum)/(n*xSqSum - xSum*xSum);
  vector<double> ab;
  ab.push_back(a);
  ab.push_back(b);
  return ab;
}
