#include <vector>
#include <algorithm>
#include <iostream>
#include "votingModelCPI.h"

using namespace std;

votingModelCPI::votingModelCPI(double projectionInterval): projectionInterval(projectionInterval) {};

void votingModelCPI::collectData(const vector<double> data) {
  _data.push_back(data);
}

vector<double> *votingModelCPI::project() {
  sort((_data[0]).begin(), (_data[0]).end());
  sort((_data[1]).begin(), (_data[1]).end());
  sort((_data[2]).begin(), (_data[2]).end());
  double minorityVsTimeSlope = getLeastSquaresSlope(_data[2], _data[0]);
  double conflictsVsTimeSlope = getLeastSquaresSlope(_data[2], _data[1]);
  vector<double> *line = new vector<double>;
  double yNew;
  if(xNew > 0) {
    yNew = a*xNew + b;
    line->push_back(xNew);
  }
  else {
    yNew = ymin;
    line->push_back(xmin);
  }
  line->push_back(yNew);
  dataPts.clear();
  intervals.clear();
  return line;
}

