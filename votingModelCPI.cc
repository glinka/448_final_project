#include <vector>
#include <algorithm>
#include <iostream>
#include "votingModelCPI.h"

using namespace std;

votingModelCPI::votingModelCPI(double projectionInterval): projectionInterval(projectionInterval) {};

void votingModelCPI::collectData(const double data, const double interval) {
  dataPts.push_back(data);
  intervals.push_back(interval);
}

vector<double> *votingModelCPI::project() {
  sort(dataPts.begin(), dataPts.end());
  sort(intervals.begin(), intervals.end());
  double ymin = intervals.front();
  double ymax = intervals.back();
  double yOffset = (ymax - ymin)/2.0;
  double xmin = dataPts.front();
  double xmax = dataPts.back();
  double xOffset = (xmax - xmin)/2.0;
  double xnorm = 0;
  unsigned int i;
  for(i = 0; i < dataPts.size(); i++) {
    dataPts[i] -= xOffset;
    intervals[i] -= yOffset;
    xnorm += dataPts[i]*dataPts[i];
  }
  int innerProd = 0;
  for(i = 0; i < dataPts.size(); i++) {
    innerProd += dataPts[i]*intervals[i];
  }
  double a = innerProd / xnorm;
  double b = yOffset - a*xOffset;
  double xNew = xmin - projectionInterval;
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

