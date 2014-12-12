#ifndef CALCGRAPHPROPS_H
#define CALCGRAPHPROPS_H

class calcGraphProps {
 private:
  calcGraphProps() {};
  ~calcGraphProps() {};
 public:
  static int getTriangles(int **A, int *opns, const int n);
  static int get_conflict_cherries(const int *conflicts, const int total_conflicts, const int n);
  static int get_conflict_squares(int** A, const int* opns, const int* opnCounts, const int n);
  static int *getDegrees(int **A, const int n);
  static double *getAdjEigVals(int **A, const int n);
  static double **getAdjEigVects(int **A, const int n);
  static double *getLaplEigVals(int **A, const int n);
  static double **getLaplEigVects(int **A, const int n);
};

#endif
