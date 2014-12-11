#include <cmath>
#include "gaussian_kernel.h"

template <typename T>
Gaussian_Kernel<T>::Gaussian_Kernel(const double epsilon): epsilon_(epsilon) {}

template <typename T>
double Gaussian_Kernel<T>::kernel(const std::vector<T>& x1, const std::vector<T>& x2) {
  // if(x1.size() != x2.size()) {
  //   std::cout << "array dimensions do not match" << std::endl;
  //   exit(1);
  // }
  int n = x1.size();
  double norm = 0;
  for(int i = 0; i < n; i++) {
    norm += std::pow(x1[i] - x2[i], 2);
  }
  return std::exp(-norm/epsilon_);
}
