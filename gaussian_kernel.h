#ifndef GAUSSIAN_KERNEL_FN_H
#define GAUSSIAN_KERNEL_FN_H
#include "kernel_function.h"
#include <vector>

template <typename T>
class Gaussian_Kernel : public Kernel_Function< std::vector<T> > {
 public:
  Gaussian_Kernel(const double epsilon);
  ~Gaussian_Kernel() {}
  virtual double kernel(const std::vector<T>& x1, const std::vector<T>& x2);
 private:
    const double epsilon_;
};

#include "gaussian_kernel.tpp"

#endif
