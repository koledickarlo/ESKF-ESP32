#include "estimator.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>

int main(int argc, char *argv[]) {
  Estimator estimator{arv[1], argv[2]};
  estimator.estimate();

  return 0;
}