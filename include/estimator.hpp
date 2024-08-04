#pragma once
#include "csv_parser.hpp"
#include "eskf.hpp"
#include "measurement_types.hpp"
#include <Eigen/Core>

#include <cmath>

/**
 * @brief Just a helper class that simulates real-time measurements on embedded
 * device. Similar logic will have to be applied on the embbeded device.
 *
 */
class Estimator {
public:
  Estimator(const std::string &inputFilename,
            const std::string &outputFilename);
  void estimate();

private:
  enum COLS { T_IMU, AX, AY, AZ, GX, GY, GZ, T_OF, DX, DY, T_RANGE, R };

  void saveState(float timestamp);

  ESKF filter;
  std::ifstream inputFile;
  std::ofstream outputFile;

  //  TODO Do not duplicate constants as in eskf.hpp
  static constexpr float MICROSECOND_TO_SECOND = 1e-6f;
  static constexpr float GRAVITY = 9.812f;
  static constexpr float PI = 3.14159265358979323846f;
  static constexpr float DEG_TO_RAD = PI / 180.0f;
  static constexpr float MM_to_M = 1e-3f;
};