#include "estimator.hpp"

Estimator::Estimator(const std::string &inputFilename,
                     const std::string &outputFilename)
    : inputFile(inputFilename), outputFile(outputFilename),
      filter(FilterParams{}) {}

void Estimator::estimate() {
  // Skip first line
  CsvParser row{};
  inputFile >> row;
  inputFile >> row;
  uint32_t last_timestamp = std::stoul(row[T_IMU]);
  uint32_t last_flow_timestamp = last_timestamp;
  Eigen::Vector3f gyro;

  while (inputFile >> row) {
    uint32_t t_imu = std::stoul(row[T_IMU]);
    if (t_imu) {
      Eigen::Vector3f acc{std::stof(row[AX]), std::stof(row[AY]),
                          std::stof(row[AZ])};
      gyro = {std::stof(row[GX]), std::stof(row[GY]), std::stof(row[GZ])};

      filter.predict(acc * GRAVITY, gyro * DEG_TO_RAD,
                     (t_imu - last_timestamp) * MICROSECOND_TO_SECOND, false);
      saveState(t_imu * MICROSECOND_TO_SECOND);
      last_timestamp = t_imu;
    }
    uint32_t t_of = std::stoul(row[T_OF]);
    if (t_of) {
      FlowMeasurement flow{std::stoi(row[DX]), std::stoi(row[DY]),
                           (t_of - last_flow_timestamp) *
                               MICROSECOND_TO_SECOND};
      last_flow_timestamp = t_of;
      filter.updateFlow(flow, gyro);
    }

    uint32_t t_range = std::stoul(row[T_RANGE]);
    if (t_range) {
      TofMeasurement tof{std::stof(row[R]) * MM_to_M};
      filter.updateTof(tof);
    }
    filter.reset();
  }
}

void Estimator::saveState(float timestamp) {
  auto position = filter.getPosition();
  auto rotation = filter.getRollPitchYaw();
  auto velocity = filter.getRotationMatrix() * filter.getVelocity();

  outputFile << timestamp << " ";
  outputFile << position.transpose() << " ";
  outputFile << velocity.transpose() << " ";
  outputFile << rotation;
  outputFile << std::endl;
}
