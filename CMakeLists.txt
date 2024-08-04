cmake_minimum_required(VERSION 3.0.0)
project(crazyflie_ekf VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_FLAGS "-O3 -march=native -DNDEBUG -fopenmp")

add_executable(estimator src/main.cpp src/eskf.cpp src/csv_parser.cpp src/estimator.cpp)

set_property(TARGET estimator PROPERTY CXX_STANDARD 17)

target_include_directories(estimator PRIVATE ${PROJECT_SOURCE_DIR}/include)

#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#target_link_libraries (estimator Eigen3::Eigen)

set(EIGEN_ROOT_DIR "${PROJECT_SOURCE_DIR}/external/Eigen/eigen-3.4.0/")
target_include_directories(estimator PRIVATE ${EIGEN_ROOT_DIR})
