#pragma once
#include <cmath>
#include <string>
#include <vector>

struct TofMeasurement {
    TofMeasurement(float range) : range(range) {
        stdDev = expStdA * (1.0f + std::exp(ExpCoeff() * (range - expPointA)));
    };
    float range;
    float stdDev;

  private:
    static constexpr float ExpCoeff() { return std::log(expStdB / expStdA) / (expPointB - expPointA); };
    static constexpr float expPointA = 1.0f;
    static constexpr float expStdA = 0.0025f; // STD at elevation expPointA [m]
    static constexpr float expPointB = 1.3f;
    static constexpr float expStdB = 0.2f; // STD at elevation expPointB [m]
};

struct FlowMeasurement {
    FlowMeasurement(int16_t dx, int16_t dy, float dt) : dx(dx), dy(dy), dt(dt) {}
    int16_t dx;
    int16_t dy;
    float dt;

    static constexpr float stdDevx = 2.0f;
    static constexpr float stdDevy = 2.0f;
};
