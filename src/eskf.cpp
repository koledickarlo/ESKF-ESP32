#include "eskf.hpp"
#include <chrono>
#include <cmath>
#include <iostream>

constexpr float EPS = 1e-6f;

FilterParams::FilterParams(float stdDevInitialPosition_xy, float stdDevInitialPosition_z, float stdDevInitialVelocity,
                           float stdDevInitialAttitude_rollpitch, float stdDevInitialAttitude_yaw,
                           float measNoiseAcc_xy, float measNoiseAcc_z, float measNoiseGyro_rollpitch,
                           float measNoiseGyro_yaw, float initialX, float initialY, float initialZ, float initialYaw)
    : stdDevInitialPosition_xy(stdDevInitialPosition_xy), stdDevInitialPosition_z(stdDevInitialPosition_z),
      stdDevInitialVelocity(stdDevInitialVelocity), stdDevInitialAttitude_rollpitch(stdDevInitialAttitude_rollpitch),
      stdDevInitialAttitude_yaw(stdDevInitialAttitude_yaw), measNoiseAcc_xy(measNoiseAcc_xy),
      measNoiseAcc_z(measNoiseAcc_z), measNoiseGyro_rollpitch(measNoiseGyro_rollpitch),
      measNoiseGyro_yaw(measNoiseGyro_yaw), initialX(initialX), initialY(initialY), initialZ(initialZ),
      initialYaw(initialYaw){};

ESKF::ESKF(const FilterParams &params)
    : mState(initState(params)), mCovariance(initCovariance(params)), mQuaternion(initRotation(params)),
      mInitialQuaternion(mQuaternion), mRotationMatrix(mQuaternion.toRotationMatrix()), isUpdated(false),
      params(params){};

void ESKF::predict(const Eigen::Vector3f &acc, const Eigen::Vector3f &gyro, float dt, bool quadIsFlying) {
    Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM> A = Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM>::Zero();

    float dt2 = dt * dt;

    A.block<3, 3>(KF_STATE_X, KF_STATE_X).diagonal() = Eigen::Vector3f::Ones();
    A.block<3, 3>(KF_STATE_VX, KF_STATE_VX).diagonal() = Eigen::Vector3f::Ones();
    A.block<3, 3>(KF_STATE_D0, KF_STATE_D0).diagonal() = Eigen::Vector3f::Ones();

    // position from body-frame velocity
    A.block<3, 3>(KF_STATE_X, KF_STATE_VX) = mRotationMatrix * dt;

    // position from attitude error
    A.block<3, 3>(KF_STATE_X, KF_STATE_D0) = -mRotationMatrix * skew(mState.velocity) * dt;

    // body-frame velocity from body-frame velocity
    A.block<3, 3>(KF_STATE_VX, KF_STATE_VX) = I_3 - skew(gyro) * dt;

    // body-frame velocity from attitude error
    A(KF_STATE_VX, KF_STATE_D0) = 0;
    A(KF_STATE_VY, KF_STATE_D0) = -GRAVITY_MAGNITUDE * mRotationMatrix(2, 2) * dt;
    A(KF_STATE_VZ, KF_STATE_D0) = GRAVITY_MAGNITUDE * mRotationMatrix(2, 1) * dt;

    A(KF_STATE_VX, KF_STATE_D1) = GRAVITY_MAGNITUDE * mRotationMatrix(2, 2) * dt;
    A(KF_STATE_VY, KF_STATE_D1) = 0;
    A(KF_STATE_VZ, KF_STATE_D1) = -GRAVITY_MAGNITUDE * mRotationMatrix(2, 0) * dt;

    A(KF_STATE_VX, KF_STATE_D2) = -GRAVITY_MAGNITUDE * mRotationMatrix(2, 1) * dt;
    A(KF_STATE_VY, KF_STATE_D2) = GRAVITY_MAGNITUDE * mRotationMatrix(2, 0) * dt;
    A(KF_STATE_VZ, KF_STATE_D2) = 0;

    // attitude error from attitude error
    /**
     * As derived in "Covariance Correction Step for Kalman Filtering with an Attitude"
     * http://arc.aiaa.org/doi/abs/10.2514/1.G000848
     */
    float d0 = gyro(0) * dt / 2;
    float d1 = gyro(1) * dt / 2;
    float d2 = gyro(2) * dt / 2;

    A(KF_STATE_D0, KF_STATE_D0) = 1 - d1 * d1 / 2 - d2 * d2 / 2;
    A(KF_STATE_D0, KF_STATE_D1) = d2 + d0 * d1 / 2;
    A(KF_STATE_D0, KF_STATE_D2) = -d1 + d0 * d2 / 2;

    A(KF_STATE_D1, KF_STATE_D0) = -d2 + d0 * d1 / 2;
    A(KF_STATE_D1, KF_STATE_D1) = 1 - d0 * d0 / 2 - d2 * d2 / 2;
    A(KF_STATE_D1, KF_STATE_D2) = d0 + d1 * d2 / 2;

    A(KF_STATE_D2, KF_STATE_D0) = d1 + d0 * d2 / 2;
    A(KF_STATE_D2, KF_STATE_D1) = -d0 + d1 * d2 / 2;
    A(KF_STATE_D2, KF_STATE_D2) = 1 - d0 * d0 / 2 - d1 * d1 / 2;

    mCovariance = A * mCovariance * A.transpose();

    if (quadIsFlying) // only acceleration in z direction
    {
        Eigen::Vector3f az{0, 0, acc(2)};
        // Position propagation, velocity and acceleration are rotated to the world frame
        mState.position += mRotationMatrix * (mState.velocity * dt + az * dt2 / 2.0f) - GRAVITY_VECTOR * dt2 / 2.0f;
        // Velocity propagation
        mState.velocity += (az - skew(gyro) * mState.velocity - mRotationMatrix.transpose() * GRAVITY_VECTOR) * dt;
    } else {
        // Position propagation, velocity and acceleration are rotated to the world frame
        mState.position += mRotationMatrix * (mState.velocity * dt + acc * dt2 / 2.0f) - GRAVITY_VECTOR * dt2 / 2.0f;
        // Velocity propagation
        mState.velocity += (acc - skew(gyro) * mState.velocity - mRotationMatrix.transpose() * GRAVITY_VECTOR) * dt;
    }

    const Eigen::Vector3f dt_gyro = dt * gyro;
    const Eigen::Quaternionf delta_q{Eigen::AngleAxisf{dt_gyro.norm(), dt_gyro.normalized()}};

    // rotate the quaternion by infitesimal gyro rotation
    mQuaternion *= delta_q;

    if (!quadIsFlying) {
        const float keep = 1.0f - ROLLPITCH_ZERO_REVERSION;
        mQuaternion.coeffs() = mQuaternion.coeffs() * keep + ROLLPITCH_ZERO_REVERSION * mInitialQuaternion.coeffs();
    }
    mQuaternion.normalize();

    addProcessNoise(dt);
    isUpdated = true;
}

Eigen::Matrix3f ESKF::skew(const Eigen::Vector3f &v) {
    Eigen::Matrix3f skew;

    skew << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;
    return skew;
}

void ESKF::correctState(const Eigen::Vector<float, KF_STATE_DIM> &correction) {
    mState.position += correction.segment(KF_STATE_X, 3);
    mState.velocity += correction.segment(KF_STATE_VX, 3);
    mState.delta_q += correction.segment(KF_STATE_D0, 3);
}

void ESKF::addProcessNoise(float dt) {
    Eigen::Vector<float, KF_STATE_DIM> Q_diag{
        std::pow(params.measNoiseAcc_xy * dt * dt, 2),
        std::pow(params.measNoiseAcc_xy * dt * dt, 2),
        std::pow(params.measNoiseAcc_z * dt * dt, 2),
        std::pow(params.measNoiseAcc_xy * dt, 2),
        std::pow(params.measNoiseAcc_xy * dt, 2),
        std::pow(params.measNoiseAcc_z * dt, 2),
        std::pow(params.measNoiseGyro_rollpitch * dt, 2),
        std::pow(params.measNoiseGyro_rollpitch * dt, 2),
        std::pow(params.measNoiseGyro_yaw * dt, 2)};

    mCovariance.diagonal() += Q_diag;

    mCovariance = 0.5 * (mCovariance + mCovariance.transpose());
    for (int i = 0; i < KF_STATE_DIM; i++) {
        for (int j = i; j < KF_STATE_DIM; j++) {
            float p = mCovariance(i, j);
            if (std::isnan(p) || p > MAX_COVARIANCE) {
                mCovariance(i, j) = mCovariance(j, i) = MAX_COVARIANCE;
            } else if (i == j && p < MIN_COVARIANCE) {
                mCovariance(i, j) = mCovariance(j, i) = MIN_COVARIANCE;
            } else {
                mCovariance(i, j) = mCovariance(j, i) = p;
            }
        }
    }
}

void ESKF::scalarUpdate(const Eigen::Matrix<float, 1, KF_STATE_DIM> &H, float error, float stdMeasNoise) {
    float R = stdMeasNoise * stdMeasNoise;
    float S = H * mCovariance * H.transpose() + R;
    const Eigen::Matrix<float, KF_STATE_DIM, 1> K = mCovariance * H.transpose() / S;
    const Eigen::Vector<float, KF_STATE_DIM> state_correction = K * error;

    correctState(state_correction);

    const Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM> KH_I = K * H - I_9;

    mCovariance = KH_I * mCovariance * KH_I.transpose() + K * R * K.transpose();
    mCovariance = 0.5 * (mCovariance + mCovariance.transpose());
    for (int i = 0; i < KF_STATE_DIM; i++) {
        for (int j = i; j < KF_STATE_DIM; j++) {
            float p = mCovariance(i, j);
            if (std::isnan(p) || p > MAX_COVARIANCE) {
                mCovariance(i, j) = mCovariance(j, i) = MAX_COVARIANCE;
            } else if (i == j && p < MIN_COVARIANCE) {
                mCovariance(i, j) = mCovariance(j, i) = MIN_COVARIANCE;
            } else {
                mCovariance(i, j) = mCovariance(j, i) = p;
            }
        }
    }
    isUpdated = true;
}

Eigen::Quaternionf ESKF::initRotation(const FilterParams &params) {
    return Eigen::Quaternionf{std::cos(params.initialYaw / 2), 0, 0, std::sin(params.initialYaw / 2)};
}

Eigen::Vector<float, KF_STATE_DIM> ESKF::initState(const FilterParams &params) {
    return Eigen::Vector<float, KF_STATE_DIM>{params.initialX, params.initialY, params.initialZ, 0, 0, 0, 0, 0, 0};
}

Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM> ESKF::initCovariance(const FilterParams &params) {
    return Eigen::Vector<float, KF_STATE_DIM>{
        std::pow(params.stdDevInitialPosition_xy, 2),        std::pow(params.stdDevInitialPosition_xy, 2),
        std::pow(params.stdDevInitialPosition_z, 2),         std::pow(params.stdDevInitialVelocity, 2),
        std::pow(params.stdDevInitialVelocity, 2),           std::pow(params.stdDevInitialVelocity, 2),
        std::pow(params.stdDevInitialAttitude_rollpitch, 2), std::pow(params.stdDevInitialAttitude_rollpitch, 2),
        std::pow(params.stdDevInitialAttitude_yaw, 2)}
        .asDiagonal();
}

void ESKF::updateTof(const TofMeasurement &tof) {
    Eigen::Matrix<float, 1, KF_STATE_DIM> H = Eigen::Matrix<float, 1, KF_STATE_DIM>::Zero();

    if (std::abs(mRotationMatrix(2, 2)) > 0.1 && mRotationMatrix(2, 2) > 0) {
        float angle = std::abs(std::acos(mRotationMatrix(2, 2)) - DEG_TO_RAD_ * 7.5);
        if (angle < 0) {
            angle = 0;
        }
        float predictedDistance = mState.position(2) / std::cos(angle);
        float measuredDistance = tof.range;

        H(0, KF_STATE_Z) = 1 / std::cos(angle);

        scalarUpdate(H, measuredDistance - predictedDistance, tof.stdDev);
    }
}

void ESKF::updateFlow(const FlowMeasurement &flow, const Eigen::Vector3f &gyroLatest) {
    Eigen::Matrix<float, 1, KF_STATE_DIM> H = Eigen::Matrix<float, 1, KF_STATE_DIM>::Zero();

    // TODO: move this to sensor definition
    constexpr float Npix = 35.0f;
    // 2*sin(42/2); 42degree is the agnle of aperture, here we computed the corresponding ground length
    constexpr float thetapix = 0.71674f;
    constexpr float FLOW_RESOLUTION = 0.1f;

    const float omegax_b = gyroLatest(0) * DEG_TO_RAD_;
    const float omegay_b = gyroLatest(1) * DEG_TO_RAD_;

    const float dx_g = mState.velocity(0);
    const float dy_g = mState.velocity(1);

    float z_g;
    if (mState.position(2) < 0.1f) {
        z_g = 0.1;
    } else {
        z_g = mState.position(2);
    }

    const float predictedNX = (flow.dt * Npix / thetapix) * ((dx_g * mRotationMatrix(2, 2) / z_g) - omegay_b);
    const float measuredNX = flow.dx * FLOW_RESOLUTION;

    H(0, KF_STATE_Z) = (Npix * flow.dt / thetapix) * ((mRotationMatrix(2, 2) * dx_g) / (-z_g * z_g));
    H(0, KF_STATE_VX) = (Npix * flow.dt / thetapix) * (mRotationMatrix(2, 2) / z_g);

    // X update
    scalarUpdate(H, (measuredNX - predictedNX), flow.stdDevx * FLOW_RESOLUTION);

    H = Eigen::Matrix<float, 1, KF_STATE_DIM>::Zero();

    const float predictedNY = (flow.dt * Npix / thetapix) * ((dy_g * mRotationMatrix(2, 2) / z_g) + omegax_b);
    const float measuredNY = flow.dy * FLOW_RESOLUTION;

    H(0, KF_STATE_Z) = (Npix * flow.dt / thetapix) * ((mRotationMatrix(2, 2) * dy_g) / (-z_g * z_g));
    H(0, KF_STATE_VY) = (Npix * flow.dt / thetapix) * (mRotationMatrix(2, 2) / z_g);

    // Y update
    scalarUpdate(H, (measuredNY - predictedNY), flow.stdDevy * FLOW_RESOLUTION);
}

void ESKF::reset() {
    if (!isUpdated) {
        return;
    }

    // Incorporate the attitude error (Kalman filter state) with the attitude
    float v0 = mState.delta_q(0);
    float v1 = mState.delta_q(1);
    float v2 = mState.delta_q(2);

    // Move attitude error into attitude if any of the angle errors are large enough
    if ((fabsf(v0) > 0.1e-3f || fabsf(v1) > 0.1e-3f || fabsf(v2) > 0.1e-3f) &&
        (fabsf(v0) < 10 && fabsf(v1) < 10 && fabsf(v2) < 10)) {

        const Eigen::Quaternionf delta_q{Eigen::AngleAxisf{mState.delta_q.norm(), mState.delta_q.normalized()}};

        // rotate the quaternion by infitesimal gyro rotation
        mQuaternion *= delta_q;
        mQuaternion.normalize();

        float d0 = v0 / 2; // the attitude error vector (v0,v1,v2) is small,
        float d1 = v1 / 2; // so we use a first order approximation to d0 = tan(|v0|/2)*v0/|v0|
        float d2 = v2 / 2;

        Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM> A = Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM>::Zero();

        A.block<3, 3>(KF_STATE_X, KF_STATE_X).diagonal() = Eigen::Vector3f::Ones();
        A.block<3, 3>(KF_STATE_VX, KF_STATE_VX).diagonal() = Eigen::Vector3f::Ones();

        /** Rotate the covariance, since we've rotated the body
         * As derived in "Covariance Correction Step for Kalman Filtering with an Attitude"
         * http://arc.aiaa.org/doi/abs/10.2514/1.G000848
         */

        A(KF_STATE_D0, KF_STATE_D0) = 1 - d1 * d1 / 2 - d2 * d2 / 2;
        A(KF_STATE_D0, KF_STATE_D1) = d2 + d0 * d1 / 2;
        A(KF_STATE_D0, KF_STATE_D2) = -d1 + d0 * d2 / 2;

        A(KF_STATE_D1, KF_STATE_D0) = -d2 + d0 * d1 / 2;
        A(KF_STATE_D1, KF_STATE_D1) = 1 - d0 * d0 / 2 - d2 * d2 / 2;
        A(KF_STATE_D1, KF_STATE_D2) = d0 + d1 * d2 / 2;

        A(KF_STATE_D2, KF_STATE_D0) = d1 + d0 * d2 / 2;
        A(KF_STATE_D2, KF_STATE_D1) = -d0 + d1 * d2 / 2;
        A(KF_STATE_D2, KF_STATE_D2) = 1 - d0 * d0 / 2 - d1 * d1 / 2;

        mCovariance = A * mCovariance * A.transpose();
    }

    // convert the new attitude to a rotation matrix, such that we can rotate body-frame velocity and acc
    mRotationMatrix = mQuaternion.toRotationMatrix();
    // reset the attitude error

    mState.delta_q = Eigen::Vector3f::Zero();

    // enforce symmetry of the covariance matrix, and ensure the values stay bounded
    mCovariance = 0.5 * (mCovariance + mCovariance.transpose());
    for (int i = 0; i < KF_STATE_DIM; i++) {
        for (int j = i; j < KF_STATE_DIM; j++) {
            float p = mCovariance(i, j);
            if (std::isnan(p) || p > MAX_COVARIANCE) {
                mCovariance(i, j) = mCovariance(j, i) = MAX_COVARIANCE;
            } else if (i == j && p < MIN_COVARIANCE) {
                mCovariance(i, j) = mCovariance(j, i) = MIN_COVARIANCE;
            } else {
                mCovariance(i, j) = mCovariance(j, i) = p;
            }
        }
    }
    isUpdated = false;
}
