#pragma once
#include "measurement_types.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>

/**
 * @brief Structure holding various ESKF filter parameters
 *
 */
struct FilterParams {
    // TODO: Restucture this, e.g., initial position can be given by initial state.position vector etc.
    // TODO: Rethink about xy, z separation

    /**
     * @brief Construct a new FilterParams object.
     *
     *
     *
     * @param stdDevInitialPosition_xy standard deviation of xy initial position
     * @param stdDevInitialPosition_z standard deviation of z initial position
     * @param stdDevInitialVelocity standard deviation of initial velocity (all directions)
     * @param stdDevInitialAttitude_rollpitch standard deviation of initial rotation in roll and pitch axes
     * @param stdDevInitialAttitude_yaw standard deviation of initial rotation in yaw axis
     * @param measNoiseAcc_xy standard deviation of measurement noise for accelerometer measurements in xy_direction (Q
     * matrix).
     * @param measNoiseAcc_z standard deviation of measurement noise for accelerometer measurements in z_direction (Q
     * matrix).
     * @param measNoiseGyro_rollpitch standard deviation of measurement noise for accelerometer measurements in
     * roll/pitch (Q matrix).
     * @param measNoiseGyro_yaw standard deviation of measurement noise for accelerometer measurements in yaw (Q matrix)
     * @param initialX initial x position
     * @param initialY initial y position
     * @param initialZ initial z position
     * @param initialYaw initial yaw
     */
    FilterParams(float stdDevInitialPosition_xy = 0, float stdDevInitialPosition_z = 0,
                 float stdDevInitialVelocity = 0.01, float stdDevInitialAttitude_rollpitch = 0.01,
                 float stdDevInitialAttitude_yaw = 0.01, float measNoiseAcc_xy = 0.5, float measNoiseAcc_z = 1.0,
                 float measNoiseGyro_rollpitch = 0.1, float measNoiseGyro_yaw = 0.1, float initialX = 0,
                 float initialY = 0, float initialZ = 0, float initialYaw = 0);

    float stdDevInitialPosition_xy;
    float stdDevInitialPosition_z;
    float stdDevInitialVelocity;
    float stdDevInitialAttitude_rollpitch;
    float stdDevInitialAttitude_yaw;

    float measNoiseAcc_xy;
    float measNoiseAcc_z;
    float measNoiseGyro_rollpitch;
    float measNoiseGyro_yaw;

    float initialX;
    float initialY;
    float initialZ;

    float initialYaw;
};

/**
 * @brief Enum for kalman filter states used in covariance indexing.
 *
 */
enum State_ids {
    KF_STATE_X,  // x position in world
    KF_STATE_Y,  // y position in world
    KF_STATE_Z,  // z position in world
    KF_STATE_VX, // x linear velocity in world
    KF_STATE_VY, // y linear velocity in world
    KF_STATE_VZ, // z linear velocity in world
    KF_STATE_D0, // error-state rotation
    KF_STATE_D1, // error-state rotation
    KF_STATE_D2, // error-state rotation
    KF_STATE_DIM // state dimension
};

/**
 * @brief Class encapsulating error-state kalman filter, its nominal and error-state states. Error-state calculcaitons
 * are inspired by https://arxiv.org/pdf/1711.02508.
 *
 */
class ESKF {
  private:
    /**
     * @brief Structure holding ESKF state.
     *
     */
    struct State {
        /**
         * @brief Construct a new State object.
         *
         * @param position position in world.
         * @param velocity velocity in world.
         * @param delta_q eror-state rotation vector.
         */
        State(const Eigen::Vector3f &position, const Eigen::Vector3f &velocity, const Eigen::Vector3f &delta_q)
            : position(position), velocity(velocity), delta_q(delta_q){};
        /**
         * @brief Construct a new State object
         *
         * @param init vector holding state values, indexed as in @see State_ids
         */
        State(const Eigen::Vector<float, KF_STATE_DIM> &init)
            : position(init.segment(0, 3)), velocity(init.segment(3, 3)), delta_q(init.segment(6, 3)){};

        Eigen::Vector3f position; // world frame position
        Eigen::Vector3f velocity; // body frame velocity
        Eigen::Vector3f delta_q;  // error-state rotation
    };

  public:
    /**
     * @brief Construct a new ESKF object.
     *
     * @param params filter parameters.
     */
    ESKF(const FilterParams &params);

    /**
     * @brief ESKF prediction function. Propagates both nominal states (position and velocty) and rotation error-state.
     *
     * @param acc accelerometer measurement.
     * @param gyro gyro measurement.
     * @param dt delta_t until last prediction.
     * @param quadIsFlying is quad flying or moved with hands.
     */
    void predict(const Eigen::Vector3f &acc, const Eigen::Vector3f &gyro, float dt, bool quadIsFlying);

    /**
     * @brief Update function for time-of-flight measurement.
     *
     * @param tof time-of-flight measurement.
     */
    void updateTof(const TofMeasurement &tof);

    /**
     * @brief Updata function for optical-flow measurement.
     *
     * @param flow optical flow measurement.
     * @param gyroLatest latesy gyro measurement.
     */
    void updateFlow(const FlowMeasurement &flow, const Eigen::Vector3f &gyroLatest);

    /**
     * @brief Reset the ESKF, i.e., update the error-state back into nominal state. Usually done after update with
     * external measurement.
     *
     */
    void reset();

    // We can return matrices (not references) due to copy ellision
    inline Eigen::Vector3f getPosition() { return mState.position; };
    inline Eigen::Vector3f getVelocity() { return mState.velocity; };
    inline Eigen::Matrix3f getRotationMatrix() { return mRotationMatrix; };
    inline Eigen::Vector3f getRollPitchYaw() { return mRotationMatrix.eulerAngles(2, 1, 0); };

  private:
    /**
     * @brief Skew-symmetric matrix calculation.
     *
     * @param v 3D vector.
     * @return  Skew-symmetric matrix.
     */
    static Eigen::Matrix3f skew(const Eigen::Vector3f &v);

    /**
     * @brief Correction function called during state update.
     *
     * @param correction correction vector
     */
    void correctState(const Eigen::Vector<float, KF_STATE_DIM> &correction);

    /**
     * @brief Adds process noise to current covariance. This is currently a very simple implementation, needs additional
     * considerations.
     * @todo Implement proper cross-correlation between states
     * @todo Implement rotation variant process noise, due to different xy and z deviations
     * @param dt
     */
    void addProcessNoise(float dt);

    /**
     * @brief Updates the ESKF state with scalar measurement.
     *
     * @param H measurement jacobian
     * @param error error between state and measurement
     * @param stdMeasNoise measurement noise stdandard deviation
     */
    void scalarUpdate(const Eigen::Matrix<float, 1, KF_STATE_DIM> &H, float error, float stdMeasNoise);

    /**
     * @brief Initialie EKF state.
     *
     * @param params filter params
     * @return initialized state as N-dimensional vector
     */
    Eigen::Vector<float, KF_STATE_DIM> initState(const FilterParams &params);

    /**
     * @brief Initialize ESKF covariance
     *
     * @param params filter params
     * @return covariance as NxN matrix
     */
    Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM> initCovariance(const FilterParams &params);

    /**
     * @brief Initialize rotation. Needed due to error-state formulation
     *
     * @param params filter params
     * @return rotation as quaternion
     */
    Eigen::Quaternionf initRotation(const FilterParams &params);

    Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM> mCovariance;
    State mState;

    Eigen::Quaternionf mQuaternion;
    Eigen::Quaternionf mInitialQuaternion;
    Eigen::Matrix3f
        mRotationMatrix; // we also keep the rotation matrix in memory to avoid unecessary quat -> rot_mat conversions

    FilterParams params;
    bool isUpdated;

    const Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM> I_9 =
        Eigen::Matrix<float, KF_STATE_DIM, KF_STATE_DIM>::Identity();
    const Eigen::Matrix<float, 3, 3> I_3 = Eigen::Matrix<float, 3, 3>::Identity();
    const Eigen::Vector3f GRAVITY_VECTOR = {0, 0, GRAVITY_MAGNITUDE};

    static constexpr float GRAVITY_MAGNITUDE = 9.81f;
    static constexpr float PI = 3.14159265358979323846f;
    static constexpr float DEG_TO_RAD_ = PI / 180.0f;
    static constexpr float ROLLPITCH_ZERO_REVERSION = 0.001f;
    static constexpr float MAX_COVARIANCE = 100.f;
    static constexpr float MIN_COVARIANCE = 1e-6f;
};