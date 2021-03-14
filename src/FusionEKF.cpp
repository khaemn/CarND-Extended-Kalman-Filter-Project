#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;


const std::map<MeasurementPackage::SensorType, long> MeasurementPackage::SIZES {
    {MeasurementPackage::LASER, 2},
    {MeasurementPackage::RADAR, 3}
};

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_      = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // H_radar is Hj, should be recomputed on each step.

  noise_ax_ = 9.0;
  noise_ay_ = 9.0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // Initial uncertainty of the position, it does not
    // depend on which sensor type is used.
    MatrixXd P_initial(4, 4);
    P_initial << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1000, 0,
                 0, 0, 0, 1000;


    const auto data_size = measurement_pack.raw_measurements_.size();
    if (data_size != MeasurementPackage::SIZES.at(measurement_pack.sensor_type_)) {
      cout << "Invalid measurement package (type " << measurement_pack.sensor_type_
           << "): data size wrong " << data_size << endl;
      return;
    }
cout << "LOG data_size" << endl;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      VectorXd x(4);
      x(0) = measurement_pack.raw_measurements_(0);
      x(1) = measurement_pack.raw_measurements_(1);
      x(2) = 0.0;
      x(3) = 0.0;
cout << "LOG measurement_pack.raw_measurements_; " << endl;

      auto Q_initial = build_Q(noise_ax_, noise_ay_, 0.0);

      // the initial transition matrix F_
      MatrixXd F_initial(4, 4);
      F_initial  << 1, 0, 1, 0,
                    0, 1, 0, 1,
                    0, 0, 1, 0,
                    0, 0, 0, 1;

      ekf_.Init(x, P_initial, F_initial, H_laser_, R_laser_, Q_initial);
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  const auto dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Modify the F matrix so that the time is integrated
  // K is of shape [4,4], and the time is at positions (row, col):
  // [0, 2] and [1,3]
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Set the process covariance matrix Q
  ekf_.Q_ = build_Q(noise_ax_, noise_ay_, dt);

  Hj_ = tools.CalculateRadarJacobian(ekf_.x_);

cout << "LOG Predict" << endl;
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
cout << "LOG updates" << endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    return;
  } else {
    // TODO: Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }
cout << "LOG output" << endl;
  // print the output
//  cout << "x_ = " << ekf_.x_ << endl;
//  cout << "P_ = " << ekf_.P_ << endl;
}

// Builds the Q matrix using X and Y acceleration noise
// and a given delta-t
Eigen::MatrixXd FusionEKF::build_Q(double noise_ax, double noise_ay, double dt)
{

  // 2. Set the process covariance matrix Q
  const double dt_4_4 = pow(dt, 4.) / 4.;
  const double dt_3_2 = pow(dt, 3.) / 2.;
  const double dt_2_1 = pow(dt, 2.) / 1.;

  MatrixXd q = MatrixXd(4, 4);
  // clang-format off
  q <<  dt_4_4 * noise_ax, 0., dt_3_2 * noise_ax, 0.,
        0., dt_4_4 * noise_ay, 0., dt_3_2 * noise_ay,
        dt_3_2 * noise_ax, 0., dt_2_1 * noise_ax, 0.,
        0., dt_3_2 * noise_ay, 0., dt_2_1 * noise_ay;
  // clang-format on
  return q;
}
