#include "FusionEKF.h"

#include "Eigen/Dense"
#include "tools.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

const std::map<MeasurementPackage::SensorType, long> MeasurementPackage::SIZES{
    {MeasurementPackage::LASER, 2}, {MeasurementPackage::RADAR, 3}};

/**
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  // measurement covariance matrix - laser
  // clang-format off
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  // clang-format on

  noise_ax_ = 9.0;
  noise_ay_ = 9.0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF()
{}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  // Initialization

  if (!is_initialized_)
  {
    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // Initial uncertainty of the position, it does not
    // depend on which sensor type is used.
    MatrixXd P_initial(4, 4);
    // clang-format off
    P_initial << 100, 0, 0, 0,
                 0, 100, 0, 0,
                 0, 0, 100, 0,
                 0, 0, 0, 100;
    // clang-format on

    auto Q_initial = Tools::BuildNoiseMatrix(noise_ax_, noise_ay_, 0.05);

    // the initial transition matrix F_
    // implements a constant linear (!) motion model
    MatrixXd F_initial(4, 4);
    // clang-format off
    F_initial  << 1, 0, 1, 0,
                  0, 1, 0, 1,
                  0, 0, 1, 0,
                  0, 0, 0, 1;
    // clang-format on

    const auto data_size = measurement_pack.raw_measurements_.size();
    if (data_size != MeasurementPackage::SIZES.at(measurement_pack.sensor_type_))
    {
      cout << "Invalid measurement package (type " << measurement_pack.sensor_type_
           << "): data size wrong " << data_size << endl;
      return;
    }

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      VectorXd carthesian = Tools::ToCarthesianXY(measurement_pack.raw_measurements_);
      VectorXd x(4);
      x(0) = carthesian(0);
      x(1) = carthesian(1);
      x(2) = 0.0;
      x(3) = 0.0;

      ekf_.Init(x, P_initial, F_initial, H_laser_, R_laser_, Q_initial);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      VectorXd x(4);
      x(0) = measurement_pack.raw_measurements_(0);
      x(1) = measurement_pack.raw_measurements_(1);
      x(2) = 0.0;
      x(3) = 0.0;

      ekf_.Init(x, P_initial, F_initial, H_laser_, R_laser_, Q_initial);
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_     = true;
    return;
  }

  // Prediction

  const auto dt       = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  static constexpr auto dt_tolerance = 0.0001;
  if (abs(last_dt_ - dt) > dt_tolerance)
  {
    // Performing costly matri operations only if the DT was
    // changed, otherwise we still can re-use the former ones.
    // Set proper time gap in the state transition matrix F
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // Set the process covariance matrix Q
    ekf_.Q_ = Tools::BuildNoiseMatrix(noise_ax_, noise_ay_, dt);
  }
  last_dt_ = dt;

  ekf_.Predict();

  // Update

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    ekf_.H_ = Tools::CalculateRadarJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}
