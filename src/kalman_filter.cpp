#include "kalman_filter.h"

#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter()
{}

KalmanFilter::~KalmanFilter()
{}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &H_in,
                        MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  MatrixXd Ft = F_.transpose();

  x_ = F_ * x_;
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  VectorXd y = z - (H_ * x_);

  UpdateImpl(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  VectorXd y = z - Tools::ToPolar(x_);

  // Normalize y(1) (e.g. phi) to [-pi .. pi] range
  y(1) = atan2(sin(y(1)), cos(y(1)));

  UpdateImpl(y);
}

void KalmanFilter::UpdateImpl(const Eigen::VectorXd &y)
{
  MatrixXd Ht = H_.transpose();
  MatrixXd S  = (H_ * P_ * Ht) + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K  = P_ * Ht * Si;

  // new estimate
  x_              = x_ + (K * y);
  long     x_size = x_.size();
  MatrixXd I      = MatrixXd::Identity(x_size, x_size);
  P_              = (I - K * H_) * P_;
}
