#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  if (estimations.empty() or estimations.size() != ground_truth.size())
  {
    cout << "Error!" << endl;
  }

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  VectorXd residuals(4);
  residuals << 0, 0, 0, 0;
  for (int i = 0; i < estimations.size(); ++i)
  {
    VectorXd diff;
    diff                  = estimations[i] - ground_truth[i];
    VectorXd squared_diff = diff.array() * diff.array();
    residuals             = residuals + squared_diff;
  }

  VectorXd mean = residuals / estimations.size();
  rmse = mean.array().sqrt();
  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
  MatrixXd Hj(3, 4);

  // recover state parameters
  const double px = x_state(0);
  const double py = x_state(1);
  const double vx = x_state(2);
  const double vy = x_state(3);

  // check division by zero
  if ((px == 0.) && (py == 0.))
  {
    cout << "Error! Both Px and Py are 0" << endl;
    return Hj;
  }

  const double squared_sum   = px * px + py * py;
  const double kappa         = 1. / sqrt(squared_sum);
  const double recipr_sq_sum = 1. / squared_sum;
  const double beta          = pow(squared_sum, 3. / 2.);
  const double kappa_px      = kappa * px;
  const double kappa_py      = kappa * py;

  Hj(0, 0) = kappa_px;
  Hj(0, 1) = kappa_py;
  Hj(0, 2) = 0.;
  Hj(0, 3) = 0.;

  Hj(1, 0) = -1. * py * recipr_sq_sum;
  Hj(1, 1) = px * recipr_sq_sum;
  Hj(1, 2) = 0.;
  Hj(1, 3) = 0.;

  Hj(2, 0) = (vx * py - vy * px) * py / beta;
  Hj(2, 1) = (vy * px - vx * py) * px / beta;
  Hj(2, 2) = kappa_px;
  Hj(2, 3) = kappa_py;

  return Hj;
}
