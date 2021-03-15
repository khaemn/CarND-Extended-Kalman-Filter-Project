#include "tools.h"

#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

using std::cout;
using std::endl;

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  if (estimations.empty() or estimations.size() != ground_truth.size())
  {
    cout << "Error! Empty estimations vector or unmatched ground truth vector size." << endl;
  }

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  VectorXd residuals(4);
  residuals << 0, 0, 0, 0;
  for (size_t i = 0; i < estimations.size(); ++i)
  {
    VectorXd diff;
    diff                  = estimations[i] - ground_truth[i];
    VectorXd squared_diff = diff.array() * diff.array();
    residuals             = residuals + squared_diff;
  }

  VectorXd mean = residuals / estimations.size();
  rmse          = mean.array().sqrt();
  // return the result
  return rmse;
}

MatrixXd Tools::CalculateRadarJacobian(const VectorXd &x_state)
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

  const double squared_sum = px * px + py * py;
  const double inv_rho     = 1. / sqrt(squared_sum);
  const double beta        = pow(squared_sum, 3. / 2.);
  const double inv_rho_px  = inv_rho * px;
  const double inv_rho_py  = inv_rho * py;

  Hj(0, 0) = inv_rho_px;
  Hj(0, 1) = inv_rho_py;
  Hj(0, 2) = 0.;
  Hj(0, 3) = 0.;

  Hj(1, 0) = -1. * py / squared_sum;
  Hj(1, 1) = px / squared_sum;
  Hj(1, 2) = 0.;
  Hj(1, 3) = 0.;

  Hj(2, 0) = (vx * py - vy * px) * py / beta;
  Hj(2, 1) = (vy * px - vx * py) * px / beta;
  Hj(2, 2) = inv_rho_px;
  Hj(2, 3) = inv_rho_py;

  return Hj;
}

Eigen::MatrixXd Tools::BuildNoiseMatrix(double noise_ax, double noise_ay, double dt)
{
  // Builds the Q matrix using X and Y acceleration noise
  // and a given delta-t

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

Eigen::VectorXd Tools::ToPolar(const Eigen::VectorXd &carthesian)
{
  VectorXd polar(3);
  polar << 0, 0, 0;
  if (carthesian.size() != 4)
  {
    cout << "Error: Carthesian position and velocities vector must be of size 4 while is "
         << carthesian.size() << endl;
    return polar;
  }

  // recover state parameters
  const double px = carthesian(0);
  const double py = carthesian(1);
  const double vx = carthesian(2);
  const double vy = carthesian(3);

  const double rho     = sqrt(px * px + py * py);
  double       phi     = atan2(py, px);
  const double rho_dot = (px * vx + py * vy) / (rho + 0.001);

  polar << rho, phi, rho_dot;
  return polar;
}

Eigen::VectorXd Tools::ToCarthesian(const Eigen::VectorXd &polar)
{
  // As the rho_dot does not provide enought information about
  // true values of vx and vy, it is only possible to convert
  // the x and y coordinates, but not velocities.

  const double rho = polar(0);
  const double phi = polar(1);
  const double x   = rho * cos(phi);
  const double y   = rho * sin(phi);

  VectorXd carthesian(4);
  carthesian << x, y, 0, 0;

  return carthesian;
}
