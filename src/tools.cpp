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

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
}
