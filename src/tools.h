#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
 public:
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                       const std::vector<Eigen::VectorXd> &ground_truth);

  static Eigen::MatrixXd CalculateRadarJacobian(const Eigen::VectorXd& x_state);

  static Eigen::MatrixXd BuildNoiseMatrix(double noise_ax, double noise_ay, double dt);


  static Eigen::VectorXd ToPolar(const Eigen::VectorXd& carthesian);
  static Eigen::VectorXd ToCarthesianXY(const Eigen::VectorXd& polar);

};

#endif  // TOOLS_H_
