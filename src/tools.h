#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
 public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, 
                                const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  Eigen::MatrixXd CalculateRadarJacobian(const Eigen::VectorXd& x_state);

  static Eigen::VectorXd ToPolar(const Eigen::VectorXd& carthesian);
  static Eigen::VectorXd ToCarthesianXY(const Eigen::VectorXd& polar);

};

#endif  // TOOLS_H_
