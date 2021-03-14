#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"
#include <map>

class MeasurementPackage {
 public:
  enum SensorType{
    LASER,
    RADAR
  } sensor_type_;

  static const std::map<SensorType, long> SIZES;

  long long timestamp_;

  Eigen::VectorXd raw_measurements_;
};


#endif // MEASUREMENT_PACKAGE_H_
