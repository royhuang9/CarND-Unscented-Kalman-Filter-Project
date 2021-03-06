#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"

enum class SensorType{
    LASER,
    RADAR
};


class MeasurementPackage {
public:
  long long timestamp_;

  SensorType sensor_type_;

  Eigen::VectorXd raw_measurements_;

};

#endif /* MEASUREMENT_PACKAGE_H_ */
