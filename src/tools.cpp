#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  
    /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here
  if (estimations.size() == 0) {
    return rmse;
  }
  if (ground_truth.size() != estimations.size()) {
    return rmse;
  }

  vector<VectorXd> diff(estimations.size());
  
  //Calculate the difference between estimation and ground truth
  std::transform(estimations.begin(), estimations.end(), ground_truth.begin(), \
                diff.begin(), [](VectorXd d1, VectorXd d2) { return d1 - d2;});
                
  //Calculate the square of every element
  std::transform(diff.begin(), diff.end(), diff.begin(), 
                [](VectorXd d){return (VectorXd)(d.array()*d.array());});
  
  // Calculate the sum of all the element
  rmse = std::accumulate(diff.begin(), diff.end(), rmse, 
                [](VectorXd d1, VectorXd d2) { return d1 + d2;} );

  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}
