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
  cout << "RMES" << endl;
  VectorXd rmse(4);
  rmse.fill(0.0);
  // check if the inputs are valid
  if(estimations.size() != ground_truth.size() || estimations.size() == 0)
    return rmse;

  // calculate square residuals
  for(unsigned int i=0; i<estimations.size(); i++){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // final rmse
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}
