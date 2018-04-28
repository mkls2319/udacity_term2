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
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  unsigned int n = estimations.size();
  if((estimations.size() == 0) || (ground_truth.size() == 0) || (n != ground_truth.size())){
  return rmse;
  }

  for(unsigned int i=0; i < estimations.size(); ++i){
  VectorXd res = estimations[i] - ground_truth[i];
  res = res.array()*res.array();
  rmse += res;
  }

  rmse = rmse / n;
  rmse = rmse.array().sqrt();
  return rmse;
}
