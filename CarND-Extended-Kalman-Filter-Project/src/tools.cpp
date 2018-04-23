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

  if( (estimations.size() != 0) && (estimations.size() == ground_truth.size()))
  {

    for(unsigned int i=0; i < estimations.size(); ++i)
    {
        VectorXd residual = estimations[i] - ground_truth[i];

        residual = residual.array()*residual.array();
        rmse += residual;
    }


    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
  }
  else
  {
    return rmse;
  }
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);

  float x = x_state(0);
  float y = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);


  float eq1 = x*x + y*y;
  float eq2 = sqrt(eq1);
  float eq3 = (eq1*eq2);

  //Check divide by zero
  if (fabs(eq1) < 0.0001){
      return Hj;
  }


  Hj << (x / eq2), (y / eq2), 0, 0,
       -(y / eq1), (x / eq1), 0, 0,
        y*(vx*y - vy*x) / eq3, x*(x*vy - y*vx) / eq3, x / eq2, y / eq2;

  return Hj;
}
