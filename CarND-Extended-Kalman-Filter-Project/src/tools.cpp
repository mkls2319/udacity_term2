#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// This code is more or less the same as explained on the conferences.
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
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

// This code is more or less the same as explained on the conferences.
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);

	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);


	double c1 = px*px+py*py;
	double c2 = sqrt(c1);
	double c3 = (c1*c2);


	if(fabs(c1) < 0.0001){
		return Hj;
	}

	Hj << (px/c2), (py/c2), 0, 0,
		  - (py/c1), (px/c1), 0, 0,
		    py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
