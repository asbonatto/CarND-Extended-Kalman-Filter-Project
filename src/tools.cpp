#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd rmse(4);
  
  // Initialize with zeros, otherwise it can blow upper_bound
  rmse << 0,0,0,0;
  
  // Exception handling
  if (estimations.size() == 0){
      std::cout << "No measurements found" << std::endl;
      return rmse;

  }
  if (estimations.size() != ground_truth.size()){
      std::cout << "Dimension mismatch" << std::endl;
      return rmse;
  }
  
  
  
  VectorXd diff;
  for (int i = 0; i < estimations.size(); i++){
      
      diff = estimations[i] - ground_truth[i];
      diff = diff.array()*diff.array();
      rmse += diff;
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  
  // Follows the equations given in the class lectures
  
  float eps = 1E-3;
  
  MatrixXd Hj(3, 4);
 
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  float rho_sq = pow(px, 2) + pow(py,2);
  float rho = sqrt(rho_sq);
  float rho_32 = rho*rho_sq;
  
  // Exception handling
  if (fabs(rho_sq) < eps){
      std::cout << "ERROR: Division by zero at CalculateJacobian" << std::endl;
      return Hj;
  }
  
  // Otherwise, just repeat the equations given in the lectures
  Hj <<  px/rho   , py/rho   , 0, 0,
        -py/rho_sq, px/rho_sq, 0, 0,
         py*(vx*py - vy*px) / rho_32, px*(px*vy - py*vx) / rho_32, px / rho, py / rho;

  return Hj;
}