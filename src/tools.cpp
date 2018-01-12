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
  
  // Exception handling
  if (estimations.size() == 0){
      std::cout << "No measurements found" << std::endl;
      return rmse;

  }
  if (estimations.size() != ground_truth.size()){
      std::cout << "Dimension mismatch" << std::endl;
      return rmse;
  }
  
  // Final computation
  rmse << 0,0,0,0;
  
  VectorXd diff;
  for (int i = 0; i < estimations.size(); i++){
      
      diff = estimations[i] - ground_truth[i];
      rmse += diff.array()*diff.array();
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  
  // Follows the equations given in the class lectures
  
  MatrixXd Hj(3, 4);
  
  float eps = 1E-3;
  
  float rho_sq = x_state[0] * x_state[0] + x_state[1] * x_state[1]
  float rho = sqrt(rho_sq)
  float rho_32 = rho*rho_sq;
  
  // Exception handling
  if (fabs(rho_sq) < eps){
      std::cout << "ERROR: Division by zero at CalculateJacobian" << std::endl;
      return Hj;
  }
  
  // Otherwise, just repeat the equations given in the lectures
  Hj <<  x_state[0]/rho   , x_state[1]/rho   , 0, 0,
        -x_state[1]/rho_sq, x_state[0]/rho_sq, 0, 0,
        x_state[1]*(x_state[2]*x_state[1] - x_state[3]*x_state[0])/rho_32, x_state[0]*(x_state[3]*x_state[0] - x_state[2]*x_state[1])/rho_32, x_state[0]/rho   , x_state[1]/rho;

  return Hj;
}