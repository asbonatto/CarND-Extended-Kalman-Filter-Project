#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /* NOTE : 
        Bu = 0, since we are not modeling directly the acceleration
  */
  x_ = F_*x_;
  // A posteriori error covariance 
  MatrixXd F_transpose = F_.transpose(); 
  // NOTE : could I avoid this declaration or would it modify the original F? RTFM
  P_ = F_ * P_ * F_transpose + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Receives the measurement and updates the state
  VectorXd y = z - H_ * x_;
  AdjustPrediction(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  // First step : to create the h'(x) function
  double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double phi = atan2(x_(1), x_(0));
  double rho_d;
  double eps = 1E-3;
  
  // Avoiding division by 0
   if (fabs(rho) < eps){
       rho_d = 0;
       phi = 0;
   }else{
        rho_d= (x_(0)*x_(2) + x_(1)*x_(3)) / rho;
    }    
  
  VectorXd h(3);
  h << rho, phi, rho_d;
  
  VectorXd y = z - h;
    // Normalization to [-pi; pi]
  float cosy1 = cos(y(1));
  float siny1 = sin(y(1));
  
  if (fabs(cosy1) < eps) {
      if(siny1 > 0){
           y(1) = M_PI/2;
      }else{
          y(1) = -M_PI/2;
      }
  }else{
      y(1) = atan2(siny1, cosy1);
  }
  
  AdjustPrediction(y);
  
}

void KalmanFilter::AdjustPrediction(const VectorXd &y) {
  /* Generic method to update the state given the error
     measurements to avoid duplication between 
   Update and UpdateEKF
  */ 
  
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  
  MatrixXd H_t = H_.transpose();
  MatrixXd S = H_ * P_ * H_t + R_;
  MatrixXd S_i = S.inverse();
  MatrixXd K =  P_ * H_t * S_i;
  
  x_ = x_ + (K * y);
  
  P_ = (I - K * H_) * P_;
}
