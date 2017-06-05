#include "kalman_filter.h"

#include <iostream>
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;


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
	/**
		 TODO:
		 * predict the state
		 */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::DoUpdate(const VectorXd &y) 
{
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;	
	
	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}


void KalmanFilter::Update(const VectorXd &z) {
	/**
		 TODO:
		* update the state by using Kalman Filter equations
		*/
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	DoUpdate(y);
}

/**
   Normalizing Angles
	 
   In C++, atan2() returns values between -pi and pi. When
   calculating phi in y = z - h(x) for radar measurements, the
   resulting angle phi in the y vector should be adjusted so that it
   is between -pi and pi. The Kalman filter is expecting small angle
   values between the range -pi and pi. HINT: when working in
   radians, you can add 2π or subtract 2π until the angle is within
   the desired range.
*/

float KalmanFilter::Normalize(float theta) 
{	
	while (theta < -M_PI || theta > M_PI) {
		if (theta < -M_PI) {
	    theta += 2*M_PI;
		}
		else {
	    theta -= 2*M_PI;
		}
    }
    return theta;
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
		 Convert to plar coordinates.
		 */
	double px = x_(0);	
	double py = x_(1);
	double vx = x_(2);
	double vy = x_(3);

	
	double rho = sqrt(px*px + py*py);
	double theta = atan2(py, px);
	double rho_dot;
    
	// avoid devision by 0
	if (fabs(rho) < 0.0001) {
		rho_dot = 0;
	}
	else {
		rho_dot = (px*vx+py*vy) / rho;
	}

	VectorXd z_pred = VectorXd(3);
	z_pred << rho, theta, rho_dot;
    
	VectorXd y = z - z_pred;
	y(1) = Normalize(y(1));
	DoUpdate(y);
}
