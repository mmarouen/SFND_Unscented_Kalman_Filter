#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  is_initialized_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 3;


  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_= 1;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  n_x_=5;

  // Augmented state dimension
  n_aug_=7;

  // Sigma point spreading parameter
  lambda_=3-n_aug_;

  weights_ = VectorXd(2*n_aug_+1);
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_){
    //use first reading as initialization
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2]; 
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
  	  double vy = rho_dot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);
      if(cos(phi)<0) v*=-1;
      x_ << x, y, 0, 0, 0;
    } 
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)  {
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      x_ << x, y, 0, 0, 0;
    }

    P_ <<
        1.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 1.0;
    // initalize time_us
    time_us_ = meas_package.timestamp_ ;

    is_initialized_ = true;
    return;
  }

  dt_ = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  PredictMotion(meas_package);
  
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    ProcessRadarMeasurement(meas_package);
  }
  if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    ProcessLidarMeasurement(meas_package);
  }

}

///// Motion model prediction
void UKF::PredictMotion(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */


  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_ = VectorXd(2*n_aug_+1);

  GenerateSigmaPoints();
  SigmaPointPrediction(dt_);
  PredictMeanAndCovariance();
}

void UKF::GenerateSigmaPoints()
{

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  // create augmented sigma points
  Xsig_aug_.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug_.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

}

void UKF::SigmaPointPrediction(double delta_t)
{

  // predict sigma points
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);
    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + (v/yawd) * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + (v/yawd) * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

}

void UKF::PredictMeanAndCovariance()
{

  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; ++i) {  // 2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // predicted state mean
  x.fill(0.0);
  //x=x_;
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    x = x + weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P.fill(0.0);
  //P=P_;
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }
  x_ = x;
  P_ = P;
}

///// lidar sensor part
void UKF::ProcessLidarMeasurement(MeasurementPackage meas_package)
{
  n_z_= 2;
  S_ = MatrixXd(n_z_,n_z_);
  z_pred_ = VectorXd(n_z_);
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  H_ = MatrixXd(n_z_,n_x_);
  R_ = MatrixXd(n_z_,n_z_);

  H_ << 
    1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0;

  R_ << 
    std_laspx_*std_laspx_, 0.0,
    0.0, std_laspy_*std_laspy_;

  PredictLidarMeasurement();
  UpdateLidar(meas_package);
  //UpdateRadar(meas_package);

}

void UKF::UpdateLidar(MeasurementPackage meas_package) 
{

  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  VectorXd z = VectorXd(n_z_);
  z=meas_package.raw_measurements_;
  VectorXd x1 = VectorXd(n_x_);
  x1 = x_;
  
  VectorXd y = z - z_pred_;
  MatrixXd Ht = H_.transpose();
  MatrixXd Si = S_.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x1 = x1 + (K * y);
  x_=x1;
  //x_(0)=x1(0);
  //x_(1)=x1(1);

  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H_)*P_;   
}

void UKF::PredictLidarMeasurement() 
{
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  z_pred_ = H_ * x_;
  Zsig_ = H_ * Xsig_pred_;
  S_ = H_ * P_ * H_.transpose() + R_;
}

///// radar sensor part
void UKF::ProcessRadarMeasurement(MeasurementPackage meas_package)
{
  n_z_= 3;
  S_ = MatrixXd(n_z_,n_z_);
  z_pred_ = VectorXd(n_z_);
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  R_ = MatrixXd(n_z_,n_z_);

  R_ <<  
    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0,std_radrd_*std_radrd_;

  PredictRadarMeasurement();
  UpdateRadar(meas_package);

}

void UKF::PredictRadarMeasurement() 
{
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
  }

  // mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  S = S + R_;

  // write result
  Zsig_=Zsig;
  z_pred_ = z_pred;
  S_ = S;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) 
{
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */


  //incoming radar measurement
  VectorXd z = VectorXd(n_z_);
  z=meas_package.raw_measurements_;

  VectorXd x1 = VectorXd(n_x_);
  x1 = x_;

  // cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  // residual
  VectorXd z_diff = z - z_pred_;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x1 = x1 + K * z_diff;
  x_= x1;
  //x_(2)=x1(2);
  //x_(3)=x1(3);
  //x_(4)=x1(4);

  P_ = P_ - K*S_*K.transpose();
}

