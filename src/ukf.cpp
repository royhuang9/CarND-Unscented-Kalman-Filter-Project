#include "ukf.h"
#include "tools.h"
#include <Eigen/Dense>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  is_initialized_ = false;
  
  /* dimension of x */
  n_x_ = 5;
  
  /* augmented dimension of x */
  n_aug_ = 7;
  
  
  lambda_ = 3 - n_aug_;
  
  /* initialize weights */
  weights_ = VectorXd::Zero(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i=1; i< 2*n_aug_ + 1; ++i) {
      weights_(i) = 0.5/(lambda_ + n_aug_);
  }
  
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2*n_aug_+1);
  
  H_ = MatrixXd::Zero(2, n_x_);
  H_(0,0) = 1.0;
  H_(1,1) = 1.0;
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  /* if not initalized, use the first measurement to initialize x */
  if (!is_initialized_) {
      x_.fill(0.0);
      
      if (meas_package.sensor_type_ == SensorType::RADAR) {
          VectorXd raMeas = meas_package.raw_measurements_;
          
          if ((raMeas(0) < 0.001) && (raMeas(1) < 0.001) && (raMeas(2) < 0.001))
                return;

          x_(0) = raMeas(0) * cos(raMeas(1));
          x_(1) = raMeas(0) * sin(raMeas(1));
          
      } else if (meas_package.sensor_type_ == SensorType::LASER) {
          VectorXd laMeas = meas_package.raw_measurements_;
          if ((laMeas(0) < 0.001) && (laMeas(1)  < 0.001))
                return;
          x_(0) = laMeas(0);
          x_(1) = laMeas(1);
      }
      
      previous_timestamp_ = meas_package.timestamp_;
      
      is_initialized_ = true;
  }
  
  /* calculate the time difference */
  double dt = (meas_package.timestamp_ - previous_timestamp_)/1.0e6;
  
  /*If delta time is too small or zero, skip prediction */
  if (dt > 1e-6) {
    Prediction(dt);
  }
  
  /*update in term of sensor type */
  if (meas_package.sensor_type_ == SensorType::LASER) {
      UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == SensorType::RADAR) {
      UpdateRadar(meas_package);
  }
}
/* 
 * Generate sigma points
 */
void UKF::GenerateSigmaPoints(MatrixXd &Xsig_aug) {
      /* Create x_aug */
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  
  /* Create P_aug */
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;
 
  /* Calculate square root of P_aug */
  MatrixXd A = P_aug.llt().matrixL();
  
  MatrixXd m = sqrt(lambda_ + n_aug_) * A;
  

  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = m.colwise() + x_aug;
  Xsig_aug.block(0, n_aug_+1, n_aug_, n_aug_) = (-1.0) * (m.colwise() - x_aug);
}

void UKF::OnePointPrediction(VectorXd &X_in, VectorXd &X_out, double dt)
{
    double px = X_in(0);
    double py = X_in(1);
    double vel = X_in(2);
    double phi = X_in(3);
    double phi_dot = X_in(4);
    double nu_a = X_in(5);
    double nu_phi_dot_dot = X_in(6);
    
    if (phi_dot < 1e-6) {
        X_out(0) = px + vel*cos(phi)*dt + 0.5*dt*dt*cos(phi)*nu_a;
        X_out(1) = py + vel*sin(phi)*dt + 0.5*dt*dt*sin(phi)*nu_a;
    } else {
        X_out(0) = px + vel/phi_dot*(sin(phi+phi_dot*dt) - sin(phi)) + 0.5*dt*dt*cos(phi)*nu_a;
        X_out(1) = py + vel/phi_dot*(cos(phi) - cos(phi+phi_dot*dt)) + 0.5*dt*dt*sin(phi)*nu_a;
    }
    X_out(2) = vel + dt*nu_a;
    X_out(3) = phi + phi_dot * dt + 0.5*dt*dt*nu_phi_dot_dot;
    X_out(4) = phi_dot + dt*nu_phi_dot_dot;
}


void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, MatrixXd &Xsig_pred, double delta_t) 
{

  for(int i=0; i < Xsig_aug.cols(); i++) {
      VectorXd X_onevec(7);
      X_onevec = Xsig_aug.col(i);
      VectorXd X_onepred(5);
      OnePointPrediction(X_onevec, X_onepred, delta_t);
      Xsig_pred.col(i) = X_onepred;
  }

}

/* calaculate mean and covariance*/
void UKF::PredictMeanAndCovariance(MatrixXd &Xsig_pred, VectorXd &x_out, MatrixXd &P_out)
{
    for (int i=0; i < 2*n_aug_+1; i++) {
        x_out = weights_[i] * Xsig_pred.col(i);
    }
    
    for (int i=0; i < 2*n_aug_+1; i++) {
        VectorXd temp = Xsig_pred.col(i) - x_;
        P_ += weights_[i]*(temp*temp.transpose());
    }
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  /* sigma points matrix */
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, 2*n_aug_ + 1);
  GenerateSigmaPoints(Xsig_aug);
  
  /* Predict sigma points */
  SigmaPointPrediction(Xsig_aug, Xsig_pred_, delta_t);
  
  /* mean and covariance of predicted sigma points */
  PredictMeanAndCovariance(Xsig_pred_, x_, P_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_prime = H_ * x_;
  VectorXd y = z - z_prime;
  
  MatrixXd S = H_ * P_ * H_.transpose();
  S(0,0) += std_laspx_ * std_laspx_;
  S(1,1) += std_laspy_ * std_laspy_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  
  x_ = x_ + K*y;
  
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
  /* calculate the lidar NIS */
  NIS_laser_ = (z -  z_prime).transpose()* S.inverse() * (z - z_prime);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;
  
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2*n_aug_+1);
  VectorXd z_pred = VectorXd::Zero(n_z);
  
  MatrixXd S=MatrixXd::Zero(n_z, n_z);
    for(int i = 0; i < 2*n_aug_ + 1; i++) {
      double px = Xsig_pred_.col(i)(0);
      double py = Xsig_pred_.col(i)(1);
      double vel = Xsig_pred_.col(i)(2);
      double phi = Xsig_pred_.col(i)(3);
      double phi_dot = Xsig_pred_.col(i)(4);
      
      double rho = sqrt(px*px + py*py);
      double psi = atan2(py, px);
      double psi_dot = (px*cos(phi) + py*sin(phi))*vel/rho;
      
      Zsig.col(i)(0) = rho;
      Zsig.col(i)(1) = psi;
      Zsig.col(i)(2) = psi_dot;
  }
  
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1; i++) {
      z_pred += weights_(i) * Zsig.col(i);
  }
  
  for (int i=0; i < 2*n_aug_ + 1; i++) {
      VectorXd temp = Zsig.col(i) - z_pred;
      
      S += weights_(i) * temp *  temp.transpose();
  }
  
  S(0,0) += std_radr_ * std_radr_;
  S(1,1) += std_radphi_ * std_radphi_;
  S(2,2) += std_radrd_ * std_radrd_;
  
  
  VectorXd z = meas_package.raw_measurements_;
  
  /* Cross correlation */
  MatrixXd T = MatrixXd::Zero(n_x_, n_z);
  for (int i=0; i<2*n_aug_+1; i++) {
      T += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }
  
  /* Kalman Gain */
  MatrixXd Kg = T * S.inverse();
  
  /* state update */
  x_ = x_ + Kg * (z - z_pred);
  
  /* Covariance matrix update */
  P_ = P_ - Kg * S * Kg.transpose();
  
  /* calculate NIS of radar */
  NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
}
