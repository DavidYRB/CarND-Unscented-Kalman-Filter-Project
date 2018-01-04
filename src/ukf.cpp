#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.32;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  time_us_ = 0.0;

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_x_;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

  weights_ = VectorXd(2*n_aug_ + 1);

  is_initialized_ = false;
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

  // initialization the filter
  if (!is_initialized_) {
    cout << "UKF: " << endl;

    x_ << 1, 1, 1, 1, 1;
    // initialize state vector with first measurement
    if(meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }
    else{
      float phi = meas_package.raw_measurements_[1];
      float px = meas_package.raw_measurements_[0]*cos(phi);
      float py = meas_package.raw_measurements_[0]*sin(phi);
      x_(0) = px;
      x_(1) = py;
    }

    // initialize covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    // initialize timestamp
    time_us_ = meas_package.timestamp_;

    // initialize lidar H function
    // H_ << 1, 0, 0, 0, 0,
    //       0, 1, 0, 0, 0;

    // initialize weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for(int i=1; i<(2*n_aug_ + 1); i++){
      weights_(i) = 0.5 / (lambda_ + n_aug_);
    }

    // set initialized flag
    is_initialized_ = true;
    cout << "initialization done!" << endl;
    return;

  }

  // prediction
  float d_t = (meas_package.timestamp_ - time_us_ )/1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(d_t);

  // Update prediction according to the measurements
  if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    use_laser_ = true;
    UpdateLidar(meas_package);
    cout << "Lidar update done" << endl;
  }
  else{
    use_radar_ = true;
    // cout << "test before update" <<Xsig_pred_ << endl;
    UpdateRadar(meas_package);
    cout << "Radar update done" << endl;
  }
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
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

  /******** Generate sigma points ***********/
  // 1. get augmented sigma points matrix

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
  // Augment the state vector with noises
  // cout << "Generating augmented sigma points" << endl;
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  // calcualte square root of the covariance Matrix
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  MatrixXd A = P_aug.llt().matrixL();

  // Get augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i=0; i<n_aug_ ; i++){
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)*A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i);
  }
  // cout << "Generating done" << endl;

  // 2. get predicted sigma points matrix
  float px, py, v, yaw, yawd, nu_a, nu_yawdd;
  for(int i=0; i<2*n_aug_+1; i++){
    px = Xsig_aug.col(i)(0);
    py = Xsig_aug.col(i)(1);
    v = Xsig_aug.col(i)(2);
    yaw = Xsig_aug.col(i)(3);
    yawd = Xsig_aug.col(i)(4);
    nu_a = Xsig_aug.col(i)(5);
    nu_yawdd = Xsig_aug.col(i)(6);

    double px_p, py_p;
    if(fabs(yawd) > 0.001){
      px_p = px + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = py + v/yawd * (-cos(yaw + yawd*delta_t) + cos(yaw));
    }
    else{
      px_p= px + v * cos(yaw) * delta_t;
      py_p= py + v * sin(yaw) * delta_t;
    }
    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    px_p += 0.5 * (delta_t*delta_t) * cos(yaw) * nu_a;
    py_p += 0.5 * (delta_t*delta_t) * sin(yaw) * nu_a;
    v_p += delta_t * nu_a;
    yaw_p += 0.5 * (delta_t*delta_t) * nu_yawdd;
    yawd_p += delta_t * nu_yawdd;
    // cout << temp << endl;
    // cout << "problem here" << Xsig_pred_.size() << '\t' << weights_.size() << endl;
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
    // cout << "problem finish" << endl;
  }

  // cout << "prediction of sigmapoints done" << endl;

  // 3. Get mean and covariance matrix from predicted sigma point matrix
  // mean
  VectorXd x_temp = VectorXd(n_x_);
  // cout << "Printing predicted sigmaponts" << endl;
  x_temp.setZero(n_x_);
  for(int i=0; i<(2*n_aug_+1); i++){
    x_temp += weights_(i) * Xsig_pred_.col(i);
    // cout << x_temp << endl;
  }
  x_ = x_temp;
  // cout << x_temp << endl;
  // covariance
  P_.fill(0.0);
  for(int i=0; i<(2*n_aug_+1); i++){
    // residual
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
  cout << Xsig_pred_ << endl;
  cout << "prediction done" << endl;
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
  /***************** Predict Measurement mean and covariance from Sigma points **********/
  cout << "updading lidar" << endl;
  // predict measurement sigma points according to the predicted sigma points
  VectorXd z = meas_package.raw_measurements_;

  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  VectorXd z_pred = VectorXd(n_z);
  Zsig.fill(0.0); z_pred.fill(0.0);

  for(int i=0; i< 2*n_aug_+1; i++){
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }
  // mean
  for(int i=0; i<2*n_aug_+1; i++){
    z_pred += weights_(i) * Zsig.col(i);
  }
  // covariance
  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);
  for(int i=0; i<(2*n_aug_+1); i++){
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;

  S = S + R;
  // cout << "Covariance calculation done" << S << endl;
  /***************** End of Predict Measurement mean and covariance from Sigma points **********/

  /***************** Update the prediction based on the predicted measurement mean and covariance **********/
  // cross correlation matrix
  MatrixXd T = MatrixXd(n_x_, n_z);
  T.fill(0.0);
  for(int i=0; i< 2*n_x_+1; i++){
    // Angle normalization for measurements
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Angle normalization for predicted value
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    T += weights_(i) * x_diff * z_diff.transpose();
  }
  // Kalman gain K
  MatrixXd K = MatrixXd(n_x_, 2);
  K = T * S.inverse();
  // residual
  VectorXd z_diff = z - z_pred;

  // Update state vector
  x_ = x_ + K * z_diff;
  // Update covariance matrix
  P_ = P_ - K * S * K.transpose();

  /***************** End of Update the prediction based on the predicted measurement mean and covariance **********/

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

  /***************** Predict Measurement mean and covariance from Sigma points **********/
  cout << "updading radar" << endl;

  VectorXd z = meas_package.raw_measurements_;

  int n_z = 3;
  // predict measurement sigma points according to the predicted sigma points
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  // cout << Xsig_pred_ << endl;
  Zsig.fill(0.0);
  double px, py, v, yaw;

  for(int i=0; i<(2*n_aug_+1); i++){
    px = Xsig_pred_(0, i);
    py = Xsig_pred_(1, i);
    v = Xsig_pred_(2, i);
    yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    Zsig(0, i) = sqrt(px*px + py*py);
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px*v1 + py*v2)/sqrt(px*px + py*py);
  }
  // get mean of measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.setZero(n_z);
  for(int i=0; i<(2*n_aug_+1); i++){
    z_pred += weights_(i) * Zsig.col(i);
  }
  // get covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for(int i=0; i<(2*n_aug_+1); i++){
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Angle normalization for measurement z
    while(z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    S += weights_(i)* z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;

  S = S + R;
  /***************** end of Predict Measurement mean and covariance from Sigma points **********/

  /***************** Update the prediction based on the predicted measurement mean and covariance **********/
  // Cross correlation matrix
  MatrixXd T = MatrixXd(n_x_, n_z);
  T.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++){

    // Angle normalization for measurements
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while(z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    // Angle normalization for predicted value
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    T += weights_(i) * x_diff * z_diff.transpose();

  }
  // Kalman gain K
  MatrixXd K = T * S.inverse();
  // angle normalization
  VectorXd z_diff = z - z_pred;
  while(z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while(z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
  // Update state vector
  x_ = x_ + K * z_diff;
  // Update covariance matrix
  P_ = P_ - K * S * K.transpose();
}

// void UKF::GetAugmentedSigmaPoints(MatrixXd *Xsig_out){
//   // Augment the state vector with noises
//   cout << "Generating augmented sigma points" << endl;
//   VectorXd x_aug = VectorXd(n_aug_);
//   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
//   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
//   x_aug.head(5) = x_;
//   x_aug(5) = 0;
//   x_aug(6) = 0;
//   // calcualte square root of the covariance Matrix
//   P_aug.topLeftCorner(5,5) = P_;
//   P_aug(5,5) = std_a_ * std_a_;
//   P_aug(6,6) = std_yawdd_ * std_yawdd_;
//   P_aug(5,6) = P_aug(6,5) = 0;
//   MatrixXd A = P_aug.llt().matrixL();
//
//   // Get augmented sigma points
//   Xsig_aug.col(0) = x_aug;
//   for(int i=0; i<n_aug_ ; i++){
//     Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)*A.col(i);
//     Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i);
//   }
//
//   *Xsig_out = Xsig_aug;
//   cout << "Generating done" << endl;
// }
//
// void UKF::PredictSigmaPoints(MatrixXd Xsig_in, double delta_t){
//   float px, py, v, phi, phi_dot;
//   VectorXd temp = VectorXd(n_x_);
//   Xsig_pred_.fill(0.0);
//   cout << Xsig_in << endl;
//   cout << Xsig_pred_ << endl;
//   for(int i=0; i<2*n_aug_+1; i++){
//     temp.setZero(n_x_);
//     px = Xsig_in.col(i)(0);
//     py = Xsig_in.col(i)(1);
//     v = Xsig_in.col(i)(2);
//     phi = Xsig_in.col(i)(3);
//     phi_dot = Xsig_in.col(i)(4);
//
//     if(phi_dot == 0){
//       temp(0) = v * cos(phi) * delta_t;
//       temp(1) = v * sin(phi) * delta_t;
//     }
//     else{
//       temp(0) = v/phi_dot * (sin(phi + phi_dot*delta_t) - sin(phi));
//       temp(1) = v/phi_dot * (-cos(phi + phi_dot*delta_t) + cos(phi));
//       temp(3) = phi_dot * delta_t;
//     }
//
//     temp(0) += 0.5 * (delta_t*delta_t) * cos(phi) * Xsig_in.col(i)(5);
//     temp(1) += 0.5 * (delta_t*delta_t) * sin(phi) * Xsig_in.col(i)(5);
//     temp(2) += delta_t * Xsig_in.col(i)(5);
//     temp(3) += 0.5 * (delta_t*delta_t) * Xsig_in.col(i)(6);
//     temp(4) += delta_t * Xsig_in.col(i)(6);
//     cout << temp << endl;
//     Xsig_pred_.col(i) = temp + Xsig_in.col(i).head(n_x_);
//   }
//   cout << Xsig_pred_ << endl;
// }
//
// void UKF::PredictedMeanAndCov(){
//   // mean
//   VectorXd x_temp = VectorXd(n_x_);
//   cout << Xsig_pred_ << endl;
//   cout << "Printing predicted sigmaponts" << endl;
//   x_temp.setZero(n_x_);
//   for(int i=0; i<(2*n_aug_+1); i++){
//     x_temp += weights_(i) * Xsig_pred_.col(i);
//     cout << x_temp << endl;
//   }
//   x_ = x_temp;
//   cout << x_temp << endl;
//   // covariance
//   P_.setZero(n_x_, n_x_);
//   for(int i=0; i<(2*n_aug_+1); i++){
//     // residual
//     VectorXd x_diff = Xsig_pred_.col(i) - x_;
//
//     // angle normalization
//     while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
//     while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
//
//     P_ += weights_(i) * x_diff * x_diff.transpose();
//   }
//
// }
