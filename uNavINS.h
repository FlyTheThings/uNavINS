/*
uNavINS.h

Original Author:
Adhika Lie
2012-10-08
University of Minnesota
Aerospace Engineering and Mechanics
Copyright 2011 Regents of the University of Minnesota. All rights reserved.

Updated to be a class, use Eigen, and compile as an Arduino library.
Added methods to get gyro and accel bias. Added initialization to
estimated angles rather than assuming IMU is level. Added method to get psi,
rather than just heading, and ground track.
Brian R Taylor
brian.taylor@bolderflight.com
2017-12-20
Bolder Flight Systems
Copyright 2017 Bolder Flight Systems

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

// XXX - accel and gyro bias not being updated.
// XXX - add set methods for sensor characteristics.
// XXX - is psi heading or track? I think heading.
// XXX - add outputs for filter covariance to measure convergence.
// XXX - incorporate magnetometers.


#ifndef UNAVINS_H
#define UNAVINS_H
#if defined(ARDUINO)
  #include "Arduino.h"
  #include "Eigen.h"
  #include <Eigen/Dense>
#else
  #include <sys/time.h>
  #include <stdint.h>
  #include <math.h>
  #include <Eigen/Core>
  #include <Eigen/Dense>

  uint64_t micros() {
      struct timeval tv;
      gettimeofday(&tv,NULL);
      return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
  }

  class elapsedMicros
  {
  private:
  	unsigned long us;
  public:
  	elapsedMicros(void) { us = micros(); }
  	elapsedMicros(unsigned long val) { us = micros() - val; }
  	elapsedMicros(const elapsedMicros &orig) { us = orig.us; }
  	operator unsigned long () const { return micros() - us; }
  	elapsedMicros & operator = (const elapsedMicros &rhs) { us = rhs.us; return *this; }
  	elapsedMicros & operator = (unsigned long val) { us = micros() - val; return *this; }
  	elapsedMicros & operator -= (unsigned long val)      { us += val ; return *this; }
  	elapsedMicros & operator += (unsigned long val)      { us -= val ; return *this; }
  	elapsedMicros operator - (int val) const           { elapsedMicros r(*this); r.us += val; return r; }
  	elapsedMicros operator - (unsigned int val) const  { elapsedMicros r(*this); r.us += val; return r; }
  	elapsedMicros operator - (long val) const          { elapsedMicros r(*this); r.us += val; return r; }
  	elapsedMicros operator - (unsigned long val) const { elapsedMicros r(*this); r.us += val; return r; }
  	elapsedMicros operator + (int val) const           { elapsedMicros r(*this); r.us -= val; return r; }
  	elapsedMicros operator + (unsigned int val) const  { elapsedMicros r(*this); r.us -= val; return r; }
  	elapsedMicros operator + (long val) const          { elapsedMicros r(*this); r.us -= val; return r; }
  	elapsedMicros operator + (unsigned long val) const { elapsedMicros r(*this); r.us -= val; return r; }
  };
#endif

class uNavINS {
  public:
    void update(unsigned long TOW,double vn,double ve,double vd,double lat,double lon,double alt,float p,float q,float r,float ax,float ay,float az,float hx,float hy, float hz);
    float getPitch_rad();
    float getRoll_rad();
    float getYaw_rad();
    float getHeading_rad();
    double getLatitude_rad();
    double getLongitude_rad();
    double getAltitude_m();
    double getVelNorth_ms();
    double getVelEast_ms();
    double getVelDown_ms();
    float getGroundTrack_rad();
    float getGyroBiasX_rads();
    float getGyroBiasY_rads();
    float getGyroBiasZ_rads();
    float getAccelBiasX_mss();
    float getAccelBiasY_mss();
    float getAccelBiasZ_mss();
  private:
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // error characteristics of navigation parameters
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Std dev of Accelerometer Wide Band Noise (m/s^2)
    const float SIG_W_A = 1.0f;       // 1 m/s^2
    // Std dev of gyro output noise (rad/s)
    const float SIG_W_G = 0.00524f;   // 0.3 deg/s
    // Std dev of Accelerometer Markov Bias
    const float SIG_A_D = 0.1f;       // 5e-2*g
    // Correlation time or time constant
    const float TAU_A = 100.0f;
    // Std dev of correlated gyro bias
    const float SIG_G_D = 0.00873f;   // 0.1 deg/s
    // Correlation time or time constant
    const float TAU_G = 50.0f;
    // GPS measurement noise std dev (m)
    const float SIG_GPS_P_NE = 3.0f;
    const float SIG_GPS_P_D = 5.0f;
    // GPS measurement noise std dev (m/s)
    const float SIG_GPS_V = 0.5f;
    // Initial set of covariance
    const float P_P_INIT = 10.0f;
    const float P_V_INIT = 1.0f;
    const float P_A_INIT = 0.34906f;     // 20 deg
    const float P_HDG_INIT = 3.14159f;   // 180 deg
    const float P_AB_INIT = 0.9810f;     // 0.5*g
    const float P_GB_INIT = 0.01745f;    // 5 deg/s
    // acceleration due to gravity
    const float G = 9.807f;
    // major eccentricity squared
    const double ECC2 = 0.0066943799901;
    // earth semi-major axis radius (m)
    const double EARTH_RADIUS = 6378137.0;
    // initialized
    bool initialized = false;
    // timing
    elapsedMicros _t;
    float _dt;
    unsigned long previousTOW;
    // estimated attitude
    float phi, theta, psi, heading;
    // initial heading angle
    float psi_initial;
    // estimated NED velocity
    double vn_ins, ve_ins, vd_ins;
    // estimated location
    double lat_ins, lon_ins, alt_ins;
    // magnetic heading corrected for roll and pitch angle
    float Bxc, Byc;
    // accelerometer bias
    float abx, aby, abz;
    // gyro bias
    float gbx, gby, gbz;
    // earth radius at location
    double Re, Rn, denom;
    // State matrix
    Eigen::Matrix<float,15,15> Fs = Eigen::Matrix<float,15,15>::Identity();
    // State transition matrix
    Eigen::Matrix<float,15,15> PHI = Eigen::Matrix<float,15,15>::Zero();
    // Covariance matrix
    Eigen::Matrix<float,15,15> P = Eigen::Matrix<float,15,15>::Zero();
    // For process noise transformation
    Eigen::Matrix<float,15,12> Gs = Eigen::Matrix<float,15,12>::Zero();
    Eigen::Matrix<float,12,12> Rw = Eigen::Matrix<float,12,12>::Zero();
    // Process noise matrix
    Eigen::Matrix<float,15,15> Q = Eigen::Matrix<float,15,15>::Zero();
    // Gravity model
    Eigen::Matrix<float,3,1> grav = Eigen::Matrix<float,3,1>::Zero();
    // Rotation rate
    Eigen::Matrix<float,3,1> om_ib = Eigen::Matrix<float,3,1>::Zero();
    // Specific force
    Eigen::Matrix<float,3,1> f_b = Eigen::Matrix<float,3,1>::Zero();
    // DCM
    Eigen::Matrix<float,3,3> C_N2B = Eigen::Matrix<float,3,3>::Zero();
    // DCM transpose
    Eigen::Matrix<float,3,3> C_B2N = Eigen::Matrix<float,3,3>::Zero();
    // Temporary to get dxdt
    Eigen::Matrix<float,3,1> dx = Eigen::Matrix<float,3,1>::Zero();
    Eigen::Matrix<double,3,1> dxd = Eigen::Matrix<double,3,1>::Zero();
    // NED velocity INS
    Eigen::Matrix<double,3,1> V_ins = Eigen::Matrix<double,3,1>::Zero();
    // LLA INS
    Eigen::Matrix<double,3,1> lla_ins = Eigen::Matrix<double,3,1>::Zero();
    // NED velocity GPS
    Eigen::Matrix<double,3,1> V_gps = Eigen::Matrix<double,3,1>::Zero();
    // LLA GPS
    Eigen::Matrix<double,3,1> lla_gps = Eigen::Matrix<double,3,1>::Zero();
    // Position ECEF INS
    Eigen::Matrix<double,3,1> pos_ecef_ins = Eigen::Matrix<double,3,1>::Zero();
    // Position NED INS
    Eigen::Matrix<double,3,1> pos_ned_ins = Eigen::Matrix<double,3,1>::Zero();
    // Position ECEF GPS
    Eigen::Matrix<double,3,1> pos_ecef_gps = Eigen::Matrix<double,3,1>::Zero();
    // Position NED GPS
    Eigen::Matrix<double,3,1> pos_ned_gps = Eigen::Matrix<double,3,1>::Zero();
    // Quat
    Eigen::Matrix<float,4,1> quat = Eigen::Matrix<float,4,1>::Zero();
    // dquat
    Eigen::Matrix<float,4,1> dq = Eigen::Matrix<float,4,1>::Zero();
    // difference between GPS and INS
    Eigen::Matrix<float,6,1> y = Eigen::Matrix<float,6,1>::Zero();
    // GPS measurement noise
    Eigen::Matrix<float,6,6> R = Eigen::Matrix<float,6,6>::Zero();
    Eigen::Matrix<float,15,1> x = Eigen::Matrix<float,15,1>::Zero();
    // Kalman Gain
    Eigen::Matrix<float,15,6> K = Eigen::Matrix<float,15,6>::Zero();
    Eigen::Matrix<float,6,15> H = Eigen::Matrix<float,6,15>::Zero();
    // skew symmetric
    Eigen::Matrix<float,3,3> sk(Eigen::Matrix<float,3,1> w);
    // lla rate
    Eigen::Matrix<double,3,1> llarate(Eigen::Matrix<double,3,1> V,Eigen::Matrix<double,3,1> lla);
    // lla to ecef
    Eigen::Matrix<double,3,1> lla2ecef(Eigen::Matrix<double,3,1> lla);
    // ecef to ned
    Eigen::Matrix<double,3,1> ecef2ned(Eigen::Matrix<double,3,1> ecef,Eigen::Matrix<double,3,1> pos_ref);
    // quaternion to dcm
    Eigen::Matrix<float,3,3> quat2dcm(Eigen::Matrix<float,4,1> q);
    // quaternion multiplication
    Eigen::Matrix<float,4,1> qmult(Eigen::Matrix<float,4,1> p, Eigen::Matrix<float,4,1> q);
    // maps angle to +/- 180
    float constrainAngle180(float dta);
    // maps angle to 0-360
    float constrainAngle360(float dta);
};

#endif
