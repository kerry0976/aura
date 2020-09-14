/*! \file nav_functions.h
 *	\brief Auxiliary functions for nav filter header file
 *
 *	\details
 *     Module:          navfunc.h
 *     Modified:        Brian Taylor (convert to eigen3)
 *						Gokhan Inalhan (remaining) 
 *                      Demoz Gebre (first three functions)
 *                      Adhika Lie
 *                      Jung Soon Jang
 *     Description:     navfunc.h contains all the variable, 
 *                      constants and function prototypes that are 
 *                      used with the inertial navigation software.
 *	\ingroup nav_fcns
 *
 * \author University of Minnesota
 * \author Aerospace Engineering and Mechanics
 * \copyright Copyright 2011 Regents of the University of Minnesota. All rights reserved.
 *
 * $Id: nav_functions.h 922 2012-10-17 19:14:09Z joh07594 $
 */

#pragma once

#include "structs.h"
#include <stdint.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
using namespace Eigen;

// Constants
const double EarthRadius = 6378137.0;        // earth semi-major axis radius (m)
const double ECC2 = 0.0066943799901;         // major eccentricity squared

// Constants for tightly coupled EKF
const double MU = 3.986005e14; //m^3 / sec^2
const double OMEGA_DOT_EARTH = 7.2921151467e-5; // rad/sec

// Constants that are no longer used
// const double EarthRate = 0.00007292115;      // rotation rate of earth (rad/sec)
// const double Eccentricity = 0.0818191908426; // major eccentricity of earth ellipsoid
// const double Flattening = 0.0033528106650;   // flattening of the ellipsoid
// const double Gravity0 = 9.7803730;           // zeroth coefficient for gravity model
// const double Gravity1 = 0.0052891;           // first coefficient for the gravity model
// const double Gravity2 = 0.0000059;           // second coefficient for the gravity model
// const double GravityNom = 9.81;              // nominal gravity
// const double Schuler2 = 1.533421593170545E-06; // Schuler Frequency (rad/sec) Squared

// This function calculates the rate of change of latitude, longitude,
// and altitude using WGS-84.
Vector3f llarate(Vector3f V, Vector3d lla);

// This function calculates the angular velocity of the NED frame,
// also known as the navigation rate using WGS-84.
Vector3d navrate(Vector3d V, Vector3d lla);

// This function calculates the ECEF Coordinate given the
// Latitude, Longitude and Altitude.
Vector3d lla2ecef(Vector3d lla);

// This function calculates the Latitude, Longitude and Altitude given
// the ECEF Coordinates.
Vector3d ecef2lla( Vector3d ecef_pos );
    
// This function converts a vector in ecef to ned coordinate centered
// at pos_ref.
Vector3f ecef2ned(Vector3d ecef, Vector3d pos_ref);

// Return a quaternion rotation from the earth centered to the
// simulation usual horizontal local frame from given longitude and
// latitude.  The horizontal local frame used in simulations is the
// frame with x-axis pointing north, the y-axis pointing eastwards and
// the z axis pointing downwards.  (Returns the ecef2ned
// transformation as a quaternion.)
Quaterniond lla2quat(double lon_rad, double lat_rad);

// This function gives a skew symmetric matrix from a given vector w
Matrix3f sk(Vector3f w);

// This function gives a skew symmetric matrix from a given vector w; using double
Matrix3d sk2(Vector3d w);

// Quaternion to euler angle: returns phi, the, psi as a vector
Vector3f quat2eul(Quaternionf q);

// Computes a quaternion from the given euler angles
Quaternionf eul2quat(float phi, float the, float psi);

// Quaternion to C_N2B
Matrix3f quat2dcm(Quaternionf q);

// EphemerisData (subframe1,2,3) to Satellite ecef x, y, z in meter, vx, vy, vz in m/s
VectorXd EphemerisData2Satecef(float t,
                               uint32_t TOW, uint8_t L2, uint16_t week_No, uint8_t L2_Flag, uint8_t SV_Acc, uint8_t SV_Hlth,
                               double T_GD, uint16_t IODC, double t_OC, int8_t a_f2, double a_f1, double a_f0,
                               uint8_t IODE, double C_rs, double delta_n, double M_0, double C_uc, double ecc, double C_us,
                               double sqrt_A, double t_OE, double C_ic, double Omega_0, double C_is, double i_0, double C_rc,
                               double omega, double Omega_dot, double IDOT);


// Compute direction cosine matric C from a Euler vector 
// eul = [yaw,pitch,roll]. (i.e., 3-2-1 rotation convention)
Matrix3d eul2dcm(Vector3d eul);

// ned 2 ecef centered at the coordinate given by lla
Vector3d ned2ecef(Vector3d ned, Vector3d lla);