/*! \file nav_functions.c
 *	\brief Auxiliary functions for nav filter
 *
 *	\details
 *     Module:          Navfuncs.c
 *     Modified:        Brian Taylor (convert to eigen3)
 *						Adhika Lie (revamp all functions)
 * 						Gokhan Inalhan (remaining)
 *                      Demoz Gebre (first three functions)
 *                      Jung Soon Jang
 *
 *     Description:     navfunc.c contains the listing for all the
 *                      real-time inertial navigation software.
 *
 *		Note: all the functions here do not create memory without
 *			  clearing it.
 *	\ingroup nav_fcns
 *
 * \author University of Minnesota
 * \author Aerospace Engineering and Mechanics
 * \copyright Copyright 2011 Regents of the University of Minnesota. All rights reserved.
 *
 * $Id: nav_functions.c 922 2012-10-17 19:14:09Z joh07594 $
 */

/*     Include Pertinent Header Files */

#include <math.h>
#include <stdint.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
using namespace Eigen;

#include "nav_functions.h"
#include "structs.h"

// This function calculates the rate of change of latitude, longitude,
// and altitude using WGS-84.
Vector3f llarate(Vector3f V, Vector3d lla)
{
    double lat = lla(0, 0);
    double h = lla(2, 0);

    double denom = fabs(1.0 - (ECC2 * sin(lat) * sin(lat)));
    double sqrt_denom = sqrt(denom);

    double Rew = EarthRadius / sqrt_denom;
    double Rns = EarthRadius * (1 - ECC2) / (denom * sqrt_denom);

    Vector3f lla_dot;
    lla_dot(0, 0) = V(0, 0) / (Rns + h);
    lla_dot(1, 0) = V(1, 0) / ((Rew + h) * cos(lat));
    lla_dot(2, 0) = -V(2, 0);

    return lla_dot;
}

// This function calculates the angular velocity of the NED frame,
// also known as the navigation rate using WGS-84.
Vector3d navrate(Vector3d V, Vector3d lla)
{
    double lat = lla(0, 0);
    double h = lla(2, 0);

    double denom = fabs(1.0 - (ECC2 * sin(lat) * sin(lat)));
    double sqrt_denom = sqrt(denom);

    double Rew = EarthRadius / sqrt_denom;
    double Rns = EarthRadius * (1 - ECC2) / (denom * sqrt_denom);

    Vector3d nr;
    nr(0, 0) = V(1, 0) / (Rew + h);
    nr(1, 0) = -V(0, 0) / (Rns + h);
    nr(2, 0) = -V(1, 0) * tan(lat) / (Rew + h);

    return nr;
}

// This function calculates the ECEF Coordinate given the
// Latitude, Longitude and Altitude.
Vector3d lla2ecef(Vector3d lla)
{
    double sinlat = sin(lla(0, 0));
    double coslat = cos(lla(0, 0));
    double coslon = cos(lla(1, 0));
    double sinlon = sin(lla(1, 0));
    double alt = lla(2, 0);

    double denom = fabs(1.0 - (ECC2 * sinlat * sinlat));

    double Rew = EarthRadius / sqrt(denom);

    Vector3d ecef;
    ecef(0, 0) = (Rew + alt) * coslat * coslon;
    ecef(1, 0) = (Rew + alt) * coslat * sinlon;
    ecef(2, 0) = (Rew * (1.0 - ECC2) + alt) * sinlat;

    return ecef;
}

// This function calculates the Latitude, Longitude and Altitude given
// the ECEF Coordinates.
Vector3d ecef2lla(Vector3d ecef_pos)
{
    const double Squash = 0.9966471893352525192801545;
    const double ra2 = 1.0 / (EarthRadius * EarthRadius);
    const double e2 = fabs(1 - Squash * Squash);
    const double e4 = e2 * e2;

    // according to
    // H. Vermeille,
    // Direct transformation from geocentric to geodetic ccordinates,
    // Journal of Geodesy (2002) 76:451-454
    Vector3d lla;
    double X = ecef_pos(0);
    double Y = ecef_pos(1);
    double Z = ecef_pos(2);
    double XXpYY = X * X + Y * Y;
    if (XXpYY + Z * Z < 25)
    {
        // This function fails near the geocenter region, so catch
        // that special case here.  Define the innermost sphere of
        // small radius as earth center and return the coordinates
        // 0/0/-EQURAD. It may be any other place on geoide's surface,
        // the Northpole, Hawaii or Wentorf. This one was easy to code
        // ;-)
        lla(0) = 0.0;
        lla(1) = 0.0;
        lla(2) = -EarthRadius;
        return lla;
    }

    double sqrtXXpYY = sqrt(XXpYY);
    double p = XXpYY * ra2;
    double q = Z * Z * (1 - e2) * ra2;
    double r = 1 / 6.0 * (p + q - e4);
    double s = e4 * p * q / (4 * r * r * r);
    /* 
       s*(2+s) is negative for s = [-2..0]
       slightly negative values for s due to floating point rounding errors
       cause nan for sqrt(s*(2+s))
       We can probably clamp the resulting parable to positive numbers
    */
    if (s >= -2.0 && s <= 0.0)
        s = 0.0;
    double t = pow(1 + s + sqrt(s * (2 + s)), 1 / 3.0);
    double u = r * (1 + t + 1 / t);
    double v = sqrt(u * u + e4 * q);
    double w = e2 * (u + v - q) / (2 * v);
    double k = sqrt(u + v + w * w) - w;
    double D = k * sqrtXXpYY / (k + e2);
    lla(1) = 2 * atan2(Y, X + sqrtXXpYY);
    double sqrtDDpZZ = sqrt(D * D + Z * Z);
    lla(0) = 2 * atan2(Z, D + sqrtDDpZZ);
    lla(2) = (k + e2 - 1) * sqrtDDpZZ / k;
    return lla;
}

// This function converts a vector in ecef to ned coordinate centered
// at pos_ref.
Vector3f ecef2ned(Vector3d ecef, Vector3d pos_ref)
{
    double lat = pos_ref(0, 0);
    double lon = pos_ref(1, 0);
    double sin_lat = sin(lat);
    double sin_lon = sin(lon);
    double cos_lat = cos(lat);
    double cos_lon = cos(lon);

    Vector3f ned;
    ned(2, 0) = -cos_lat * cos_lon * ecef(0, 0) - cos_lat * sin_lon * ecef(1, 0) - sin_lat * ecef(2, 0);
    ned(1, 0) = -sin_lon * ecef(0, 0) + cos_lon * ecef(1, 0);
    ned(0, 0) = -sin_lat * cos_lon * ecef(0, 0) - sin_lat * sin_lon * ecef(1, 0) + cos_lat * ecef(2, 0);

    return ned;
}

// Return a quaternion rotation from the earth centered to the
// simulation usual horizontal local frame from given longitude and
// latitude.  The horizontal local frame used in simulations is the
// frame with x-axis pointing north, the y-axis pointing eastwards and
// the z axis pointing downwards.  (Returns the ecef2ned
// transformation as a quaternion.)
Quaterniond lla2quat(double lon_rad, double lat_rad)
{
    Quaterniond q;
    double zd2 = 0.5 * lon_rad;
    double yd2 = -0.25 * M_PI - 0.5 * lat_rad;
    double Szd2 = sin(zd2);
    double Syd2 = sin(yd2);
    double Czd2 = cos(zd2);
    double Cyd2 = cos(yd2);
    q.w() = Czd2 * Cyd2;
    q.x() = -Szd2 * Syd2;
    q.y() = Czd2 * Syd2;
    q.z() = Szd2 * Cyd2;
    return q;
}

// This function gives a skew symmetric matrix from a given vector w
Matrix3f sk(Vector3f w)
{
    Matrix3f C;

    C(0, 0) = 0.0;
    C(0, 1) = -w(2, 0);
    C(0, 2) = w(1, 0);
    C(1, 0) = w(2, 0);
    C(1, 1) = 0.0;
    C(1, 2) = -w(0, 0);
    C(2, 0) = -w(1, 0);
    C(2, 1) = w(0, 0);
    C(2, 2) = 0.0;

    return C;
}

// This function gives a skew symmetric matrix from a given vector w; using double
Matrix3d sk2(Vector3d w)
{
    Matrix3d C;

    C(0, 0) = 0.0;
    C(0, 1) = -w(2, 0);
    C(0, 2) = w(1, 0);
    C(1, 0) = w(2, 0);
    C(1, 1) = 0.0;
    C(1, 2) = -w(0, 0);
    C(2, 0) = -w(1, 0);
    C(2, 1) = w(0, 0);
    C(2, 2) = 0.0;

    return C;
}


// Quaternion to euler angle: returns phi, the, psi as a vector
Vector3f quat2eul(Quaternionf q)
{
    float q0, q1, q2, q3;
    float m11, m12, m13, m23, m33;

    q0 = q.w();
    q1 = q.x();
    q2 = q.y();
    q3 = q.z();

    m11 = 2 * (q0 * q0 + q1 * q1) - 1;
    m12 = 2 * (q1 * q2 + q0 * q3);
    m13 = 2 * (q1 * q3 - q0 * q2);
    m23 = 2 * (q2 * q3 + q0 * q1);
    m33 = 2 * (q0 * q0 + q3 * q3) - 1;

    Vector3f result;
    result(2) = atan2(m12, m11);
    result(1) = asin(-m13);
    result(0) = atan2(m23, m33);

    return result;
}

// Computes a quaternion from the given euler angles
Quaternionf eul2quat(float phi, float the, float psi)
{
    float sin_psi = sin(psi * 0.5);
    float cos_psi = cos(psi * 0.5);
    float sin_the = sin(the * 0.5);
    float cos_the = cos(the * 0.5);
    float sin_phi = sin(phi * 0.5);
    float cos_phi = cos(phi * 0.5);

    Quaternionf q;
    q.w() = cos_psi * cos_the * cos_phi + sin_psi * sin_the * sin_phi;
    q.x() = cos_psi * cos_the * sin_phi - sin_psi * sin_the * cos_phi;
    q.y() = cos_psi * sin_the * cos_phi + sin_psi * cos_the * sin_phi;
    q.z() = sin_psi * cos_the * cos_phi - cos_psi * sin_the * sin_phi;

    return q;
}

// Quaternion to C_N2B
Matrix3f quat2dcm(Quaternionf q)
{
    float q0, q1, q2, q3;
    Matrix3f C_N2B;

    q0 = q.w();
    q1 = q.x();
    q2 = q.y();
    q3 = q.z();

    C_N2B(0, 0) = 2 * (q0 * q0 + q1 * q1) - 1;
    C_N2B(1, 1) = 2 * (q0 * q0 + q2 * q2) - 1;
    C_N2B(2, 2) = 2 * (q0 * q0 + q3 * q3) - 1;

    C_N2B(0, 1) = 2 * (q1 * q2 + q0 * q3);
    C_N2B(0, 2) = 2 * (q1 * q3 - q0 * q2);

    C_N2B(1, 0) = 2 * (q1 * q2 - q0 * q3);
    C_N2B(1, 2) = 2 * (q2 * q3 + q0 * q1);

    C_N2B(2, 0) = 2 * (q1 * q3 + q0 * q2);
    C_N2B(2, 1) = 2 * (q2 * q3 - q0 * q1);

    return C_N2B;
}

// EphemerisData (subframe1,2,3) to Satellite ecef x, y, z in meter, vx, vy, vz in m/s
VectorXd EphemerisData2Satecef(float t,
                               uint32_t TOW, uint8_t L2, uint16_t week_No, uint8_t L2_Flag, uint8_t SV_Acc, uint8_t SV_Hlth,
                               double T_GD, uint16_t IODC, double t_OC, int8_t a_f2, double a_f1, double a_f0,
                               uint8_t IODE, double C_rs, double delta_n, double M_0, double C_uc, double ecc, double C_us,
                               double sqrt_A, double t_OE, double C_ic, double Omega_0, double C_is, double i_0, double C_rc,
                               double omega, double Omega_dot, double IDOT)
{
    // All the equations are based ON: https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf
    // pg. 104-105, Also Grove  p335-338

    // Process subframe 1,2,3 information
    double A_semiMajorAxis;        // Semi-major axis
    double n_0_computedMeanMotion; // Computed mean motion
    double n_correctedMeanMotion;  // Corrected mean motion
    double e_eccentricity;         // Eccentricity
    //double phi_k_argumentOfLattitude;   // Argument of latitude
    double M_0_trueAnomalyAtRef;
    double omega0_longitudeofAscendingNodeofOrbitPlane;
    double omega_argumentOfPerigee;
    double omegaDot_argumentOfPerigee;
    double i_0_inclinationAtRef;
    double iDot_rateOfInclination;

    A_semiMajorAxis = pow(sqrt_A, 2);
    n_0_computedMeanMotion = sqrt(MU / pow(A_semiMajorAxis, 3));
    n_correctedMeanMotion = n_0_computedMeanMotion + delta_n;
    e_eccentricity = ecc;
    M_0_trueAnomalyAtRef = M_0;
    omega0_longitudeofAscendingNodeofOrbitPlane = Omega_0;
    omega_argumentOfPerigee = omega;
    omegaDot_argumentOfPerigee = Omega_dot;
    i_0_inclinationAtRef = i_0;
    iDot_rateOfInclination = IDOT;

    // Compute the time from the ephemeris reference epoch
    double t_k_timeFromReferenceEpoch = t - t_OE;
    // Correct that time for end-of-week crossovers
    if (t_k_timeFromReferenceEpoch > 302400)
    {
        t_k_timeFromReferenceEpoch -= 604800;
    }
    if (t_k_timeFromReferenceEpoch < -302400)
    {
        t_k_timeFromReferenceEpoch += 604800;
    }

    // Compute the mean anomaly
    double M_k_meanAnomaly = M_0_trueAnomalyAtRef + n_correctedMeanMotion * t_k_timeFromReferenceEpoch;

    // Below, we iteratively solve for E_k_eccentricAnomaly using Newton-Raphson method
    double solutionError = 1000000.;
    double E_k_eccentricAnomaly = 1.;
    double currentDerivative = 0;
    int iterationCount = 0;

    solutionError = (E_k_eccentricAnomaly -
                     (e_eccentricity * sin(E_k_eccentricAnomaly)) -
                     M_k_meanAnomaly);

    while ((fabs(solutionError) > 1.0e-6) &&
           iterationCount < 1000)
    {
        currentDerivative = (1.0 - (e_eccentricity * cos(E_k_eccentricAnomaly)));
        E_k_eccentricAnomaly = E_k_eccentricAnomaly - solutionError / currentDerivative;

        solutionError = (E_k_eccentricAnomaly -
                         (e_eccentricity * sin(E_k_eccentricAnomaly)) -
                         M_k_meanAnomaly);
        iterationCount += 1;
        //   if (VERBOSE)
        //   {
        //     std::cout<< "Iteration #: " << iterationCount << " Error: " << solutionError << std::endl;
        //   }
    }
    double cos_E_k = cos(E_k_eccentricAnomaly);
    double sin_E_k = sin(E_k_eccentricAnomaly);
    double nu_k_trueAnomaly = atan2(
        (sqrt(1.0 - pow(e_eccentricity, 2)) * sin_E_k) /
            (1.0 - (e_eccentricity * cos_E_k)),
        (cos_E_k - e_eccentricity) /
            (1.0 - e_eccentricity * cos_E_k));

    double phi_k_argumentOfLatitude = nu_k_trueAnomaly + omega_argumentOfPerigee;

    // Compute the corrective 2nd order terms
    double sin2PhiK = sin(2.0 * phi_k_argumentOfLatitude);
    double cos2PhiK = cos(2.0 * phi_k_argumentOfLatitude);

    double deltaU_argumentOfLatCorrection = (C_us * sin2PhiK) + (C_uc * cos2PhiK);
    double deltaR_radiusCorrection = (C_rs * sin2PhiK) + (C_rc * cos2PhiK);
    double deltaI_inclinationCorrection = (C_is * sin2PhiK) + (C_ic * cos2PhiK);

    // Now compute the updated corrected orbital elements
    double u_argumentOfLat = phi_k_argumentOfLatitude + deltaU_argumentOfLatCorrection;
    double r_radius = (A_semiMajorAxis * (1.0 - (e_eccentricity * cos_E_k))) + deltaR_radiusCorrection;
    double i_inclination =
        i_0_inclinationAtRef +
        (iDot_rateOfInclination * t_k_timeFromReferenceEpoch) +
        deltaI_inclinationCorrection;

    // Compute the satellite position within the orbital plane
    double xPositionOrbitalPlane = r_radius * cos(u_argumentOfLat);
    double yPositionOrbitalPlane = r_radius * sin(u_argumentOfLat);
    double omegaK_longitudeAscendingNode =
        omega0_longitudeofAscendingNodeofOrbitPlane +
        ((omegaDot_argumentOfPerigee - OMEGA_DOT_EARTH) * t_k_timeFromReferenceEpoch) -
        (OMEGA_DOT_EARTH * t_OE);

    double sinOmegaK = sin(omegaK_longitudeAscendingNode);
    double cosOmegaK = cos(omegaK_longitudeAscendingNode);

    double sinIK = sin(i_inclination);
    double cosIK = cos(i_inclination);
    // Earth-fixed coordinates:
    double x = (xPositionOrbitalPlane * cosOmegaK) - (yPositionOrbitalPlane * cosIK * sinOmegaK);
    double y = (xPositionOrbitalPlane * sinOmegaK) + (yPositionOrbitalPlane * cosIK * cosOmegaK);
    double z = (yPositionOrbitalPlane * sinIK);


    
    // ECEF velocity calculation:
    double E_dot_k_eccentricAnomaly = n_correctedMeanMotion/(1.0 - (e_eccentricity * cos_E_k)); // Eq.(8.21)
    double phi_dot_k_argumentOfLatitude = sin(nu_k_trueAnomaly)/sin_E_k*E_dot_k_eccentricAnomaly; // Eq.(8.22)
    double r_dot_o_os = (A_semiMajorAxis*e_eccentricity*sin_E_k)*E_dot_k_eccentricAnomaly +
                        2*((C_rs * cos2PhiK) - (C_rc * sin2PhiK))*phi_dot_k_argumentOfLatitude; // Eq.(8.23a)
    double u_dot_o_os = (1+2*C_us * cos2PhiK - 2*C_uc * sin2PhiK)*phi_dot_k_argumentOfLatitude; // Eq.(8.23b)                    
     
    double x_dot_o_os = r_dot_o_os*cos(u_argumentOfLat) - r_radius*u_dot_o_os*sin(u_argumentOfLat); // Eq.(8.24a)
    double y_dot_o_os = r_dot_o_os*sin(u_argumentOfLat) + r_radius*u_dot_o_os*cos(u_argumentOfLat); // Eq.(8.24b)

    double omega_dot_K_longitudeAscendingNode = omegaDot_argumentOfPerigee - OMEGA_DOT_EARTH; // Eq. (8.25)
    double i_dot_inclination = iDot_rateOfInclination + 2*((C_is * cos2PhiK) - (C_ic * sin2PhiK))*phi_dot_k_argumentOfLatitude; // Eq. (8.26)
    
    // Eq. (8.27)
    double vx = (x_dot_o_os*cosOmegaK - y_dot_o_os*cosIK*sinOmegaK + i_dot_inclination*yPositionOrbitalPlane*sinIK*sinOmegaK) -
                omega_dot_K_longitudeAscendingNode*(xPositionOrbitalPlane*sinOmegaK + yPositionOrbitalPlane*cosIK*cosOmegaK);
    double vy = (x_dot_o_os*sinOmegaK + y_dot_o_os*cosIK*cosOmegaK - i_dot_inclination*yPositionOrbitalPlane*sinIK*cosOmegaK) -
                omega_dot_K_longitudeAscendingNode*(-xPositionOrbitalPlane*cosOmegaK +yPositionOrbitalPlane*cosIK*sinOmegaK);
    double vz = (y_dot_o_os*sinIK + i_dot_inclination*yPositionOrbitalPlane*cosIK);
    
    VectorXd pos_vel_Sat_ecef(6);
    pos_vel_Sat_ecef << x, y, z, vx, vy, vz;
    // cout << x << y << y << endl;


    return pos_vel_Sat_ecef;
}

// Compute direction cosine matric C from a Euler vector 
// eul = [yaw,pitch,roll]. (i.e., 3-2-1 rotation convention)
Matrix3d eul2dcm(Vector3d eul)
{
    Matrix3d C, C1, C2, C3;
    double ps=eul(0);
    double th=eul(1); 
    double ph=eul(2);
     
    C1(0, 0) = 1;       C1(0, 1) = 0;         C1(0, 2) = 0;
    C1(1, 0) = 0;       C1(1, 1) = cos(ph);   C1(1, 2) = sin(ph);
    C1(2, 0) = 0;       C1(2, 1) = -sin(ph);  C1(2, 2) = cos(ph);
   
    C2(0, 0) = cos(th);       C2(0, 1) = 0;   C2(0, 2) = -sin(th);
    C2(1, 0) = 0;             C2(1, 1) = 1;   C2(1, 2) = 0;
    C2(2, 0) = sin(th);       C2(2, 1) = 0;   C2(2, 2) = cos(th);
    
    C3(0, 0) = cos(ps);   C3(0, 1) = sin(ps);   C3(0, 2) = 0;
    C3(1, 0) = -sin(ps);  C3(1, 1) = cos(ps);   C3(1, 2) = 0;
    C3(2, 0) = 0;         C3(2, 1) = 0;         C3(2, 2) = 1;
   
    C=C1*C2*C3;
    return C;
}


// ned 2 ecef centered at the coordinate given by lla
Vector3d ned2ecef(Vector3d ned, Vector3d lla)
{
    double lat = lla(0, 0);
    double lon = lla(1, 0);
    double pitch = abs(lat) + M_PI/2;
    Vector3d eul;
    if (lat >=0)
    {
        eul(0,0) = lon; 
        eul(1,0) = -pitch;
        eul(2,0) = 0;
    }
    else
    {
        eul(0,0) = lon; 
        eul(1,0) = pitch;
        eul(2,0) = 0;
    }

    Matrix3d C_ecef2ned = eul2dcm(eul);
    Vector3d ecef  = C_ecef2ned.transpose()*ned;


    return ecef;
}
