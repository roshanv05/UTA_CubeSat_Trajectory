/*
Sun Synchronous Transfer Model (2D)
Author: Karthik Nuti
University of Texas- Arlington AIAA CubeSat
Date: 12/05/2025

Purpose: 
    The following code is to calculate [will add the summary]



Units:
    - All units are in SI (meters, seconds, radians)
    - Distance: KM
    - Time: s
    - Velocity: km/s
    - Mass: kg
    - Isp: s
    - g0: km/s^2
*/
#include <iostream>
#include <cmath>
#include <array>
#include <limits>

//outputs and overwrites the file
#include <iomanip>
#include <fstream>

using std::cout; using std::cin; using std::endl;

// ===================== CONSTANTS ===================== //
const double MU_EARTH = 398600.4418;        // km^3/s^2
const double R_EARTH = 6378.1363;           // km
const double J2_EARTH = 1.082638e-3;        // Earth's J2 coefficient
const double PI = 3.14159265358979323846;   // Pi
const double DEG2RAD = PI / 180.0;          // Degrees to Radians conversion
const double RAD2DEG = 180.0 / PI;          // Radians to Degrees conversion
const double g0 = 9.80665/1000.0;            // km/s^2  

// ===================== Basic STRUCTURES ===================== //
class System{
    public:
        double time; // simulation time [s]
};

class Orbit {
    public:
        double a;          // semi-major axis [km]
        double e;          // eccentricity [-]
        double inc;        // inclination [rad] (for sso, computed from J2)
        double ang_mom;    // specific angular momentum [km^2/s]

        Orbit() : a(0.0), e(0.0), inc(0.0), ang_mom(0.0) {}

        double rPeri() const { return a * (1.0-e); } // periapsis radius
        double rApop() const { return a * (1.0+e); } // apoapsis radius
};

class Satellite{
    public:
        std::array<double,2> pos;       // 2D Cart. pos [km]
        std::array<double,2> vel;       // 2D Cart. vel [km/s]  
        double ang_mom;
        double mass0;                   // initial mass [kg]

        Satellite () : ang_mom(0.0), mass0(0.0) {
            pos = {0.0, 0.0};
            vel = {0.0, 0.0};
        }
};



class TransferOrbit : public Orbit {
    public:
        double dv1;        // first impulse [km/s] 
        double dv2;        // second impulse [km/s]
        double dv_tot;     // total impulse [km/s]

        TransferOrbit() : Orbit(), dv1(0.0), dv2(0.0), dv_tot(0.0) {}
};


// ===================== ASTRODYNAMICS ===================== //

//vis-viva equation: v = sqrt(mu*(2/r - 1/a))
double visViva(double r, double a) {
    return std::sqrt(MU_EARTH * (2.0 / r - 1.0 / a));
}

// compute specific angular momentum for a Keplerian orbit
double specificAngularMomentum (double a, double e) {
    return std::sqrt(MU_EARTH * a * (1.0-e*e));
}

/*
- J2 NODAL Precession rate for given orbit [rad/s]
- rho = domega/dt = -(3/2) * J2 * (R_EARTH^2)/ p^2) * sqrt(mu/p) *cos(i)
- for the sun-synchronous orbit, we will us eht equivalent form in terms 
of a, e and solve for cos(i)
*/

double sunSyncInc (double a_km, double e) {
    //target nodal precession: app. 360 deg/sidereal year
    const double day_sec = 86400.0;
    const double year_sec = 365.2422 *day_sec;
    const double rho_target = 2.0 * PI / year_sec; // [rad/s]

    //semi-latus rectum p = a(1-e^2)
    double p = a_km * (1.0 -e*e);

    /* 
    rho ≈ -(3*J2*R_E^2*sqrt(mu)*cos(i)) / (2 * a^(7/2) * (1 - e^2)^2)
    so, cos(i) = -2 * rho * a^(7/2) * (1 - e^2)^2 / (3 * J2 * R_E^2 * sqrt(mu))
    */

    double a_72 = std::pow(a_km, 3.5); // a^(7/2)
    double denom = 3.0 * J2_EARTH * R_EARTH * R_EARTH * std::sqrt(MU_EARTH);
    double numer = -2.0 * rho_target *a_72 * (1.0-e*e) * (1.0-e*e);

    double cos_i = numer/denom;


    // clamp cos(i) to [-1,1] to avoid NaN
    if (cos_i > 1.0) cos_i = 1.0;
    if (cos_i < -1.0) cos_i = -1.0;

    return std::acos(cos_i); //return inclination in radians
}

Orbit buildOrbit(double alt_km, double e, bool makeSunSync) {
    Orbit orb;
    orb.e = e;
    orb.a = R_EARTH + alt_km; //semi-major axis [km]
    orb.ang_mom = specificAngularMomentum(orb.a, orb.e);

    if (makeSunSync){
        orb.inc = sunSyncInc(orb.a, orb.e);
    } else {
        orb.inc = 0.0; //equatorial orbit  
    }
    return orb;
}


/*
Build Hohmann transfer orbit between two coplanar tangent radii/
Assumptions:
   - first impulse at periapsis radius of LEO orbit
   - second impulse at periapsis radius of SSO orbit
*/

TransferOrbit buildTransferOrbit(const Orbit& leo, const Orbit& sso) {
    TransferOrbit t;

    //tangent radii (use both periapsis radii for both orbits for simplicity for now)
    double r1 = leo.rPeri();
    double r2 = sso.rPeri();

    //semi major acis of transfer orbit
    t.a = 0.5 * (r1+r2);
    t.e = (r2-r1) / (r1+r2); // eccentricity of hohmann ellipse
    t.ang_mom = specificAngularMomentum(t.a, t.e);
    t.inc = leo.inc; //coplanar assumtion

    // Velocities at first tanget
    double v_leo = visViva(r1, leo.a);
    double v_trans1 = visViva(r1, t.a);
    t.dv1 = std::fabs(v_trans1 - v_leo);

    // Velocities at second tangent
    double v_sso = visViva(r2, sso.a);
    double v_trans2 = visViva(r2, t.a);
    t.dv2 = std::fabs(v_sso - v_trans2);

    t.dv_tot = t.dv1 + t.dv2;

    return t;
}

/*
Rocket equation
Solve for Mf and prop mass
*/

void rocketEq(double deltaV_km_s, double Isp_s, double m0_kg, double& mf_kg, double& mprop_kg) {
    if (Isp_s <= 0.0 || m0_kg <= 0.0){
        mf_kg = 0.0;
        mprop_kg = 0.0;
        return;
    }

    double exponent = deltaV_km_s / (Isp_s * g0);
    mf_kg = m0_kg / std::exp(exponent);
    mprop_kg = m0_kg - mf_kg;
}

// ===================== MAIN ===================== //

int main() {
    std::cout << "==== 2D Sun-Synchronous Transfer Model (C++) ====\n\n";

    double h_leo, e_leo;
    double h_sso, e_sso;
    double Isp, m0;

    std::cout << "Enter LEO altitude above Earth [km]: ";
    std::cin >> h_leo;

    std::cout << "Enter LEO eccentricity (0 ~ near-circular): ";
    std::cin >> e_leo;

    std::cout << "Enter SSO altitude above Earth [km]: ";
    std::cin >> h_sso;

    std::cout << "Enter SSO eccentricity (0 ~ near-circular): ";
    std::cin >> e_sso;

    std::cout << "Enter engine specific impulse Isp [s]: ";
    std::cin >> Isp;

    std::cout << "Enter initial spacecraft mass m0 [kg]: ";
    std::cin >> m0;

    // ---- Build orbits ----
    Orbit leo  = buildOrbit(h_leo,  e_leo,  false);
    Orbit sso  = buildOrbit(h_sso,  e_sso,  true);

    // ---- Transfer orbit ----
    TransferOrbit trans = buildTransferOrbit(leo, sso);

    // ---- Rocket equation ----
    double mf, mprop;
    rocketEq(trans.dv_tot, Isp, m0, mf, mprop);
    double Itot = mprop * Isp * (g0 * 1000.0); // N·s

    
    // ---- Console output (for you) ----
    cout << "\n=== ORBIT PARAMETERS ===\n";
    cout << "LEO a [km]: " << leo.a << "  e: " << leo.e << "\n";
    cout << "SSO a [km]: " << sso.a << "  e: " << sso.e << "  i_sso [deg]: " << sso.inc * RAD2DEG << "\n";
    cout << "Transfer a [km]: " << trans.a << "  e: " << trans.e << "\n";
    cout << "dv1 [km/s]: " << trans.dv1 << "  dv2 [km/s]: " << trans.dv2 << "  dv_total [km/s]: " << trans.dv_tot << "\n";
    cout << "m0 [kg]: " << m0 << "  mf [kg]: " << mf << "  m_prop [kg]: " << mprop << "\n";
    cout << "Total impulse [N·s]: " << Itot << "\n";

    // ---- Console output (sanity check) ----
    std::cout << "\n=== SUMMARY ===\n";
    std::cout << "LEO a [km]: " << leo.a << "  e: " << leo.e << "\n";
    std::cout << "SSO a [km]: " << sso.a << "  e: " << sso.e
              << "  i_sso [deg]: " << sso.inc * RAD2DEG << "\n";
    std::cout << "Transfer a [km]: " << trans.a << "  e: " << trans.e << "\n";
    std::cout << "dv_total [km/s]: " << trans.dv_tot << "\n";
    std::cout << "m0 [kg]: " << m0 << "  mf [kg]: " << mf
              << "  m_prop [kg]: " << mprop << "\n";
    std::cout << "Total impulse [N·s]: " << Itot << "\n";

    // ---- CSV FILE OUTPUT ----
    std::ofstream fout("orbit_results.csv", std::ios::trunc);
    if (!fout.is_open()) {
        std::cerr << "Error: could not open orbit_results.csv for writing.\n";
        return 1;
    }

    fout << std::fixed << std::setprecision(8);

    // Header
    fout << "a_leo_km,e_leo,r_peri_leo_km,r_apo_leo_km,"
         << "a_sso_km,e_sso,inc_sso_deg,r_peri_sso_km,r_apo_sso_km,"
         << "a_trans_km,e_trans,r_peri_trans_km,r_apo_trans_km,"
         << "dv1_kms,dv2_kms,dv_total_kms,"
         << "m0_kg,mf_kg,mprop_kg,Itot_Ns\n";

    // Row
    fout << leo.a << "," << leo.e << "," << leo.rPeri() << "," << leo.rApop() << ","
         << sso.a << "," << sso.e << "," << (sso.inc * RAD2DEG) << "," << sso.rPeri() << "," << sso.rApop() << ","
         << trans.a << "," << trans.e << "," << trans.rPeri() << "," << trans.rApop() << ","
         << trans.dv1 << "," << trans.dv2 << "," << trans.dv_tot << ","
         << m0 << "," << mf << "," << mprop << "," << Itot << "\n";

    fout.close();

    std::cout << "\n[OK] Results written to orbit_results.csv in this folder.\n";
    std::cout << "\nPress Enter to exit...";
    std::cin.ignore();
    std::cin.get();
    return 0;
}
