// Libraries
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Variables
double e = 0;
double z = 500;

double dt = 1;

double earth_mass = 5.972e24;
double earth_radius = 6378;

double cubeSat_mass = 1.33;

// Constants
const double G = 6.6743e-20;

class System{
  public:
    int t; // Time

    System(){
      t = 0;
    }
};

class Body{ //Assume fixed in space
  public:
    double m;
    Vector2d pos;

    Body(double body_mass){
      m = body_mass;
      pos.setZero();
    }
};

class Satellite{
  public:
    double m;
    Vector2d pos; 
    Vector2d vel; 
    Vector2d accel; 

    double r;
    double true_anom;

    double current_mu;

    Satellite(double satellite_mass, double periapsis_altitude, double earth_radius ){
      m = satellite_mass;
      pos << periapsis_altitude + earth_radius, 0;
      r = solveRadius(pos);
    }

    double solveRadius(Vector2d pos){
      r = pos.norm();
      return r;
    }

    Vector2d solveVelocity(double mu, double h, Vector2d y){
      r = solveRadius(y);
      true_anom = atan2(y.y(), y.x());
      double v_rad = (mu/h)*e*sin(true_anom);
      double v_tan = (mu/h)*(1+e*cos(true_anom));
      vel << v_rad*cos(true_anom) - v_tan*sin(true_anom), v_rad*sin(true_anom) - v_tan*cos(true_anom);
      cout << vel.x() << "," << vel.y() << "\n";
      return vel;
    }

    Vector2d solveAccel(double mu, Vector2d y){
      r = solveRadius(y);
      Vector2d ddydt = -mu * y/pow(r, 3);
      current_mu = mu;
      return ddydt;
    };

    void propagate(double dt) {
      // K1
      Vector2d K1_v = vel;
      Vector2d K1_a = solveAccel(current_mu, pos);

      // K2
      Vector2d K2_v = vel + 0.5 * dt * K1_a;
      Vector2d K2_a = solveAccel(current_mu, pos + 0.5 * dt * K1_v);

      // K3
      Vector2d K3_v = vel + 0.5 * dt * K2_a;
      Vector2d K3_a = solveAccel(current_mu, pos + 0.5 * dt * K2_v);

      // K4
      Vector2d K4_v = vel + dt * K3_a;
      Vector2d K4_a = solveAccel(current_mu, pos + dt * K3_v);

      // Update position and velocity
      pos += (dt / 6.0) * (K1_v + 2*K2_v + 2*K3_v + K4_v);
      vel += (dt / 6.0) * (K1_a + 2*K2_a + 2*K3_a + K4_a);
    }
};

class Orbit{
  public:
    double a; // Semi-Major Axis
    double e; // Eccentricity
    double h; // Angular Momentum
    double mu; // Gravitational Parameter

    Orbit(double satellite_mass, double body_mass, double G, double periapsis_radius, double e){
      mu = G*(satellite_mass+body_mass);
      a = periapsis_radius/(1-e);
      h = sqrt(periapsis_radius*mu*(1+e));
    }

};

int main() {

  System system;
  Body earth(earth_mass);
  Satellite cubeSat(cubeSat_mass, z, earth_radius);
  Orbit orbit(cubeSat_mass, earth_mass, G, cubeSat.pos.norm(), e);

  cubeSat.current_mu = orbit.mu;
  cubeSat.vel = cubeSat.solveVelocity(orbit.mu, orbit.h, cubeSat.pos);
  int i = 0;
  while (i < 500) {
    cubeSat.propagate(1);
    cout << "Position:" << cubeSat.pos << "\n";
    i++;
  }

  return 0;
}