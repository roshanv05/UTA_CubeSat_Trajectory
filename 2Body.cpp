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

double earthMass = 5.972e24;
double earthRadius = 6378;

double cubeSatMass = 1.33;

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

    Body(double bodyMass){
      m = bodyMass;
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

    Satellite(double satelliteMass, double altitude, double earthRadius){
      m = satelliteMass;
      pos << 0, altitude + earthRadius;
      r = pos.norm();
    }

    void solveAcceleration(double mu){
      accel = -mu * pos/pow(r, 3);
    };

    void solveEOM(double t, double dt){
      for (int i = 0; i < t; i++) {
        t = i * dt;
        vel[i+1] = accel[i] * t + vel[i];
        pos[i+1] = vel[i] * t + pos[i];
        };
      };
    
  };

class Orbit{
  public:
    double a; // Semi-Major Axis
    double e; // Eccentricity
    double ang_mom; // Angular Momentum
    double mu; // Gravitational Parameter

    Orbit(double satelliteMass, double bodyMass, double G){
      mu = G*(satelliteMass+bodyMass);
    }

};

int main() {

  System system;
  Body earth(earthMass);
  Satellite cubeSat(cubeSatMass, z, earthRadius);
  Orbit orbit(earthMass, cubeSatMass, G);

  cubeSat.solveAcceleration(orbit.mu);
  cout << "Accel:" << cubeSat.accel;

  return 0;
}