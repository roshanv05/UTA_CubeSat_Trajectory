// Libraries
#include <iostream>
#include <array>
using namespace std;

// Variables
double earthMass = 5.972e24;
double cubeSatMass = 1.33;

// Constants
const double G = 6.6743e-11;

class System{
  public:
    int t; // Time

    System(){
      t = 0;
    }
};

class Body{ //Assume fixed in space
  public:
    double m; // Mass
    array<int, 2> pos; // Position

    Body(double bodyMass){
      m = bodyMass;
    }
};

class Satellite{
  public:
    double m; // Mass
    array<int,2> pos; // Position
    array<int, 2> vel; // Velocity
    int ang_mom;

    Satellite(double satelliteMass){
      m = satelliteMass;
    }
};

class Orbit{
  public:
    int a; // Semi-Major Axis
    int e; // Eccentricity
    int ang_mom; // Angular Momentum
    double mu;

    Orbit(double satelliteMass, double bodyMass, double G){
      mu = G*(satelliteMass+bodyMass);
    }
};

int main() {
  // Create Objects
  System system;

  Body earth( earthMass);
  Satellite cubeSat(cubeSatMass);

  Orbit orbit(earthMass, cubeSatMass, G);

  cout << "Gravitational Parameter is " << orbit.mu;

  return 0;
}