#include <iostream>
#include <array>
using namespace std;

class System{
  public:
    int time;
};

class Satellite{
  public:
    array<int,2> pos;
    array<int, 2> vel;
    int ang_mom;
};

class Orbit{
  public:
    int a;
    int e;
    int ang_mom;
};


int main() {
  Satellite cubeSat;
  cubeSat.pos = {1, 2};
  cout << "Position X is " << cubeSat.pos[0];
}