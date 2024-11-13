#include "simulation.h"

int main(){
  global g;
  hardsphere h;
  asakuraoosawa a;
  yukawa y;

  std::cout << "Simulating...\n";

  simulation(g, h, a, y);
  std::cout << "Simulation complete.\n";
  return 0;
}