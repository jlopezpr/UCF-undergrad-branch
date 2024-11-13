#ifndef GLOBAL_H
#define GLOBAL_H

#include <vector>
#include <random>
#include "vec3.h"

const double pi = 3.141592;
const int N = 500;
const double phi = 0.01;
const int timeSteps = 100000;
const double dt = 0.0001;
const int writingRate = 100;
const double maxOverlap = 100;

struct global{

  double L; // length of box
  int timeStep = 0;

  std::vector<vec3> pos;
  std::vector<vec3> forces;

  std::mt19937 gen;
  std::normal_distribution<> gauss;
};

#endif