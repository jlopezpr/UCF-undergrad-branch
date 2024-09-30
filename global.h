#ifndef GLOBAL_H
#define GLOBAL_H

#include<cmath> 
#include<vector>
#include "vec3.h"

constexpr double Pi = 3.14159265358979323846;

constexpr int N = 1000; // numero de particulas
constexpr double phi = 0.01;
constexpr double Dt = 0.00001;
constexpr int numberOfTimeSteps = 10000;

constexpr int maxOverlapp = 100;
constexpr int writtingRate = 100;

struct global{

    double L; // longitud de la caja de simulacion
    int timeStep = 0;

    std::vector<vec3> pos;
};

#endif