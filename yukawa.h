#ifndef YUKAWA_H
#define YUKAWA_H

#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include "vec3.h"
#include "global.h"
#include "init.h"

struct yukawa {

    double cutoff = 25.0;
    const double A = 8.0;
    const double kappa = 0.8;

    std::vector<std::vector<int> > neighborCells;
    std::vector<std::vector<int> > neighborList;
    
    std::vector<vec3> force;

    std::mt19937 gen;
    std::normal_distribution<> gauss;
};

void yukawaNeighborCells(yukawa& y, global& g){

    const int cellsPerDim = g.L / y.cutoff;
    const double cellSize = g.L / cellsPerDim;
    const int numberOfCells = cellsPerDim * cellsPerDim * cellsPerDim;

    y.neighborCells.resize(numberOfCells);

    for (int i = 0; i < cellsPerDim; i++){
        for (int j = 0; j < cellsPerDim; j++){
            for (int k = 0; k < cellsPerDim; k++){ 
                int label = k + j * cellsPerDim + i * cellsPerDim * cellsPerDim;

                std::vector<int> neighbors;
                
                for (int x = -1; x < 2; ++x){
                    for (int y = -1; y < 2; ++y){
                        for (int z = -1; z < 2; ++z){
                            int nx = i + x;
                            int ny = j + y;
                            int nz = k + z;

                            nx += (nx < 0) * cellsPerDim - (nx == cellsPerDim) * cellsPerDim;
                            ny += (ny < 0) * cellsPerDim - (ny == cellsPerDim) * cellsPerDim;
                            nz += (nz < 0) * cellsPerDim - (nz == cellsPerDim) * cellsPerDim;

                            int idx = nz + ny * cellsPerDim + nx * cellsPerDim * cellsPerDim;
                            neighbors.push_back(idx);
                        }
                    }
                }

                y.neighborCells[label] = neighbors;
            }
        }
    }
}

void yukawaRandom(yukawa& y){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::normal_distribution<> gauss(0, 1);

    y.gen = gen;
    y.gauss = gauss;
}

void yukawaNeighborList(yukawa& k, global& g){
    
    const int cellsPerDim = g.L / k.cutoff;
    const double cellSize = g.L / cellsPerDim;
    const int numberOfCells = cellsPerDim * cellsPerDim * cellsPerDim;

    std::vector<std::vector<int> > verletList(numberOfCells);

    for (int i = 0; i < N; ++ i){
        double x = g.pos[i].x;
        double y = g.pos[i].y;
        double z = g.pos[i].z;

        int nx = (x / cellSize);
        int ny = (y / cellSize);
        int nz = (z / cellSize);

        nx += (nx < 0) * 1 - (nx == cellsPerDim) * 1;
        ny += (ny < 0) * 1 - (ny == cellsPerDim) * 1;
        nz += (nz < 0) * 1 - (nz == cellsPerDim) * 1;

        int label = (nx * cellsPerDim * cellsPerDim) + (ny * cellsPerDim) + nz;
        verletList[label].push_back(i);
    }
    

    double squareCut = k.cutoff * k.cutoff;

    for (int s = 0; s < numberOfCells; ++s){
        for (int i : verletList[s]){

            double x = g.pos[i].x;
            double y = g.pos[i].y;
            double z = g.pos[i].z;

            for (int t : k.neighborCells[s]){
                for (int j : verletList[t]){

                    if (i >= j) continue;
                    double x1 = x - g.pos[j].x;
                    double y1 = y - g.pos[j].y;
                    double z1 = z - g.pos[j].z;

                    double dist = x1 * x1 + y1 * y1 + z1 * z1;
                    if (dist > squareCut) continue;
                    k.neighborList[i].push_back(j);
                }
            }
        }
    }
}

void yukawaForces(yukawa& y, global& g){

    for (int i = 0; i < N; ++i){
        y.force[i].x = 0;
        y.force[i].y = 0;
        y.force[i].z = 0;
    }
    
    double A = y.A;
    double kappa = y.kappa;
    

    for (int i = 0; i < N; i++){
        std::vector<int> neighbors = y.neighborList[i];
        for (int j : neighbors){

            double dx = g.pos[i].x - g.pos[j].x;
            double dy = g.pos[i].y - g.pos[j].y;
            double dz = g.pos[i].z - g.pos[j].z;

            vec3 r = {dx, dy, dz};
            double r2 = r.squaremag();
            double r1 = std::sqrt(r2);
            double r_2 = 1.0 / r2;
            double r_1 = 1.0 / r1;
            double r_3 = r_2 * r_1;

            double expo = std::exp(-kappa * r1);
            double value = -A * r_3;
            value -= kappa * A * r_2; 
            value *= expo; 

            vec3 force = value * r;

            y.force[i] = y.force[i] + force;
            y.force[j] = y.force[j] - force;
            
        }
        
    }
}

void YKinit(yukawa& y, global& g){
    y.force.resize(N);
    y.neighborList.resize(N);

    yukawaNeighborCells(y, g);
}




#endif