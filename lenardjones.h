#ifndef LENARDJONES_H
#define LENARDJONES_H
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include "vec3.h"
#include "global.h"
#include "init.h"

struct lenardjones {

    double cutoff = 4.0;
    const double epsilon = 100;
    const double sigma = 2.1;

    std::vector<std::vector<int> > neighborCells;
    std::vector<std::vector<int> > neighborList;
    
    std::vector<vec3> force;

    std::mt19937 gen;
    std::normal_distribution<> gauss;
};

void generateNeighborCells(lenardjones& l, global& g){

    const int cellsPerDim = g.L / l.cutoff;
    const double cellSize = g.L / cellsPerDim;
    const int numberOfCells = cellsPerDim * cellsPerDim * cellsPerDim;

    l.neighborCells.resize(numberOfCells);

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

                l.neighborCells[label] = neighbors;
            }
        }
    }
}

void generateRandom(lenardjones& l){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::normal_distribution<> gauss(0, 1);

    l.gen = gen;
    l.gauss = gauss;
}

void generateNeighborList(lenardjones& l, global& g){
    
    const int cellsPerDim = g.L / l.cutoff;
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
    

    double squareCut = l.cutoff * l.cutoff;

    for (int s = 0; s < numberOfCells; ++s){
        for (int i : verletList[s]){

            double x = g.pos[i].x;
            double y = g.pos[i].y;
            double z = g.pos[i].z;

            for (int t : l.neighborCells[s]){
                for (int j : verletList[t]){

                    if (i >= j) continue;
                    double x1 = x - g.pos[j].x;
                    double y1 = y - g.pos[j].y;
                    double z1 = z - g.pos[j].z;

                    double dist = x1 * x1 + y1 * y1 + z1 * z1;
                    if (dist > squareCut) continue;
                    l.neighborList[i].push_back(j);
                }
            }
        }
    }
}

void generateForces(lenardjones& l, global& g){

    for (int i = 0; i < N; ++i){
        l.force[i].x = 0;
        l.force[i].y = 0;
        l.force[i].z = 0;
    }
    
    double epsilon = l.epsilon;
    double sigma = l.sigma;
    double sigma6 = std::pow(sigma, 6);
    epsilon *= sigma6;

    for (int i = 0; i < N; i++){
        std::vector<int> neighbors = l.neighborList[i];
        for (int j : neighbors){

            double dx = g.pos[i].x - g.pos[j].x;
            double dy = g.pos[i].y - g.pos[j].y;
            double dz = g.pos[i].z - g.pos[j].z;

            vec3 r = {dx, dy, dz};
            double r2 = r.squaremag();
            double rsqrt = std::sqrt(r2);
            if (r2 > 4.0) continue;
            double r6 = std::pow(r2,3);
            double r8 = r6 * r2;
            
            double outFactor = (24 * epsilon * sigma6) / (r8);
            double inFactor = ((2 * sigma6) / r6) - 1;

            vec3 force = (inFactor * outFactor) * r;
            l.force[i] = l.force[i] - force;
            l.force[j] = l.force[j] + force
        }
        
    }
}

void periodicLC(global& g){

    for (int i = 0; i < N; ++i){
        g.pos[i].x += (g.pos[i].x < 0) * g.L - (g.pos[i].x >= g.L) * g.L;
        g.pos[i].y += (g.pos[i].y < 0) * g.L - (g.pos[i].y >= g.L) * g.L;
        g.pos[i].z += (g.pos[i].z < 0) * g.L - (g.pos[i].z >= g.L) * g.L;
    }
}

void LJinit(lenardjones& l, global& g){
    l.force.resize(N);
    l.neighborList.resize(N);

    generateNeighborCells(l, g);
}

#endif