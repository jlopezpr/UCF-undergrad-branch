#ifndef ASAKURAOOSAWA_H
#define ASAKURAOOSAWA_H

#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include "vec3.h"
#include "global.h"
#include "init.h"
#include "lenardjones.h"

struct asakuraoosawa {
    double cutoff = 4.0;
    const double epsilon = 100;
    const double sigma = 2.1;

    std::vector<std::vector<int> > neighborCells;
    std::vector<std::vector<int> > neighborList;
    
    std::vector<vec3> force;

    std::mt19937 gen;
    std::normal_distribution<> gauss;
};

void asakuraoosawaNeighborCells(asakuraoosawa& a, global& g){

    const int cellsPerDim = g.L / a.cutoff;
    const double cellSize = g.L / cellsPerDim;
    const int numberOfCells = cellsPerDim * cellsPerDim * cellsPerDim;

    a.neighborCells.resize(numberOfCells);

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

                a.neighborCells[label] = neighbors;
            }
        }
    }
}

void asakuraoosawaRandom(asakuraoosawa& a){
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::normal_distribution<> gauss(0, 1);

    a.gen = gen;
    a.gauss = gauss;
}

void asakuraoosawaNeighborList(asakuraoosawa& a, global& g){
    
    const int cellsPerDim = g.L / a.cutoff;
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
    

    double squareCut = a.cutoff * a.cutoff;

    for (int s = 0; s < numberOfCells; ++s){
        for (int i : verletList[s]){

            double x = g.pos[i].x;
            double y = g.pos[i].y;
            double z = g.pos[i].z;

            for (int t : a.neighborCells[s]){
                for (int j : verletList[t]){

                    if (i >= j) continue;
                    double x1 = x - g.pos[j].x;
                    double y1 = y - g.pos[j].y;
                    double z1 = z - g.pos[j].z;

                    double dist = x1 * x1 + y1 * y1 + z1 * z1;
                    if (dist > squareCut) continue;
                    a.neighborList[i].push_back(j);
                }
            }
        }
    }
}

void asakuraoosawaForces(asakuraoosawa& a, global& g){

    for (int i = 0; i < N; ++i){
        a.force[i].x = 0;
        a.force[i].y = 0;
        a.force[i].z = 0;
    }
    
    double A = 1.0; // hydrodynamic radius of the particles
    double delta = 0.1 * A + A; 
    double Ua = -3.0; // proportional to polymer concentration
    double A2 = A * A;
    double A3 = A2 * A;
    double delta2 = delta * delta;
    double delta3 = delta2 * delta;

    double denominator = 4 * A3 * delta3 - 6 * A3 * delta3 + 2 * A3;

    for (int i = 0; i < N; i++){
        std::vector<int> neighbors = a.neighborList[i];
        for (int j : neighbors){
            double dx = g.pos[i].x - g.pos[j].x;
            double dy = g.pos[i].y - g.pos[j].y;
            double dz = g.pos[i].z - g.pos[j].z;

            vec3 r = {dx, dy, dz};
            double r2 = r.squaremag();

            if (r2 < 4.0 || r2 > delta2) continue;

            double r1 = std::sqrt(r2);
            double numerator = 12 * A2 * delta2 - 3.0 * r2;

            double inFactor = numerator / denominator; 
            double outFactor = Ua / r1;

            vec3 force = (inFactor * outFactor) * r;
            a.force[i] = a.force[i] - force;
            a.force[j] = a.force[j] + force;
        }
        
    }
}

void periodicAC(global& g){

    for (int i = 0; i < N; ++i){
        g.pos[i].x += (g.pos[i].x < 0) * g.L - (g.pos[i].x >= g.L) * g.L;
        g.pos[i].y += (g.pos[i].y < 0) * g.L - (g.pos[i].y >= g.L) * g.L;
        g.pos[i].z += (g.pos[i].z < 0) * g.L - (g.pos[i].z >= g.L) * g.L;
    }
}

void updater(asakuraoosawa& a, global& g){

    double flucDis = std::sqrt(2.0 * Dt);

    for (int i = 0; i < N; ++i){

        double randX = flucDis * a.gauss(a.gen);
        double randY = flucDis * a.gauss(a.gen);
        double randZ = flucDis * a.gauss(a.gen);

        vec3 brown = {randX, randY, randZ};
        vec3 force = a.force[i];

        g.pos[i] = g.pos[i] + Dt * force + brown;
    }

    periodicAC(g);
}

void AOinit(asakuraoosawa& a, global& g){
    a.force.resize(N);
    a.neighborList.resize(N);

    asakuraoosawaNeighborCells(a, g);
}


#endif