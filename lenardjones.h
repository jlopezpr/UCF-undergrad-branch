#ifndef LENARDJONES_H
#define LENARDJONES_H
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>

#include "global.h"
#include "init.h"

struct lenardjones {

    double cutoff = 9.55;
    const double epsilon = 100;
    const double sigma = 2.1;

    std::vector<std::vector<int>> neighborCells;
    std::vector<std::vector<int>> neighborList;
    std::vector<double> forceX;
    std::vector<double> forceY;
    std::vector<double> forceZ;

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

    std::vector<std::vector<int>> verletList(numberOfCells);

    for (int i = 0; i < N; ++ i){
        double x = g.posX[i];
        double y = g.posY[i];
        double z = g.posZ[i];

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

            double x = g.posX[i];
            double y = g.posY[i];
            double z = g.posZ[i];

            for (int t : l.neighborCells[s]){
                for (int j : verletList[t]){

                    if (i >= j) continue;
                    double x1 = x - g.posX[j];
                    double y1 = y - g.posY[j];
                    double z1 = z - g.posZ[j];

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
        l.forceX[i] = 0;
        l.forceY[i] = 0;
        l.forceZ[i] = 0;
    }
    
    double epsilon = l.epsilon;
    double sigma = l.sigma;
    double sigma6 = std::pow(sigma, 6);
    epsilon *= sigma6;

    for (int i = 0; i < N; i++){
        std::vector<int> neighbors = l.neighborList[i];
        for (int j : neighbors){

            double dx = g.posX[i] - g.posX[j];
            double dy = g.posY[i] - g.posY[j];
            double dz = g.posZ[i] - g.posZ[j];

            double r = 1.0 / (dx*dx + dy*dy + dz*dz);
            double r3 = std::pow(r, 3);
            double rsqrt = std::sqrt(r);
            
            double outFactor = (24 * epsilon * sigma6) * (r3 * rsqrt);
            double inFactor = (2 * sigma6) * r3 - 1;

            double force = outFactor * inFactor;

            l.forceX[i] -= force * dx;
            l.forceY[i] -= force * dy;
            l.forceZ[i] -= force * dz;

            l.forceX[j] += force * dx;
            l.forceY[j] += force * dy;
            l.forceZ[j] += force * dz;
        }
        
    }
}

void periodicBC(global& g){

    for (int i = 0; i < N; ++i){
        g.posX[i] += (g.posX[i] < 0) * g.L - (g.posX[i] >= g.L) * g.L;
        g.posY[i] += (g.posY[i] < 0) * g.L - (g.posY[i] >= g.L) * g.L;
        g.posZ[i] += (g.posZ[i] < 0) * g.L - (g.posZ[i] >= g.L) * g.L;
    }
}

void updater(lenardjones& l, global& g){

    double flucDis = std::sqrt(2.0 * Dt);

    for (int i = 0; i < N; ++i){

        double randX = flucDis * l.gauss(l.gen);
        double randY = flucDis * l.gauss(l.gen);
        double randZ = flucDis * l.gauss(l.gen);

        // g.posX[i] += l.forceX[i] * Dt + randX;
        // g.posY[i] += l.forceY[i] * Dt + randY;
        // g.posZ[i] += l.forceZ[i] * Dt + randZ;

        g.posX[i] += l.forceX[i] * Dt + randX;
        g.posY[i] += l.forceY[i] * Dt + randY;
        g.posZ[i] += l.forceZ[i] * Dt + randZ;
    }

    periodicBC(g);
}

void LJinit(lenardjones& l, global& g){
    l.forceX.resize(N);
    l.forceY.resize(N);
    l.forceZ.resize(N);
    l.neighborList.resize(N);

    generateNeighborCells(l, g);
}

void simulation(lenardjones& l, global& g){

    initialization(g);
    LJinit(l, g);
    generateNeighborList(l, g);
    std::cout << "Simulation has started" << std::endl;

    for (int i = 0; i < numberOfTimeSteps; ++i){
        generateForces(l, g);
        updater(l, g);
        
        if (i % 100 == 99) generateNeighborList(l, g);
        if ((i + 1) % writtingRate == 0){
            g.timeStep = i;
            updateOVITOFile(g);
        }
    }
}
#endif