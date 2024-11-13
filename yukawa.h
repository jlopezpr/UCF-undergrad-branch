#ifndef YUKAWA_H
#define YUKAWA_H

#include <iostream>
#include "global.h"

struct yukawa{

    double cutoff = 10.0;
    const double A = 8.0;
    const double kappa = 0.7;

    std::vector<std::vector<int> > neighborCells;
    std::vector<std::vector<int> > neighborList;
};

void ykNeighborCells(global& g, yukawa& y){

    const int cellsPerDim = g.L / y.cutoff;
    const int numCells = cellsPerDim * cellsPerDim * cellsPerDim;

    y.neighborCells.resize(numCells);

    for (int i = 0; i < cellsPerDim; i++){
        for (int j = 0; j < cellsPerDim; j++){
            for (int k = 0; k < cellsPerDim; k++){
                int cellIndex = k + j * cellsPerDim + i * cellsPerDim * cellsPerDim;

                std::vector<int> cellNeighbors;

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
                            cellNeighbors.push_back(idx);
                        }
                    }
                }
                y.neighborCells[cellIndex] = cellNeighbors;
            }
        }
    }
}

void ykNeighborLists(global& g, yukawa& y){

    const int cellsPerDim = g.L / y.cutoff;
    const double cellSize = g.L / cellsPerDim;
    const int numCells = cellsPerDim * cellsPerDim * cellsPerDim;

    for (int i = 0; i < N; i++){
        y.neighborList[i].clear();
    }

    std::vector<std::vector<int> > verletList(numCells);

    for (int i = 0; i < N; ++ i){
        vec3 r = g.pos[i];

        int nx = r.x / cellSize;
        int ny = r.y / cellSize;
        int nz = r.z / cellSize;

        nx += (nx < 0) * 1 - (nx == cellsPerDim) * 1;
        ny += (ny < 0) * 1 - (ny == cellsPerDim) * 1;
        nz += (nz < 0) * 1 - (nz == cellsPerDim) * 1;

        int cellIndex = (nx * cellsPerDim * cellsPerDim) + (ny * cellsPerDim) + nz;

        if (cellIndex > numCells) std::cout << "yukawa " << cellIndex << std::endl;

        verletList.at(cellIndex).push_back(i);
    }


    double squareCut = y.cutoff * y.cutoff;

    for (int s = 0; s < numCells; ++s){
        for (int i : verletList[s]){

            vec3 ri = g.pos[i];

            for (int t : y.neighborCells[s]){
                for (int j : verletList[t]){

                    if (i >= j) continue;
                    vec3 r = ri - g.pos[j];
                    double dist = r.squaremag();
                    if (dist > squareCut) continue;
                    y.neighborList[i].push_back(j);
                }
            }
        }
    }
}

void ykForces(global& g, yukawa& y){

    double A = y.A;
    double kappa = y.kappa;

    double squareCut = y.cutoff * y.cutoff;

    for (int i = 0; i < N; i++){
        std::vector<int> neighbors = y.neighborList[i];
        for (int j : neighbors){
            vec3 r = g.pos[i] - g.pos[j];
            double r2 = r.squaremag();

            if (r2 > squareCut) continue;

            double rsqrt = std::sqrt(r2);

            double expo = std::exp(-kappa * rsqrt);
            double outFactor = A * expo / rsqrt;
            double inFactor = (1 / rsqrt) + kappa;

            vec3 force = (outFactor * inFactor) * r;

            g.forces[i] = g.forces[i] - force;
            g.forces[j] = g.forces[j] + force;
        }

    }
}

void YKinit(global& g, yukawa& y){
    y.neighborList.resize(N);
    ykNeighborCells(g, y);
}

#endif