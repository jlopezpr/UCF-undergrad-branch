#ifndef ASAKURAOOSAWA_H
#define ASAKURAOOSAWA_H

#include <iostream>
#include "global.h"

struct asakuraoosawa {
    double cutoff = 4.0;
    const double epsilon = 100;
    const double sigma = 2.1;

    std::vector<std::vector<int> > neighborCells;
    std::vector<std::vector<int> > neighborList;
};

void aoNeighborCells(global& g, asakuraoosawa& a){

    const int cellsPerDim = g.L / a.cutoff;

    const int numCells = cellsPerDim * cellsPerDim * cellsPerDim;

    a.neighborCells.resize(numCells);

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
                a.neighborCells[cellIndex] = cellNeighbors;
            }
        }
    }
}

void aoNeighborLists(global& g, asakuraoosawa& a){

    const int cellsPerDim = g.L / a.cutoff;
    const double cellSize = g.L / cellsPerDim;
    const int numCells = cellsPerDim * cellsPerDim * cellsPerDim;

    for (int i = 0; i < N; i++){
        a.neighborList[i].clear();
    }

    std::vector<std::vector<int> > verletList(numCells);

    for (int i = 0; i < N; ++ i){
        vec3 r = g.pos[i];

        int nx = r.x / cellSize;
        int ny = r.y / cellSize;
        int nz = r.z / cellSize;

        nx = (nx < 0) * cellsPerDim - (nx == cellsPerDim) * cellsPerDim;
        ny = (ny < 0) * cellsPerDim - (ny == cellsPerDim) * cellsPerDim;
        nz = (nz < 0) * cellsPerDim - (nz == cellsPerDim) * cellsPerDim;

        int cellIndex = nx * cellsPerDim * cellsPerDim + ny * cellsPerDim + nz;

        if (cellIndex > numCells) std::cout << "asakura " << cellIndex << std::endl;

        verletList.at(cellIndex).push_back(i);
    }

    double squareCut = a.cutoff * a.cutoff;

    for (int s = 0; s < numCells; ++s){
        for (int i : verletList[s]){
            vec3 ri = g.pos[i];

            for (int t : a.neighborCells[s]){
                for (int j : verletList[t]){
                    if (i >= j) continue;

                    vec3 r = ri - g.pos[j];
                    double dist = r.squaremag();

                    if (dist > squareCut) continue;

                    a.neighborList[i].push_back(j);
                }
            }
        }
    }
}

void aoForces(global& g, asakuraoosawa& a){

    double delta = 0.1 + 1;
    double R = 2 * delta;
    double Ua = -3.0; // proportional to polymer concentration
    double R2 = R * R;
    double R3 = R2 * R;
    double denominator = 2 * R3 - 6 * R2 + 8;
    vec3 vec;

    for (int i = 0; i < N; i++){
        std::vector<int> neighbors = a.neighborList[i];
        for (int j : neighbors){
            vec3 r = g.pos[i] - g.pos[j];
            double d = r.squaremag();
            if (d < 4.0 || d > R2) continue;

            double rsqrt = std::sqrt(d);

            double numerator = 3 * (R2 - d);

            double inFactor = numerator / denominator;
            double outFactor = Ua / rsqrt;

            vec3 force = (inFactor * outFactor) * r;
            g.forces[i] = g.forces[i] + force;
            g.forces[j] = g.forces[j] - force;
        }

    }
}

void AOinit(global& g, asakuraoosawa& a){
    a.neighborList.resize(N);
    aoNeighborCells(g, a);
}

#endif