#ifndef HARDSPHERE_H
#define HARDSPHERE_H

#include <iostream>
#include "global.h"

struct hardsphere{
  double cutoff = 6.0;
  std::vector<std::vector<int> > neighborList;
  std::vector<std::vector<int> > neighborCells;
};

void hsNeighborCells(global& g, hardsphere& h){

  const int cellsPerDim = g.L / h.cutoff;
  const int numCells = cellsPerDim * cellsPerDim * cellsPerDim;

  h.neighborCells.resize(numCells);
  for (int i = 0; i < cellsPerDim; i++){
    for (int j = 0; j < cellsPerDim; j++){
      for (int k = 0; k < cellsPerDim; k++){
        int cellIndex = i * cellsPerDim * cellsPerDim + j * cellsPerDim + k;

        std::vector<int> cellNeighbors;
        for (int x = -1; x < 2; x++){
          for (int y = -1; y < 2; y++){
            for (int z = -1; z < 2; z++){
              int nx = i + x;
              int ny = j + y;
              int nz = k + z;

              nx += (nx < 0) * cellsPerDim - (nx == cellsPerDim) * cellsPerDim;
              ny += (ny < 0) * cellsPerDim - (ny == cellsPerDim) * cellsPerDim;
              nz += (nz < 0) * cellsPerDim - (nz == cellsPerDim) * cellsPerDim;

              int idx = nx * cellsPerDim * cellsPerDim + ny * cellsPerDim + nz;
              cellNeighbors.push_back(idx);
            }
          }
        }
        h.neighborCells[cellIndex] = cellNeighbors;
      }
    }
  }
}

void hsNeighborLists(global& g, hardsphere& h){
  const int cellsPerDim = g.L / h.cutoff;
  const double cellSize = g.L / cellsPerDim;
  const int numCells = cellsPerDim * cellsPerDim * cellsPerDim;

  h.neighborList.resize(N);

  for (int i = 0; i < N; i++){
    h.neighborList[i].clear();
  }

  std::vector<std::vector<int> > verletList(numCells);

  for (int i = 0; i < N; i++){
    vec3 r = g.pos[i];

    int nx = r.x / cellSize;
    int ny = r.y / cellSize;
    int nz = r.z / cellSize;

    nx = (nx < 0) * cellsPerDim - (nx == cellsPerDim) * cellsPerDim;
    ny = (ny < 0) * cellsPerDim - (ny == cellsPerDim) * cellsPerDim;
    nz = (nz < 0) * cellsPerDim - (nz == cellsPerDim) * cellsPerDim;

    int cellIndex = nx * cellsPerDim * cellsPerDim + ny * cellsPerDim + nz;

    if (cellIndex > numCells) std::cout << "hardsphere " << cellIndex << std::endl;

    verletList.at(cellIndex).push_back(i);
  }

  double squareCut = h.cutoff * h.cutoff;

  for (int i = 0; i < numCells; i++){
    for (int j : verletList[i]){
      vec3 rj = g.pos[j];

      for (int k : h.neighborCells[i]){
        for (int m : verletList[k]){
          if (j >= m) continue;
          vec3 r = rj - g.pos[m];
          double dist = r.squaremag();

          if (dist > squareCut) continue;
          h.neighborList[j].push_back(m);
        }
      }
    }
  }
}

void hsForces(global& g, hardsphere& h){

    for (int i = 0; i < N; i++){
        std::vector<int> neighbors = h.neighborList[i];
        for (int j : neighbors){
            vec3 r = g.pos[i] - g.pos[j];
            double d = r.squaremag();

            if (d > 4.0) continue;

            double rsqrt = std::sqrt(d);

            double outFactor = 1.0 / (2.0 * dt);
            double inFactor = (2 - rsqrt) / rsqrt;

            vec3 force = (inFactor * outFactor) * r;
            g.forces[i] = g.forces[i] + force;
            g.forces[j] = g.forces[j] - force;
        }
    }
}

void HSinit(global& g, hardsphere& h){
    h.neighborList.resize(N);
    hsNeighborCells(g, h);
}

#endif