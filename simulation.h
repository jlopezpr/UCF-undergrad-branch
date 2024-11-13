#ifndef SIMULATION_H
#define SIMULATION_H

#include "init.h"
#include "hardsphere.h"
#include "asakuraoosawa.h"
#include "yukawa.h"

void simulation(global& g, hardsphere& h, asakuraoosawa& a, yukawa& y){

    initialization(g);

    HSinit(g, h);
    AOinit(g, a);
    YKinit(g, y);

    g.forces.resize(N);

    hsNeighborLists(g, h);
    aoNeighborLists(g, a);
    ykNeighborLists(g, y);

    for (int i = 0; i < timeSteps; ++i){
        aoForces(g, a);
        ykForces(g, y);
        hsForces(g, h);

        updater(g);

        if (i % 50 == 49) {
            ykNeighborLists(g, y);
            hsNeighborLists(g, h);
            aoNeighborLists(g, a);
        }

        if ((i + 1) % writingRate == 0){
            g.timeStep = i + 1;
            updateOVITOFile(g);
        }
    }
}

#endif