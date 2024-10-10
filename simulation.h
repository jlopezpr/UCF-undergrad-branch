#ifndef SIMULATION_H
#define SIMULATION_H

#include "init.h"
#include "hardsphere.h"
#include "yukawa.h"
#include "asakuraoosawa.h"

void simulation(global& g, hardsphere& h, yukawa& y, asakuraoosawa& a){

    initialization(g);
    HJinit(h, g);
    AOinit(a, g);
    YKinit(y, g);

    generateNeighborList(h, g);
    asakuraoosawaNeighborList(a, g); 
    yukawaNeighborList(y, g);

    std::cout << "Simulation has started" << std::endl;

    for (int i = 0; i < numberOfTimeSteps; ++i){
        generateForces(h, g);
        asakuraoosawaForces(a, g);
        yukawaForces(y, g);

        updater(a, g);
        
        if (i % 100 == 99) {
            generateNeighborList(h, g);
            asakuraoosawaNeighborList(a, g);
            yukawaNeighborList(y, g);
        }
        if ((i + 1) % writtingRate == 0){
            g.timeStep = i;
            updateOVITOFile(g);
        }
    }
}
#endif