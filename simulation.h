#ifndef SIMULATION_H
#define SIMULATION_H

#include "global.h"
#include "init.h"
#include "lenardjones.h"
#include "yukawa.h"
#include "asakuraoosawa.h"

void simulation(global& g, lenardjones& l, yukawa& y, asakuraoosawa& a){

    initialization(g);
    LJinit(l, g);
    AOinit(a, g);
    YKinit(y, g);

    generateNeighborList(l, g);
    asakuraoosawaNeighborList(a, g); 
    yukawaNeighborList(y, g);

    std::cout << "Simulation has started" << std::endl;

    for (int i = 0; i < numberOfTimeSteps; ++i){
        generateForces(l, g);
        asakuraoosawaForces(a, g);
        yukawaForces(y, g);

        updater(a, g);
        
        if (i % 100 == 99) {
            generateNeighborList(l, g);
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