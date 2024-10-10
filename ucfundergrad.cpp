#include <iostream>
#include "simulation.h"
#include "init.h"


int main(){

    global g;
    hardsphere h;
    asakuraoosawa a;
    yukawa y;
    
    simulation(g, h, y, a);

    return 0;
}