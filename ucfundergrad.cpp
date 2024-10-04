#include <iostream>
#include "simulation.h"
#include "init.h"
#include "global.h"

int main(){

    global g;
    lenardjones l;
    asakuraoosawa a;
    yukawa y;
    
    simulation(g, l, y, a);

    return 0;
}