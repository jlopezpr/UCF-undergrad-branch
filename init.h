#ifndef INIT_H
#define INIT_H
#include "global.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>

int overlappingCheck(global& g){

    int output = 0; 
    double invBox = 2.0 / g.L; // 

    for (int i = 0; i < N; ++i){
        for (int j = 0; j < i; ++j){

            double dx = g.posX[i] - g.posX[j]; // distancia en x entre dos particulas
            double dy = g.posY[i] - g.posY[j]; // distancia en y
            double dz = g.posZ[i] - g.posZ[j]; // distancia en z

            //Minimum image
            dx -= g.L * int (dx * invBox); // como la caja es periodica, al estar en
            dy -= g.L * int (dy * invBox); // lados opuestos, en realidad la distancia
            dz -= g.L * int (dz * invBox); // es menos

            double distance = dx*dx + dy*dy + dz*dz;

        if (distance < 4.0) {
            output++;
            return output;
        }
        }
    }

    return output;
}

void resolveOverlapping(global& g){
    int counter = 0;
    int overlapping = overlappingCheck(g);

    while (overlapping == 1 && counter <= maxOverlapp){
        for (int i = 0; i < N; ++i){
          for (int j = 0; j < i; ++j){

            double dx = g.posX[i] - g.posX[j];
            double dy = g.posY[i] - g.posY[j];
            double dz = g.posZ[i] - g.posZ[j];

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= 4.0) continue;

            double rr = std::sqrt(r2);
            rr = (2.0 - rr) / 2.0;   

            g.posX[i] += rr * dx; 
            g.posY[i] += rr * dy;
            g.posZ[i] += rr * dz; 

            g.posX[j] -= rr * dx; 
            g.posY[j] -= rr * dy;
            g.posZ[j] -= rr * dz;   
        }
    }
        overlapping = overlappingCheck(g);
        ++counter;
    }

    if (counter == maxOverlapp){
        std::cout << "Overlapping is not resolved" << std::endl;
    }
}

void generateOVITOFile(global& g) {
    
    std::ofstream outFile("test.dump");
    outFile << "ITEM: TIMESTEP\n";
    outFile << g.timeStep << "\n";  // Timestep number
    outFile << "ITEM: NUMBER OF ATOMS\n";
    outFile << N << "\n";
    outFile << "ITEM: BOX BOUNDS pp pp pp\n";  // 'pp' means periodic boundaries in each direction
    outFile << "0 " << g.L << "\n";
    outFile << "0 " << g.L << "\n";
    outFile << "0 " << g.L << "\n";

    outFile << "ITEM: ATOMS id type x y z\n";
    for (int i = 0; i < N; ++i) {
        double x = g.posX[i];
        double y = g.posY[i];
        double z = g.posZ[i];

        outFile << i  << " " << 1 << " " << x << " " << y << " " << z << "\n";
        // pos = pos + 0.001 * ori;
        // outFile << i  << " " << 2 << " " << pos.x << " " << pos.y << " " << pos.z << "\n";
    }
    outFile.close();
}

void updateOVITOFile(global& g) {
    
    std::ofstream outFile("test.dump", std::ios::app);
    outFile << "ITEM: TIMESTEP\n";
    outFile << g.timeStep << "\n";  // Timestep number
    outFile << "ITEM: NUMBER OF ATOMS\n";
    outFile << N << "\n";
    outFile << "ITEM: BOX BOUNDS pp pp pp\n";  // 'pp' means periodic boundaries in each direction
    outFile << "0 " << g.L << "\n";
    outFile << "0 " << g.L << "\n";
    outFile << "0 " << g.L << "\n";

    outFile << "ITEM: ATOMS id type x y z\n";
    for (int i = 0; i < N; ++i) {
        double x = g.posX[i];
        double y = g.posY[i];
        double z = g.posZ[i];

        outFile << i  << " " << 1 << " " << x << " " << y << " " << z << "\n";
        // pos = pos + 0.001 * ori;
        // outFile << i  << " " << 2 << " " << pos.x << " " << pos.y << " " << pos.z << "\n";
    }
    outFile.close();
}

void initialization(global& g){

    double vol = (4.0 / 3.0) * Pi * N / phi;
    g.L = std::cbrt(vol);

    std::random_device rd; //Seed
    std::mt19937 genPos(rd()); //Mersenne twister generator
    std::uniform_real_distribution<> disPos(0.0, 1.0);

    g.posX.resize(N);
    g.posY.resize(N);
    g.posZ.resize(N);

    for (int i = 0; i < N; ++i){
      g.posX[i] = disPos(genPos) * g.L;
      g.posY[i] = disPos(genPos) * g.L;
      g.posZ[i] = disPos(genPos) * g.L;
    }
    resolveOverlapping(g);
    generateOVITOFile(g);
}   


#endif