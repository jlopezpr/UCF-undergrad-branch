#ifndef INIT_H
#define INIT_H
#include "global.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>

int overlappingCheck(global& g){
    int output = 0;
    double invBox = 2.0 / g.L;

    for (int i = 0; i < N; ++i){
        for (int j = 0; j < i; ++j){
            vec3 r = g.pos[i] - g.pos[j];

            r.x -= g.L * int (r.x * invBox);
            r.y -= g.L * int (r.y * invBox);
            r.z -= g.L * int (r.z * invBox);

            double dist = r.squaremag();
            if (dist < 4.0) {
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

    while (overlapping == 1 && counter <= maxOverlap){
        for (int i = 0; i < N; ++i){
          for (int j = 0; j < i; ++j){

            vec3 r = g.pos[i] - g.pos[j];
            double dist = r.squaremag();
            if (dist >= 4.0) continue;

            double rr = std::sqrt(dist);
            rr = (2.0 - rr) / 2.0;

            g.pos[i] = g.pos[i] + rr * r;
            g.pos[j] = g.pos[j] - rr * r;
        }
    }
        overlapping = overlappingCheck(g);
        ++counter;
    }
}

void randomDistribution(global& g){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> gauss(0, 1);

    g.gen = gen;
    g.gauss = gauss;
}

void periodic(global& g){
    for (int i = 0; i < N; ++i){
        g.pos[i].x += (g.pos[i].x < 0) * g.L - (g.pos[i].x >= g.L) * g.L;
        g.pos[i].y += (g.pos[i].y < 0) * g.L - (g.pos[i].y >= g.L) * g.L;
        g.pos[i].z += (g.pos[i].z < 0) * g.L - (g.pos[i].z >= g.L) * g.L;
    }
}

void updater(global& g){

    double flucDis = std::sqrt(2.0 * dt);
    vec3 vec;
    for (int i = 0; i < N; ++i){
        double randX = flucDis * g.gauss(g.gen);
        double randY = flucDis * g.gauss(g.gen);
        double randZ = flucDis * g.gauss(g.gen);

        vec3 brown(randX, randY, randZ);

        g.pos[i] = g.pos[i] + dt * g.forces[i] + brown;

        g.forces[i] = vec;
    }
    periodic(g);
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
        vec3 r = g.pos[i];

        outFile << i  << " " << 1 << " " << r.x << " " << r.y << " " << r.z << "\n";
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
        vec3 r = g.pos[i];
        outFile << i << " " << 1 << " " << r.x << " " << r.y << " " << r.z << "\n";
    }
    
    outFile.close();
}


void initialization(global& g){

    double vol = (4.0 / 3.0) * pi * N / phi;
    g.L = std::cbrt(vol);

    std::random_device rd; //Seed
    std::mt19937 genPos(rd()); //Mersenne twister generator
    std::uniform_real_distribution<> disPos(0.0, 1.0);

    g.pos.resize(N);

    for (int i = 0; i < N; ++i){
      g.pos[i].x = disPos(genPos) * g.L;
      g.pos[i].y = disPos(genPos) * g.L;
      g.pos[i].z = disPos(genPos) * g.L;
    }

    resolveOverlapping(g);
    randomDistribution(g);
    generateOVITOFile(g);
}



#endif