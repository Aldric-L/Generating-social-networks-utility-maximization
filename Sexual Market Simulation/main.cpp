//
//  main.cpp
//  Sexual Market Simulation
//
//  Created by SMS&co on 21/03/2023.
//

#include <iostream>
#include "Individual.hpp"
#include "SexualMarket.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";


    // Faisons des petits tests
    SexualMarket world({ 255, 255 }, 100);
    Individual* onePechno = world.getAliveIndividuals()[42];
    std::vector<Individual *> visibles = world.getVisibleIndividuals(100, onePechno->getCoord());
    for (int i(0); i < visibles.size(); i++) {
        std::cout << "Pos: " << visibles[i]->getCoord()[0] << ":" << visibles[i]->getCoord()[1] << std::endl;
    }
    return 0;
}
