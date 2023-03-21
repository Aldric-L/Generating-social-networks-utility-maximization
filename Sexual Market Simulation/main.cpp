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


    // -------- Faisons des petits tests
    // PEtit world
    SexualMarket world({ 255, 255 }, 100);
    // On prend un pechno au hasard
    Individual* onePechno = world.getAliveIndividuals()[42];
    // On liste ceux qu'il voit
    for (int i(0); i < onePechno->getVisibleIndividuals().size(); i++) {
        std::cout << "Pos: " << onePechno->getVisibleIndividuals()[i]->getCoord()[0] << ":" << onePechno->getVisibleIndividuals()[i]->getCoord()[1] << std::endl;
    }
    return 0;
}
