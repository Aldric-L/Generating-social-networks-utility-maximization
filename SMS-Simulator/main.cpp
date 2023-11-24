//
//  main.cpp
//  Sexual Market Simulation
//
//  Created by SMS&co on 21/03/2023.
//

#include <iostream>
#include "Constants.hpp"
#include "Individual.hpp"
#include "SexualMarket.hpp"

int main(int argc, const char * argv[]) {
    std::cout << "Hello, terrible World!\n";
    
    SexualMarket sm;
    sm.initializeLinks();
    //Individual* node8 = sm.getIndividuals()[7];
    //std::array<SexualMarket::Link*, GRAPH_SIZE-1> test = sm.getIndividualRelations(node8);
    //std::vector<SexualMarket::Link> testScope = node8->getScope();
    //std::array<SexualMarket::Link*, GRAPH_SIZE-1> testScope0 = sm.getIndividuals()[0]->getRelations();
    //std::vector<SexualMarket::Link> scope8 = node8->getScope();
    //node8->computeUtility();
    
    for (int i(0); i < 10; i++){
        
        Individual* node8 = sm.getIndividuals()[7];
        Individual* nodei = sm.getIndividuals()[i];
        std::cout << "\n --Individual " << nodei << std::endl;
        
        std::vector<SexualMarket::Link> scopei = nodei->getScope();
        nodei->takeAction();
    }
    
    
    return 0;
}
