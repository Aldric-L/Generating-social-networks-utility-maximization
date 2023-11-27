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

bool is_number(const std::string& s) { return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end(); }

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
    
    /*for (int i(0); i < 10; i++){
        
        Individual* node8 = sm.getIndividuals()[7];
        Individual* nodei = sm.getIndividuals()[i];
        std::cout << "\n --Individual " << nodei << std::endl;
        
        std::vector<SexualMarket::Link> scopei = nodei->getScope();
        nodei->takeAction();
    }
    */
    
    std::cout << "Welcome in the SMS-Simulator \n";
    std::string roundsinput = "";
    
    while (!is_number(roundsinput)){
        std::cout << "How many round would you like to simulate ?";
        std::cin >> roundsinput;
    }
    std::cout << "\n A look to the initialization adjacency matrix : \n";
    akml::cout_matrix(sm.asAdjacencyMatrix());
    int rounds = std::stoi(roundsinput);
    for (int i(0); i < rounds; i++){
        sm.processARound();
    }
    std::cout << "\n A look to the final adjacency matrix : \n";
    akml::cout_matrix(sm.asAdjacencyMatrix());
    
    return 0;
}
