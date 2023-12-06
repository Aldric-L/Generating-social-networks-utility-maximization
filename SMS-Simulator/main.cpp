//
//  main.cpp
//  Sexual Market Simulation
//
//  Created by SMS&co on 21/03/2023.
//

#include <iostream>
#include <filesystem>
#include "Constants.hpp"
#include "Individual.hpp"
#include "SexualMarket.hpp"

bool is_number(const std::string& s) { return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end(); }

int main(int argc, const char * argv[]) {
    std::cout << "Hello, terrible World!\n";
    
    SexualMarket sm;
    sm.initializeLinks();

    std::cout << "Welcome in the SMS-Simulator \n";
    std::string roundsinput = "";
    
    while (!is_number(roundsinput)){
        std::cout << "How many round would you like to simulate ? ";
        std::cin >> roundsinput;
    }
    std::string log = "";
    while (log != "yes" && log != "y" && log != "no" && log != "n"){
        std::cout << "\nShould I log the results (yes/no) ? ";
        std::cin >> log;
    }
    if (log == "yes" || log == "y"){
        std::cout << "\nLog mode active. Files will be saved in " << std::filesystem::current_path();
        SexualMarket::SHOULD_I_LOG = true;
    }else {
        SexualMarket::SHOULD_I_LOG = false;
    }
        
    std::cout << "\n A look to the initialization adjacency matrix : \n";
    akml::cout_matrix(sm.asAdjacencyMatrix());
    int rounds = std::stoi(roundsinput);
    unsigned short int inactive_consecutive_rounds_counter(0);
    for (int i(0); i < rounds; i++){
        if (inactive_consecutive_rounds_counter == 3){
            std::cout << "\n\n\n Inactivity detected - Stopping generation at round " << i;
            break;
        }
        
        unsigned int in(0);
        if ((in = sm.processARound()) == GRAPH_SIZE)
            inactive_consecutive_rounds_counter++;
        else
            inactive_consecutive_rounds_counter=0;
        //std::cout << "\n\n\n Inactivity = " << in;
    }
    std::cout << "\n A look to the final adjacency matrix : \n";
    akml::cout_matrix(sm.asAdjacencyMatrix());
    
    return 0;
}

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
