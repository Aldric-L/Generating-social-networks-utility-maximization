//
//  main.cpp
//  Sexual Market Simulation
//
//  Created by SMS&co on 21/03/2023.
//

#include <iostream>
#include <filesystem>
#include <thread>
#include "Constants.hpp"
#include "Individual.hpp"
#include "SexualMarket.hpp"

bool is_number(const std::string& s) { return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end(); }

int main(int argc, const char * argv[]) {
    std::cout << "    _____ __  ________    _____ _                 __      __             \n";
    std::cout << "   / ___//  |/  / ___/   / ___/(_)___ ___  __  __/ /___ _/ /_____  _____ \n";
    std::cout << "   \\__ \\/ /|_/ /\\__ \\    \\__ \\/ / __ `__ \\/ / / / / __ `/ __/ __ \\/ ___/ \n";
    std::cout << "  ___/ / /  / /___/ /   ___/ / / / / / / / /_/ / / /_/ / /_/ /_/ / /     \n";
    std::cout << " /____/_/  /_//____/   /____/_/_/ /_/ /_/\\__,_/_/\\__,_/\\__/\\____/_/      \n";
    std::cout << "\nWelcome in the SMS-Simulator ! \n";
    std::cout << "This build will generate simulations with " << GRAPH_SIZE << " individuals. \n";
    std::string roundsinput = "";
    
    while (!is_number(roundsinput)){
        std::cout << "How many rounds would you like to simulate ? ";
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
    int rounds = std::stoi(roundsinput);
    
    roundsinput = "";
    while (!is_number(roundsinput)){
        std::cout << "\nHow many simulations should be conducted in parallel ? ";
        std::cin >> roundsinput;
    }
    
    int threads_nb = std::stoi(roundsinput);
    
    roundsinput = "";
    while (!is_number(roundsinput) || std::stoi(roundsinput) < 0 || std::stoi(roundsinput) > 100){
        std::cout << "\nWhat is the share of greedy individuals (0-100) ? ";
        std::cin >> roundsinput;
    }
    Individual::GREEDY_SHARE = static_cast<unsigned short int>(std::stoi(roundsinput));
    
    if (threads_nb > 1)
        std::cout.setstate(std::ios_base::failbit);
    
    auto processGame = [](int rds) {
        SexualMarket sm;
        sm.initializeLinks();
        #if GRAPH_SIZE < 50
            std::cout << "\n A look to the initialization adjacency matrix : \n";
            akml::cout_matrix(sm.asAdjacencyMatrix());
        #endif
        unsigned short int inactive_consecutive_rounds_counter(0);
        for (int i(0); i < rds; i++){
                if (inactive_consecutive_rounds_counter == 3){
                    std::cout << "\n\n\n Inactivity detected - Stopping generation at round " << i;
                    break;
                }
                
                if ((sm.processARound()) == GRAPH_SIZE)
                    inactive_consecutive_rounds_counter++;
                else
                    inactive_consecutive_rounds_counter=0;
        }
        #if GRAPH_SIZE < 50
            std::cout << "\n A look to the final adjacency matrix : \n";
            akml::cout_matrix(sm.asAdjacencyMatrix());
        #endif
    };
    if (threads_nb > 1){
        std::vector<std::thread> threads;
        for (int th(0); th < threads_nb; th++){
            std::thread thread(std::ref(processGame), rounds);
            threads.push_back(std::move(thread));
        }
        for (int th(0); th < threads_nb; th++){
            if (threads[th].joinable())
                threads[th].join();
        }
    }else {
        processGame(rounds);
    }
    
    if (threads_nb > 1)
        std::cout.clear();
    
    return 0;
}
