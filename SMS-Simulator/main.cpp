//
//     _____ __  ________    _____ _                 __      __
//   / ___//  |/  / ___/   / ___/(_)___ ___  __  __/ /___ _/ /_____  _____
//   \__ \/ /|_/ /\__ \    \__ \/ / __ `__ \/ / / / / __ `/ __/ __ \/ ___/
//  ___/ / /  / /___/ /   ___/ / / / / / / / /_/ / / /_/ / /_/ /_/ / /
// /____/_/  /_//____/   /____/_/_/ /_/ /_/\__,_/_/\__,_/\__/\____/_/
//
//                                                SMS Associates 2023-2024
//

#include <iostream>
#include <filesystem>
#include <thread>
#include "AKML-lib/AKML.hpp"
#include "Constants.hpp"
#include "UtilityFunction.cpp"
#include "Individual.hpp"
#include "SocialMatrix.hpp"
#include "OptimalMatrix.hpp"
#include "AKML-lib/AgentBasedUtilities/CLInterface.hpp"

int main(int argc, const char * argv[]) {
    std::cout << "    _____ __  ________    _____ _                 __      __             \n";
    std::cout << "   / ___//  |/  / ___/   / ___/(_)___ ___  __  __/ /___ _/ /_____  _____ \n";
    std::cout << "   \\__ \\/ /|_/ /\\__ \\    \\__ \\/ / __ `__ \\/ / / / / __ `/ __/ __ \\/ ___/ \n";
    std::cout << "  ___/ / /  / /___/ /   ___/ / / / / / / / /_/ / / /_/ / /_/ /_/ / /     \n";
    std::cout << " /____/_/  /_//____/   /____/_/_/ /_/ /_/\\__,_/_/\\__,_/\\__/\\____/_/      \n";
    std::cout << "\nWelcome in the SMS-Simulator ! (v. "<< SMS_VERSION << ")\n";
    std::cout << "This build will generate simulations with " << GRAPH_SIZE << " individuals (MT" << ((float)MAX_THREADS_USAGE)/100 << "). \n";
    
    std::size_t rounds = 1000;
    unsigned short int simulationsNb = 1;
    bool shouldITryToFindOptimalGraph = false;
    
    auto CLOptionsTuple = std::make_tuple
        (akml::CLOption<std::size_t> (&rounds, "R", "rounds", "How many rounds?"),
         akml::CLOption<unsigned short int> (&simulationsNb, "S", "simuls", "How many simulations?"),
         akml::CLOption<bool> (&shouldITryToFindOptimalGraph, "O", "optiGraph", "Should I compute optimal graph?"),
         akml::CLOption<unsigned short int> (&Individual::GREEDY_SHARE, "s", "greedyS", "Share of greedy individuals [0,100]"),
         akml::CLOption<unsigned short int> (&Individual::GREEDY_FREQ, "f", "greedyF", "Frequency of the greedy bonus [0,100]", 10),
         akml::CLOption<float> (&Individual::DEFAULT_DELTA, "D", "delta", "Utility parameter", 2),
         akml::CLOption<float> (&Individual::DEFAULT_KAPPA, "K", "kappa", "Utility parameter"),
         akml::CLOption<float> (&Individual::GAMMA_MEAN, "G", "gamma", "Utility parameter", 9),
         akml::CLOption<bool> (&Individual::HETEROGENEOUS_P, "p", "htroP", "Enable/Disable the two groups of P", false),
         akml::CLOption<bool> (&SocialMatrix::COMPUTE_CLEARING, "c", "clearing", "Enable/Disable the clearing and decaying mechanism", true),
         akml::CLOption<bool> (&SocialMatrix::COMPUTE_CLUSTERING, "C", "clustering", "Enable/Disable the computing of clustering coefficients", false),
         akml::CLOption<bool> (&SocialMatrix::SHOULD_I_LOG, "l", "log", "Should we log results?", true),
         akml::CLOption<std::string> (&SocialMatrix::GLOBAL_LOG_PREFIX, "L", "logPrefix", "Select a subfolder for registering logs"));
    
    try {
        akml::CLManager localCLManager(argc, argv, CLOptionsTuple);
    }catch (...) {
        return -1;
    }
    
    if (SocialMatrix::SHOULD_I_LOG)
        std::cout << "Log mode active. Files will be saved in " << std::filesystem::current_path() << "\n\n";
    
    auto processGame = [&CLOptionsTuple, &shouldITryToFindOptimalGraph](std::size_t rds, std::size_t id, std::size_t batchSize) {
        SocialMatrix sm;
        std::string logPath = sm.whereWillYouLog().first;
        std::string logId = sm.whereWillYouLog().second;
        auto start = std::chrono::high_resolution_clock::now();
            sm.initializeLinks();
            if (shouldITryToFindOptimalGraph){
                std::cout << "Simulation " << id << " / " << batchSize << ": Computing the analytical optimal graph...\n";
                OptimalMatrix optiMatComputer;
                auto optiMat = optiMatComputer.compute(sm.asAdjacencyMatrix(), sm.getIndividuals());
                akml::CSV_Saver<akml::FullMatrixSave<akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE>>> optimalAdjacencyMatrixTrackersManager;
                optimalAdjacencyMatrixTrackersManager.addSave(optiMat);
                optimalAdjacencyMatrixTrackersManager.saveToCSV(logPath + "SMS-Save-OptimalGraph-" + logId + ".csv", false);
            }
        
            std::cout << "Simulation " << id << " / " << batchSize << ": Processing simulation...\n";
            unsigned short int inactive_consecutive_rounds_counter(0);
            for (std::size_t i(0); i < rds; i++){
                    if (inactive_consecutive_rounds_counter == 3){
                        std::cout << "Simulation " << id << " / " << batchSize << ": Inactivity detected - Stopping generation at round " << i << "\n";
                        break;
                    }
                    
                    if (sm.processARound(static_cast<std::size_t>(rds)) == GRAPH_SIZE)
                        inactive_consecutive_rounds_counter++;
                    else
                        inactive_consecutive_rounds_counter=0;
            }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        
        if (SocialMatrix::SHOULD_I_LOG){
            std::ofstream cout(std::string(logPath + "SMS-SimulInfos-" + logId +".txt"));
            std::ios_base::sync_with_stdio(false);
            auto *coutbuf = std::cout.rdbuf();
            std::cout.rdbuf(cout.rdbuf());
            akml::CLManager::printOptionsValues(false, CLOptionsTuple);
            std::cout << "--agentsNb=" << GRAPH_SIZE << "\n";
            std::cout << "--PDimension=" << P_DIMENSION << "\n";
            std::cout << "--executionTime=" << duration.count() << "\n";
            std::cout << "--SMSVersion=" << SMS_VERSION << "\n";
            std::cout.rdbuf(coutbuf);
        }
    };
    
    unsigned int maxThreads = (std::thread::hardware_concurrency()*MAX_THREADS_USAGE/100 > 1) ? std::thread::hardware_concurrency()*MAX_THREADS_USAGE/100 : 1;
    
    std::vector<std::thread> workers;
    workers.reserve(maxThreads);
    for (std::size_t s(0); s < simulationsNb;){
        for (unsigned int conc_s(0); s+conc_s < std::min((std::size_t)simulationsNb, s+maxThreads); conc_s++){
            std::cout << "Simulation " << s+conc_s+1 << " / " << simulationsNb << ": Initialization.\n";
            workers.emplace_back(std::ref(processGame), rounds, s+conc_s+1, simulationsNb);
        }
        for (unsigned int conc_s(0); conc_s < workers.size(); conc_s++){
            if (workers[conc_s].joinable())
                workers[conc_s].join();
        }
        s += workers.size();
        workers.clear();
    }
    
    return 0;
}

//                                                                 AA     LL
//                 pzr--}fC/][dn1wYnfJ                            AAAA    LL
//              c{?_x[v){Ur}][ft(vX(|wzrmOcQ                     AA |AA   LL
//            #][_j(L}u-[[{x[{)q|)fztfXxUQuJkUm                 AAAAAAA   LL
//        .])-1-U)t}-v})}}]L/(}hZjjmtXdYCZuCkJCnLC|Xx          |AA   AAA  LLLLLLLL
//      I{(}11_-+/-q{p}f[1Oj{){QnfL|J00kC*#aqh0zQnjvwz/fCj^
//    "tf1f(fxtujXzt(jcf{{vX/t|bOxJ8hmJ0zwLwJCqzCxYLrYXXYcxx?
//   >()vxcJpUaXMUfdtcuj{]XrtrkkQC0#WhwkhokwCCzChbdcQB@BuO|vxn[
//  "f(|xcOhW&OB@B*BMvJ//(cCzXmqzhB%Mo*qrbUZLbma0&paB@B&X*Zw)xruLx/"
//  rv/    a%Mw88(  ./jujn0nvu   mh%%B%B%ad*Qk0Y0#B@@@@@@*puQqdt/ntCujx/)
// .xu  O  oh88c  .`.1ptL|OXY   @@@@c   #mJhd_        BBqhoLYkd#ju/j|rtnfnrjxjx~'
// `YuO  doBa`    nOCOXmQhmw    @@@B"    #ozbJk,         {?{jb*QbO#W*8*#mha0LxfIl`
//  cCZwdk&,      jzzQcak*u.    %%%%B    WoddJmO_            ddddbM*p/U0dbM*p/Ueudk
//   hMW          J[fY:          wBB8    M*p/U0{                 wLwJCqzCxYwLwJCqzfdfd
//                !qv              !qv   `}czjf:                    wLwJCqzCxYwLwJCddddfd
//
