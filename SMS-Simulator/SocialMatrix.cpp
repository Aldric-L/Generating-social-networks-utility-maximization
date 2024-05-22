//
//  SocialMatrix.cpp
//  Social Matrix Simulation
//
//  Created by SMS Associates on 21/03/2023.
//
#include "SocialMatrix.hpp"
#include "Individual.hpp"

SocialMatrix::SocialMatrix(){
    links.reserve(LINKS_NB);
    edgeTrackersManager.setParameterNames({{ "round", "vertex1", "vertex2", "old_weight", "new_weight", "accepted" }});
    utilityTrackersManager.setParameterNames({{ "round", "agentid", "utility" }});
    verticesTrackersManager.setParameterNames({{ "round", "agentid", "gamma", "isgreedy", "meandist", "vardist", "maxdist", "P" }});

    SocialMatrix::VerticesSaveTrackerType* save;
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        individuals[{indiv,0}] = new Individual(*this, indiv);
    }
    
    std::size_t link_i(0);
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            Link l;
            l.first = individuals[{indiv,0}];
            l.second = individuals[{indiv_t,0}];
            l.weight = 0.f;
            SocialMatrix::links[link_i] = l;
            link_i++;
        }
    }
    
    long int t = static_cast<long int> (std::clock());
    logID = std::to_string(t);
    
    logPath = GLOBAL_LOG_PREFIX;
    #if MODE_FOLDER_LOG
    logPath += "sim_" + logID + "/";
    #endif
    std::filesystem::create_directories(logPath);
}

SocialMatrix::~SocialMatrix(){
    if (SocialMatrix::SHOULD_I_LOG){
        SocialMatrix::VerticesSaveTrackerType* save;
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> binaryadjacencymatrix = SocialMatrix::asBinaryAdjacencyMatrix();
        akml::Matrix<std::size_t, GRAPH_SIZE, 1> dijkstra_distance_mat;
        for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
            std::string p = "";
            for (std::size_t i(0); i < P_DIMENSION; i++){
                p.push_back( char( individuals[{indiv,0}]->getP()(i+1, 1) + 48) );
            }
            dijkstra_distance_mat = akml::dijkstra_distance_algorithm(binaryadjacencymatrix, indiv);
            verticesTrackersManager.addSave(SocialMatrix::currentRound, individuals[{indiv,0}]->agentid, individuals[{indiv,0}]->gamma, individuals[{indiv,0}]->is_greedy, akml::mean(dijkstra_distance_mat, false, ULONG_MAX), akml::stat_var(dijkstra_distance_mat, ULONG_MAX), akml::max(dijkstra_distance_mat), p);
        }
        
        if (SocialMatrix::COMPUTE_CLUSTERING)
            SocialMatrix::clusteringTrackersManager.addSave(SocialMatrix::computeClusteringCoefficients(&binaryadjacencymatrix));
        SocialMatrix::finalAdjacencyMatrixTrackersManager.addSave(SocialMatrix::asAdjacencyMatrix());
        
        SocialMatrix::edgeTrackersManager.saveToCSV(logPath + "SMS-Save-Edges-" + logID + ".csv", false);
        SocialMatrix::utilityTrackersManager.saveToCSV(logPath + "SMS-Save-Utility-" + logID + ".csv", false);
        SocialMatrix::verticesTrackersManager.saveToCSV(logPath + "SMS-Save-Vertices-"  + logID + ".csv", false);
        if (SocialMatrix::COMPUTE_CLUSTERING)
            SocialMatrix::clusteringTrackersManager.saveToCSV(logPath + "SMS-Save-Clustering-" + logID + ".csv", false);
        SocialMatrix::finalAdjacencyMatrixTrackersManager.saveToCSV(logPath + "SMS-Save-AdjacencyMatrix-" + logID + ".csv", false);
        
    }
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        delete individuals[{indiv,0}];
    }
}

void SocialMatrix::initializeLinks(){
    // Any link that is to be initialized should be saved
    SocialMatrix::EdgeSaveTrackerType* save;
    int link_i(0);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned short int> distribution(1,GRAPH_SIZE);
    
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            unsigned short int temp = distribution(gen);
            if (temp <= 4){
                SocialMatrix::links[link_i].weight = temp * 0.1;
            }
        
            edgeTrackersManager.addSave(SocialMatrix::currentRound, SocialMatrix::links[link_i].first->agentid, SocialMatrix::links[link_i].second->agentid, 0, SocialMatrix::links[link_i].weight, true);
            link_i++;
        }
    }
    if (SocialMatrix::SHOULD_I_LOG){
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> binaryadjacencymatrix = SocialMatrix::asBinaryAdjacencyMatrix();
        akml::Matrix<std::size_t, GRAPH_SIZE, 1> dijkstra_distance_mat;
        for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
            std::string p = "";
            for (int i(0); i < P_DIMENSION; i++){
                p.push_back( char( individuals[{indiv,0}]->getP()(i+1, 1) + 48) );
            }
            dijkstra_distance_mat = akml::dijkstra_distance_algorithm(binaryadjacencymatrix, indiv);
            verticesTrackersManager.addSave(SocialMatrix::currentRound, individuals[{indiv,0}]->agentid, individuals[{indiv,0}]->gamma, individuals[{indiv,0}]->is_greedy, akml::mean(dijkstra_distance_mat, false, ULONG_MAX), akml::stat_var(dijkstra_distance_mat, ULONG_MAX), akml::max(dijkstra_distance_mat), p);
            
        }
        if (SocialMatrix::COMPUTE_CLUSTERING)
            clusteringTrackersManager.addSave(SocialMatrix::computeClusteringCoefficients(&binaryadjacencymatrix));
        
        SocialMatrix::finalAdjacencyMatrixTrackersManager.addSave(SocialMatrix::asAdjacencyMatrix());
    }
    SocialMatrix::currentRound = 1;
    /*#if GRAPH_SIZE >= 100
        std::cout << "Thread " << std::this_thread::get_id() << " - Status: Initialized\n";
    #endif*/
}

akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> SocialMatrix::getIndividualRelations(Individual* indiv) {
    akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> linksForIndividual;
    unsigned short int incr = 0;
    for (std::size_t i(0); i < LINKS_NB; i++){
        if ((SocialMatrix::links[i].first == indiv || SocialMatrix::links[i].second == indiv)){
            linksForIndividual[{incr,0}] = &links[i];
            incr++;
        }
    }
    if (incr != GRAPH_SIZE-1)
        throw std::invalid_argument("Something went really wrong");
    return linksForIndividual;
}

akml::Matrix<Individual*, GRAPH_SIZE, 1> SocialMatrix::getIndividuals() {
    return SocialMatrix::individuals;
}

Individual* SocialMatrix::getIndividual(const std::size_t indiv_id) {
    return SocialMatrix::individuals[{indiv_id,0}];
}

std::vector<SocialMatrix::Link> SocialMatrix::getIndividualScope(Individual* indiv) {
    std::vector<SocialMatrix::Link> scope;
    akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> linksofIndividual = SocialMatrix::getIndividualRelations(indiv);
    for (std::size_t level0(0); level0 < GRAPH_SIZE-1; level0++){
        if (linksofIndividual[{level0, 0}]->weight > 0){
            Individual* target0 = nullptr;
            Individual* target1 = nullptr;
            
            if (linksofIndividual[{level0, 0}]->first == indiv)
                target0 = linksofIndividual[{level0, 0}]->second;
            else
                target0 = linksofIndividual[{level0, 0}]->first;
            
            bool redundant = false;
            for (int icheck(0); icheck < scope.size(); icheck++){
                if ((scope[icheck].first == linksofIndividual[{level0, 0}]->first && scope[icheck].second == linksofIndividual[{level0, 0}]->second)
                    || (scope[icheck].first == linksofIndividual[{level0, 0}]->second && scope[icheck].second == linksofIndividual[{level0, 0}]->first))
                    redundant = true;
            }
            if (!redundant)
                scope.push_back(*linksofIndividual[{level0, 0}]);
            
            akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> linksofTarget = SocialMatrix::getIndividualRelations(target0);
            for (std::size_t level1(0); level1 < GRAPH_SIZE-1; level1++){
                if (linksofTarget[{level1, 0}]->weight > 0){
                    if (linksofTarget[{level1, 0}]->first == target0)
                        target1 = linksofTarget[{level1, 0}]->second;
                    else
                        target1 = linksofTarget[{level1, 0}]->first;
                    
                    if (target0 != target1 && target1 != indiv && target0 != indiv){
                        Link l;
                        l.first = indiv;
                        l.second = target1;
                        l.weight = linksofTarget[{level1, 0}]->weight * linksofIndividual[{level0, 0}]->weight;
                        bool redundant = false;
                        for (int icheck(0); icheck < scope.size(); icheck++){
                            if ((scope[icheck].first == l.first && scope[icheck].second == l.second)
                                || (scope[icheck].first == l.second && scope[icheck].second == l.first))
                                redundant = true;
                        }
                        if (!redundant)
                            scope.push_back(l);
                    }
                }
            }
            
        }
    }
    return scope;
}

// Actually, I think this method should not exists : be optimized, use pointers.
void SocialMatrix::editLink(Individual* indiv1, Individual* indiv2, float newWeight, bool accepted) {
    if (indiv1 == indiv2 || indiv1 == nullptr || indiv2 == nullptr)
        throw std::invalid_argument("Attempting to edit a non-consistent link");
    
    for (std::size_t link_i(0); link_i<LINKS_NB; link_i++){
        if ((SocialMatrix::links[link_i].first == indiv1 && SocialMatrix::links[link_i].second == indiv2)
            || (SocialMatrix::links[link_i].second == indiv1 && SocialMatrix::links[link_i].first == indiv2)){
            SocialMatrix::editLink(&SocialMatrix::links[link_i], newWeight, accepted);
            break;
        }
    }
}

void SocialMatrix::editLink(SocialMatrix::Link* link, float newWeight, bool accepted) {
    if (link == nullptr)
        throw std::invalid_argument("Attempting to edit a non-consistent link");
    
    edgeTrackersManager.addSave(SocialMatrix::currentRound, link->first->agentid, link->second->agentid, link->weight, newWeight, accepted);
    if (accepted)
        link->weight = newWeight;
}


unsigned int SocialMatrix::processARound(std::size_t totalrounds) {
    #if GRAPH_SIZE < 100
        std::cout << "\n\n ---- ROUND " << SocialMatrix::currentRound;
    #endif
    if (SocialMatrix::COMPUTE_CLEARING){
        // Test implementation of a usure of weights and a clearing of small weights
        if (SocialMatrix::currentRound != 0 && SocialMatrix::currentRound != 1 && SocialMatrix::currentRound != totalrounds && SocialMatrix::currentRound % 10 == 0){
            for (std::size_t link_id(0); link_id < SocialMatrix::links.size(); link_id++){
                if (links[link_id].weight < 0.02)
                    SocialMatrix::editLink(&links[link_id], 0);
                else
                    SocialMatrix::editLink(&links[link_id], links[link_id].weight-0.01);
            }
        }
    }
        
    std::uniform_int_distribution<std::size_t> distribution(0,GRAPH_SIZE-1);
    unsigned int inactions(0);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::size_t start_indiv = distribution(gen);
    
    for (std::size_t i(start_indiv); i < GRAPH_SIZE; i++){
        Individual* nodei = SocialMatrix::getIndividual(i);
        #if GRAPH_SIZE < 100
            std::cout << "\n --Individual (" << i+1 << " / " << GRAPH_SIZE <<  ") " << nodei << std::endl;
        #endif
        inactions += nodei->takeAction() ? 0 : 1;
    }
    for (std::size_t i(0); i < start_indiv; i++){
        Individual* nodei = SocialMatrix::getIndividual(i);
        #if GRAPH_SIZE < 100
            std::cout << "\n --Individual (" << i+1 << " / " << GRAPH_SIZE <<  ") " << nodei << std::endl;
        #endif
        inactions += nodei->takeAction() ? 0 : 1;
    }
    #if GRAPH_SIZE >= 100
        if ((MODE_ECO_LOG && ((totalrounds != 0 && totalrounds > 100 && SocialMatrix::currentRound % totalrounds/10 == 0) || SocialMatrix::currentRound == 1))
            || (!MODE_ECO_LOG && SocialMatrix::currentRound % 5 == 0 ) ) {
            //std::cout << "\n Computing utility: round " << SocialMatrix::currentRound << " / " << totalrounds << " thread " << std::this_thread::get_id();
            for (std::size_t i(0); i < GRAPH_SIZE; i++){
                SocialMatrix::utilityTrackersManager.addSave(SocialMatrix::currentRound, SocialMatrix::getIndividual(i)->agentid, SocialMatrix::getIndividual(i)->computeUtility(nullptr));
            }
            
            if (!MODE_ECO_LOG && SocialMatrix::currentRound % 250 == 0){
                SocialMatrix::finalAdjacencyMatrixTrackersManager.addSave(SocialMatrix::asAdjacencyMatrix());
                if (SocialMatrix::COMPUTE_CLUSTERING){
                    auto binaryadjacencymatrix = SocialMatrix::asBinaryAdjacencyMatrix();
                    SocialMatrix::clusteringTrackersManager.addSave(SocialMatrix::computeClusteringCoefficients(&binaryadjacencymatrix));
                }
            }
        }/*else if (totalrounds != 0 && totalrounds > 100 && (SocialMatrix::currentRound+1) % (totalrounds/100) == 0 && ((SocialMatrix::currentRound+1)*100/totalrounds)%10==0 ){
            std::cout << "Thread " << std::this_thread::get_id() << " - Status: " << SocialMatrix::currentRound*100/totalrounds << "% completed\n";
        }*/
    #else
        for (std::size_t i(0); i < GRAPH_SIZE; i++){
            SocialMatrix::utilityTrackersManager.addSave(SocialMatrix::currentRound, SocialMatrix::getIndividual(i)->agentid, SocialMatrix::getIndividual(i)->computeUtility(nullptr));
        }
    #endif
    
    #if GRAPH_SIZE >= 100
    if (SocialMatrix::currentRound % 100 == 0){
        SocialMatrix::edgeTrackersManager.bufferize(false);
        SocialMatrix::utilityTrackersManager.bufferize(false);
    }
    #endif
    SocialMatrix::currentRound++;
    return inactions;
}

akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> SocialMatrix::asAdjacencyMatrix(){
    akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> output;
    for (std::size_t l(0); l<LINKS_NB; l++){
        output(SocialMatrix::links[l].first->agentid+1, SocialMatrix::links[l].second->agentid+1) = SocialMatrix::links[l].weight;
        output(SocialMatrix::links[l].second->agentid+1, SocialMatrix::links[l].first->agentid+1) = SocialMatrix::links[l].weight;
    }
    return output;
}

akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> SocialMatrix::asBinaryAdjacencyMatrix(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE>* adjacencymatrix){
    akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> localmatrix;
    if (adjacencymatrix != nullptr){
        for (std::size_t i (0); localmatrix.getStorage() != localmatrix.getStorageEnd(); i++){
            *(localmatrix.getStorage() + i) = (*(adjacencymatrix->getStorage() + i) > 0.005) ? 1 : 0;
        }
    }else {
        for (std::size_t l(0); l<LINKS_NB; l++){
            localmatrix(SocialMatrix::links[l].first->agentid+1, SocialMatrix::links[l].second->agentid+1) = (SocialMatrix::links[l].weight > 0.005) ? 1 : 0;
            localmatrix(SocialMatrix::links[l].second->agentid+1, SocialMatrix::links[l].first->agentid+1) = (SocialMatrix::links[l].weight > 0.005) ? 1 : 0;
        }
    }
    return localmatrix;
}

akml::Matrix<std::size_t, GRAPH_SIZE, GRAPH_SIZE> SocialMatrix::computeDegreesOfSeparation(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix){
    bool bam_td = false;
    if (binaryadjacencymatrix == nullptr){
        binaryadjacencymatrix = new akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>;
        *binaryadjacencymatrix = SocialMatrix::asBinaryAdjacencyMatrix();
        bam_td = true;
    }
    std::array<akml::Matrix<std::size_t, GRAPH_SIZE, 1>, GRAPH_SIZE> preresult;
    for (std::size_t i(0); i < GRAPH_SIZE; i++){
        preresult[i] = akml::dijkstra_distance_algorithm(*binaryadjacencymatrix, i);
    }
    
    if (bam_td)
        delete binaryadjacencymatrix;
    
    akml::Matrix<std::size_t, GRAPH_SIZE, GRAPH_SIZE> result(preresult);
    return result;
}

akml::Matrix<float, GRAPH_SIZE+1, 1> SocialMatrix::computeClusteringCoefficients(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix){
    bool bam_td = false;
    if (binaryadjacencymatrix == nullptr){
        binaryadjacencymatrix = new akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>;
        *binaryadjacencymatrix = SocialMatrix::asBinaryAdjacencyMatrix();
        bam_td = true;
    }
    akml::Matrix<float, GRAPH_SIZE+1, 1> result;
    std::size_t pairs(0), possiblepairs(0), global_closed_triples, global_dim_scopes;
    for (std::size_t node_id(0); node_id < GRAPH_SIZE; node_id++){
        std::size_t dimscope(0);
        for (std::size_t node_target_id(0); node_target_id < GRAPH_SIZE; node_target_id++){
            if (node_id != node_target_id && (*binaryadjacencymatrix)[{node_id, node_target_id}] == 1)
                dimscope++;
            if (node_id != node_target_id && (*binaryadjacencymatrix)[{node_id, node_target_id}] == 1){
                for (std::size_t node_target2_id(0); node_target2_id < GRAPH_SIZE; node_target2_id++){
                    if (node_target2_id != node_id && node_target2_id != node_target_id
                        && (*binaryadjacencymatrix)[{node_id, node_target2_id}] == 1){
                        if ((*binaryadjacencymatrix)[{node_target_id, node_target2_id}] == 1){
                            global_closed_triples++;
                            pairs++;
                        }
                        possiblepairs++;
                    }
                }
            }
        }
        
        if (pairs > 0 && possiblepairs >0 && GRAPH_SIZE > 1){
            result[{node_id, 0}] = ((float)pairs) / ((float)possiblepairs);
        }else {
            result[{node_id, 0}] = 0;
        }
        pairs = 0;
        possiblepairs = 0;
        global_dim_scopes += dimscope;
    }
    
    if (global_closed_triples > 0 && global_dim_scopes >0){
        result[{GRAPH_SIZE, 0}] = (float)global_closed_triples / ((float)global_dim_scopes * ((float)global_dim_scopes-1));
    }else {
        result[{GRAPH_SIZE, 0}] = 0;
    }

    if (bam_td)
        delete binaryadjacencymatrix;
    
    return result;
}


std::pair<std::string, std::string> SocialMatrix::whereWillYouLog() {
    return std::make_pair(logPath, logID);
}
