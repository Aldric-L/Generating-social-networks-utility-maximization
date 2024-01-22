//
//  SexualMarket.cpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//
#include "SexualMarket.hpp"
#include "Individual.hpp"

SexualMarket::SexualMarket(){
    #if GRAPH_SIZE >= 100
        std::cout << "\nThread " << std::this_thread::get_id() << " - Status: Creation";
    #endif
    links.reserve(LINKS_NB);
    SexualMarket::EdgeSaveTrackerType::default_parameters_name = {{ "round", "vertex1", "vertex2", "old_weight", "new_weight", "accepted" }};
    SexualMarket::UtilitySaveTrackerType::default_parameters_name = {{ "round", "agentid", "utility" }};
    SexualMarket::VerticesSaveTrackerType::default_parameters_name = {{ "round", "agentid", "gamma", "isgreedy", "meandist", "vardist", "maxdist", "P" }};
    
    SexualMarket::VerticesSaveTrackerType* save;
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        individuals[{indiv,0}] = new Individual(*this, indiv);
    }
    
    int link_i(0);
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            Link l;
            l.first = individuals[{indiv,0}];
            l.second = individuals[{indiv_t,0}];
            l.weight = 0.f;
            SexualMarket::links[link_i] = l;
            link_i++;
        }
    }
}

SexualMarket::~SexualMarket(){
    if (SexualMarket::SHOULD_I_LOG){
        SexualMarket::VerticesSaveTrackerType* save;
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> binaryadjacencymatrix = SexualMarket::asBinaryAdjacencyMatrix();
        akml::Matrix<std::size_t, GRAPH_SIZE, 1> dijkstra_distance_mat;
        for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
            std::string p = "";
            for (std::size_t i(0); i < P_DIMENSION; i++){
                p.push_back( char( individuals[{indiv,0}]->getP()(i+1, 1) + 48) );
            }
            dijkstra_distance_mat = akml::dijkstra_distance_algorithm(binaryadjacencymatrix, indiv);
            //std::cout << "\n Indiv " << indiv << " mean= " << akml::mean(dijkstra_distance_mat, false, ULONG_MAX) << "\n";
            //std::cout << dijkstra_distance_mat;
            save = new SexualMarket::VerticesSaveTrackerType(SexualMarket::currentRound, individuals[{indiv,0}]->agentid, individuals[{indiv,0}]->gamma, individuals[{indiv,0}]->is_greedy, akml::mean(dijkstra_distance_mat, false, ULONG_MAX), akml::stat_var(dijkstra_distance_mat, ULONG_MAX), akml::max(dijkstra_distance_mat), p);
            verticesTrackersManager.addSave(save);
        }
        
        #if COMPUTE_CLUSTERING
        SexualMarket::ClusteringSaveTrackerType* clsave;
        clsave = new SexualMarket::ClusteringSaveTrackerType(std::move(SexualMarket::computeClusteringCoefficients(&binaryadjacencymatrix)));
        SexualMarket::clusteringTrackersManager.addSave(clsave);
        #endif
        long int t = static_cast<long int> (std::clock());
        std::string curt = std::to_string(t);
        SexualMarket::edgeTrackersManager.saveToCSV("SMS-Save-Edges-" + curt + ".csv", false);
        SexualMarket::utilityTrackersManager.saveToCSV("SMS-Save-Utility-" + curt + ".csv", false);
        SexualMarket::verticesTrackersManager.saveToCSV("SMS-Save-Vertices-"  + curt + ".csv", false);
        SexualMarket::clusteringTrackersManager.saveToCSV("SMS-Save-Clustering-" + curt + ".csv", false);
        
    }
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        delete individuals[{indiv,0}];
    }
}

void SexualMarket::initializeLinks(){
    // Any link that is to be initialized should be saved
    SexualMarket::EdgeSaveTrackerType* save;
    int link_i(0);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned short int> distribution(1,GRAPH_SIZE);
    
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            unsigned short int temp = distribution(gen);
            if (temp <= 4){
                SexualMarket::links[link_i].weight = temp * 0.1;
            }
            
            save = new SexualMarket::EdgeSaveTrackerType(SexualMarket::currentRound, SexualMarket::links[link_i].first->agentid, SexualMarket::links[link_i].second->agentid, 0, SexualMarket::links[link_i].weight, true);
            edgeTrackersManager.addSave(save);
            link_i++;
        }
    }
    if (SexualMarket::SHOULD_I_LOG){
        SexualMarket::VerticesSaveTrackerType* vertice_save;
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> binaryadjacencymatrix = SexualMarket::asBinaryAdjacencyMatrix();
        akml::Matrix<std::size_t, GRAPH_SIZE, 1> dijkstra_distance_mat;
        for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
            std::string p = "";
            for (int i(0); i < P_DIMENSION; i++){
                p.push_back( char( individuals[{indiv,0}]->getP()(i+1, 1) + 48) );
            }
            dijkstra_distance_mat = akml::dijkstra_distance_algorithm(binaryadjacencymatrix, indiv);
            vertice_save = new SexualMarket::VerticesSaveTrackerType(SexualMarket::currentRound, individuals[{indiv,0}]->agentid, individuals[{indiv,0}]->gamma, individuals[{indiv,0}]->is_greedy, akml::mean(dijkstra_distance_mat, false, ULONG_MAX), akml::stat_var(dijkstra_distance_mat, ULONG_MAX), akml::max(dijkstra_distance_mat), p);
            verticesTrackersManager.addSave(vertice_save);
            
        }
        #if COMPUTE_CLUSTERING
        SexualMarket::ClusteringSaveTrackerType* clsave;
        clsave = new SexualMarket::ClusteringSaveTrackerType(std::move(SexualMarket::computeClusteringCoefficients(&binaryadjacencymatrix)));
        clusteringTrackersManager.addSave(clsave);
        #endif
    }
    SexualMarket::currentRound = 1;
    #if GRAPH_SIZE >= 100
        std::cout << "\nThread " << std::this_thread::get_id() << " - Status: Initialized";
    #endif
}

akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1> SexualMarket::getIndividualRelations(Individual* indiv) {
    akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1> linksForIndividual;
    unsigned short int incr = 0;
    for (std::size_t i(0); i < LINKS_NB; i++){
        if ((SexualMarket::links[i].first == indiv || SexualMarket::links[i].second == indiv)){
            linksForIndividual[{incr,0}] = &links[i];
            incr++;
        }
    }
    if (incr != GRAPH_SIZE-1)
        throw std::invalid_argument("Something went really wrong");
    return linksForIndividual;
}

akml::Matrix<Individual*, GRAPH_SIZE, 1> SexualMarket::getIndividuals() {
    return SexualMarket::individuals;
}

Individual* SexualMarket::getIndividual(const std::size_t indiv_id) {
    return SexualMarket::individuals[{indiv_id,0}];
}

std::vector<SexualMarket::Link> SexualMarket::getIndividualScope(Individual* indiv) {
    std::vector<SexualMarket::Link> scope;
    akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1> linksofIndividual = SexualMarket::getIndividualRelations(indiv);
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
            
            akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1> linksofTarget = SexualMarket::getIndividualRelations(target0);
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
void SexualMarket::editLink(Individual* indiv1, Individual* indiv2, float newWeight, bool accepted) {
    if (indiv1 == indiv2 || indiv1 == nullptr || indiv2 == nullptr)
        throw std::invalid_argument("Attempting to edit a non-consistent link");
    
    for (std::size_t link_i(0); link_i<LINKS_NB; link_i++){
        if ((SexualMarket::links[link_i].first == indiv1 && SexualMarket::links[link_i].second == indiv2)
            || (SexualMarket::links[link_i].second == indiv1 && SexualMarket::links[link_i].first == indiv2)){
            SexualMarket::editLink(&SexualMarket::links[link_i], newWeight, accepted);
            break;
        }
    }
}

void SexualMarket::editLink(SexualMarket::Link* link, float newWeight, bool accepted) {
    if (link == nullptr)
        throw std::invalid_argument("Attempting to edit a non-consistent link");
    
    SexualMarket::EdgeSaveTrackerType save(SexualMarket::currentRound, link->first->agentid, link->second->agentid, link->weight, newWeight, accepted);
    edgeTrackersManager.addSave(save);
    if (accepted)
        link->weight = newWeight;
}


unsigned int SexualMarket::processARound(std::size_t totalrounds) {
    #if GRAPH_SIZE < 100
        std::cout << "\n\n ---- ROUND " << SexualMarket::currentRound;
    #endif
    #if COMPUTE_CLEARING == 1
    // Test implementation of a usure of weights and a clearing of small weights
    if (SexualMarket::currentRound != 0 && SexualMarket::currentRound != 1 && SexualMarket::currentRound != totalrounds && SexualMarket::currentRound % 10 == 0){
        for (std::size_t link_id(0); link_id < SexualMarket::links.size(); link_id++){
            if (links[link_id].weight < 0.02)
                SexualMarket::editLink(&links[link_id], 0);
            else
                SexualMarket::editLink(&links[link_id], links[link_id].weight-0.01);
        }
    }
    #endif
    
    
    unsigned int inactions(0);
    for (std::size_t i(0); i < GRAPH_SIZE; i++){
        Individual* nodei = SexualMarket::getIndividual(i);
        #if GRAPH_SIZE < 100
            std::cout << "\n --Individual (" << i+1 << " / " << GRAPH_SIZE <<  ") " << nodei << std::endl;
        #endif
        inactions += nodei->takeAction() ? 0 : 1;
    }
    #if GRAPH_SIZE >= 100
        if ((totalrounds != 0 && totalrounds > 100 && SexualMarket::currentRound % (totalrounds/10) == 0) || SexualMarket::currentRound == 1) {
            std::cout << "\n Computing utility: round " << SexualMarket::currentRound << " / " << totalrounds << " thread " << std::this_thread::get_id();
            for (std::size_t i(0); i < GRAPH_SIZE; i++){
                SexualMarket::UtilitySaveTrackerType save (SexualMarket::currentRound, SexualMarket::getIndividual(i)->agentid, SexualMarket::getIndividual(i)->computeUtility(nullptr));
                SexualMarket::utilityTrackersManager.addSave(save);
            }
            #if COMPUTE_CLUSTERING
            SexualMarket::ClusteringSaveTrackerType* clsave;
            auto binaryadjacencymatrix = SexualMarket::asBinaryAdjacencyMatrix();
            clsave = new SexualMarket::ClusteringSaveTrackerType(std::move(SexualMarket::computeClusteringCoefficients(&binaryadjacencymatrix)));
            SexualMarket::clusteringTrackersManager.addSave(clsave);
            #endif
        }else if (totalrounds != 0 && totalrounds > 100 && SexualMarket::currentRound % (totalrounds/100) == 0){
            std::cout << "\nThread " << std::this_thread::get_id() << " - Status: " << SexualMarket::currentRound*100/totalrounds << "% completed";
        }
    #else
        for (std::size_t i(0); i < GRAPH_SIZE; i++){
            SexualMarket::UtilitySaveTrackerType save (SexualMarket::currentRound, SexualMarket::getIndividual(i)->agentid, SexualMarket::getIndividual(i)->computeUtility(nullptr));
            SexualMarket::utilityTrackersManager.addSave(save);
        }
    #endif
    
    #if GRAPH_SIZE >= 100
    if (SexualMarket::currentRound % 100 == 0){
        SexualMarket::edgeTrackersManager.bufferize(false);
        SexualMarket::utilityTrackersManager.bufferize(false);
    }
    #endif
    SexualMarket::currentRound++;
    return inactions;
}

akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> SexualMarket::asAdjacencyMatrix(){
    akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> output;
    for (std::size_t l(0); l<LINKS_NB; l++){
        output(SexualMarket::links[l].first->agentid+1, SexualMarket::links[l].second->agentid+1) = SexualMarket::links[l].weight;
        output(SexualMarket::links[l].second->agentid+1, SexualMarket::links[l].first->agentid+1) = SexualMarket::links[l].weight;
    }
    return output;
}

akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> SexualMarket::asBinaryAdjacencyMatrix(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE>* adjacencymatrix){
    akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> localmatrix;
    if (adjacencymatrix != nullptr){
        for (std::size_t i (0); localmatrix.getStorage() != localmatrix.getStorageEnd(); i++){
            *(localmatrix.getStorage() + i) = (*(adjacencymatrix->getStorage() + i) > 0.005) ? 1 : 0;
        }
    }else {
        for (std::size_t l(0); l<LINKS_NB; l++){
            localmatrix(SexualMarket::links[l].first->agentid+1, SexualMarket::links[l].second->agentid+1) = (SexualMarket::links[l].weight > 0.005) ? 1 : 0;
            localmatrix(SexualMarket::links[l].second->agentid+1, SexualMarket::links[l].first->agentid+1) = (SexualMarket::links[l].weight > 0.005) ? 1 : 0;
        }
    }
    return localmatrix;
}

akml::Matrix<std::size_t, GRAPH_SIZE, GRAPH_SIZE> SexualMarket::computeDegreesOfSeparation(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix){
    bool bam_td = false;
    if (binaryadjacencymatrix == nullptr){
        binaryadjacencymatrix = new akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>;
        *binaryadjacencymatrix = SexualMarket::asBinaryAdjacencyMatrix();
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

akml::Matrix<float, GRAPH_SIZE+1, 1> SexualMarket::computeClusteringCoefficients(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix){
    bool bam_td = false;
    if (binaryadjacencymatrix == nullptr){
        binaryadjacencymatrix = new akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>;
        *binaryadjacencymatrix = SexualMarket::asBinaryAdjacencyMatrix();
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
