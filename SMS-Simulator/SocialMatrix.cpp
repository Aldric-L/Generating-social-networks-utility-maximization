//
//  SocialMatrix.cpp
//  Social Matrix Simulation
//
//  Created by SMS Associates on 21/03/2023.
//
#include "SocialMatrix.hpp"
#include "Individual.hpp"

SocialMatrix::SocialMatrix() : gen(std::random_device{}()){
    edgeTrackersManager.setParameterNames({{ "round", "vertex1", "vertex2", "old_weight", "new_weight", "accepted", "forced" }});
    utilityTrackersManager.setParameterNames({{ "round", "agentid", "utility" }});
    verticesTrackersManager.setParameterNames({{ "round", "agentid", "gamma", "isgreedy", "meandist", "vardist", "maxdist", "P" }});

    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++)
        individuals[{indiv,0}] = new Individual(*this, indiv);
    
    std::size_t link_i(0);
    links.reserve((GRAPH_SIZE*GRAPH_SIZE-GRAPH_SIZE)/2);
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            SocialMatrix::links.emplace_back(
                                                individuals[{indiv,0}],
                                                individuals[{indiv_t,0}],
                                                0.f,
                                                akml::inner_product(individuals[{indiv,0}]->getP(), individuals[{indiv_t,0}]->getP()) / Individual::P_DIMENSION    );
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

SocialMatrix::SocialMatrix(const akml::DynamicMatrix<float>& compatibilityMatrix) : SocialMatrix::SocialMatrix() {
    if (compatibilityMatrix.getNRows() != compatibilityMatrix.getNColumns() || compatibilityMatrix.getNColumns() != GRAPH_SIZE)
        throw std::invalid_argument("Compatibility matrix is not properly sized for the graph");
    
    std::size_t link_i(0);
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            SocialMatrix::links[link_i].compatibility = compatibilityMatrix[{indiv, indiv_t}];
            link_i++;
        }
    }
};

void SocialMatrix::forceEditCompatibilityMatrix(const akml::DynamicMatrix<float>& compatibilityMatrix) {
    if (compatibilityMatrix.getNRows() != compatibilityMatrix.getNColumns() || compatibilityMatrix.getNColumns() != GRAPH_SIZE)
        throw std::invalid_argument("Compatibility matrix is not properly sized for the graph");
    
    std::size_t link_i(0);
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            SocialMatrix::links[link_i].compatibility = compatibilityMatrix[{indiv, indiv_t}];
            link_i++;
        }
    }
};

SocialMatrix::~SocialMatrix(){
    if (SocialMatrix::SHOULD_I_LOG){
        SocialMatrix::VerticesSaveTrackerType* save;
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> binaryadjacencymatrix = SocialMatrix::asBinaryAdjacencyMatrix();
        akml::Matrix<std::size_t, GRAPH_SIZE, 1> dijkstra_distance_mat;
        for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
            std::string p = "";
            for (std::size_t i(0); i < Individual::P_DIMENSION; i++){
                p.push_back( char( individuals[{indiv,0}]->getP()(i+1, 1) + 48) );
            }
            if (SocialMatrix::currentRound == 1){
                verticesTrackersManager.addSave(SocialMatrix::currentRound, individuals[{indiv,0}]->agentid, individuals[{indiv,0}]->gamma, individuals[{indiv,0}]->is_greedy, -1,-1, -1, p);
            }else {
                dijkstra_distance_mat = akml::dijkstra_distance_algorithm(binaryadjacencymatrix, indiv);
                verticesTrackersManager.addSave(SocialMatrix::currentRound, individuals[{indiv,0}]->agentid, individuals[{indiv,0}]->gamma, individuals[{indiv,0}]->is_greedy, akml::mean(dijkstra_distance_mat, false, ULONG_MAX), akml::stat_var(dijkstra_distance_mat, ULONG_MAX), akml::max(dijkstra_distance_mat), p);
            }
            
        }
        SocialMatrix::verticesTrackersManager.saveToCSV(logPath + "SMS-Save-Vertices-"  + logID + ".csv", false);
        if (SocialMatrix::currentRound != 1){
            if (SocialMatrix::COMPUTE_CLUSTERING)
                SocialMatrix::clusteringTrackersManager.addSave(SocialMatrix::computeClusteringCoefficients(&binaryadjacencymatrix));
            
            SocialMatrix::edgeTrackersManager.saveToCSV(logPath + "SMS-Save-Edges-" + logID + ".csv", false);
            SocialMatrix::utilityTrackersManager.saveToCSV(logPath + "SMS-Save-Utility-" + logID + ".csv", false);
            
            if (SocialMatrix::COMPUTE_CLUSTERING)
                SocialMatrix::clusteringTrackersManager.saveToCSV(logPath + "SMS-Save-Clustering-" + logID + ".csv", false);
        }
        SocialMatrix::finalAdjacencyMatrixTrackersManager.addSave(SocialMatrix::asAdjacencyMatrix());
        SocialMatrix::finalAdjacencyMatrixTrackersManager.saveToCSV(logPath + "SMS-Save-AdjacencyMatrix-" + logID + ".csv", false);
        
        
    }
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++)
        delete individuals[{indiv,0}];
}

SocialMatrix::Link::Link(Individual* first, Individual* second, float weight) : first(first), second(second), weight(weight) {
    if (first != nullptr && second != nullptr)
        compatibility = akml::inner_product(first->getP(), second->getP()) / Individual::P_DIMENSION;
};

void SocialMatrix::initializeLinks(){
    std::size_t link_i(0);
    
    std::uniform_int_distribution<unsigned short int> distribution(1,GRAPH_SIZE);
    
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            unsigned short int temp = distribution(gen);
            if (temp <= 4){
                SocialMatrix::links[link_i].weight = temp * 0.1;
            }
        
            // Any link that is to be initialized should be saved
            edgeTrackersManager.addSave(SocialMatrix::currentRound, SocialMatrix::links[link_i].first->agentid, SocialMatrix::links[link_i].second->agentid, 0, SocialMatrix::links[link_i].weight, true, true);
            link_i++;
        }
    }
    if (SocialMatrix::SHOULD_I_LOG){
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> binaryadjacencymatrix(SocialMatrix::asBinaryAdjacencyMatrix());
        akml::Matrix<std::size_t, GRAPH_SIZE, 1> dijkstra_distance_mat;
        for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
            std::string p = "";
            for (int i(0); i < Individual::P_DIMENSION; i++){
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
}

void SocialMatrix::initializeLinks(const akml::DynamicMatrix<float>& adjacencyMatrix) {
    if (adjacencyMatrix.getNRows() != adjacencyMatrix.getNColumns() || adjacencyMatrix.getNColumns() != GRAPH_SIZE)
        throw std::invalid_argument("Adjacency matrix is not properly sized for the graph");
    
    std::size_t link_i(0);
    
    for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (std::size_t indiv_t(0); indiv_t<indiv; indiv_t++){
            SocialMatrix::links[link_i].weight = adjacencyMatrix[{indiv, indiv_t}];
        
            // Any link that is to be initialized should be saved
            edgeTrackersManager.addSave(SocialMatrix::currentRound, SocialMatrix::links[link_i].first->agentid, SocialMatrix::links[link_i].second->agentid, 0, SocialMatrix::links[link_i].weight, true, true);
            link_i++;
        }
    }
    if (SocialMatrix::SHOULD_I_LOG){
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> binaryadjacencymatrix(SocialMatrix::asBinaryAdjacencyMatrix());
        akml::Matrix<std::size_t, GRAPH_SIZE, 1> dijkstra_distance_mat;
        for (std::size_t indiv(0); indiv<GRAPH_SIZE; indiv++){
            std::string p = "";
            for (int i(0); i < Individual::P_DIMENSION; i++){
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
}

akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> SocialMatrix::getIndividualRelations(Individual* indiv) {
    akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> linksForIndividual;
    unsigned short int incr = 0;
    for (std::size_t i(0); i < links.size(); i++){
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
    akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> linksofIndividual(SocialMatrix::getIndividualRelations(indiv));
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
                        //l.compatibility = akml::inner_product(indiv->getP(), target1->getP()) / Individual::P_DIMENSION;
                        l.compatibility = getCompatibilityBtwnIndividuals(indiv, target1);
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
void SocialMatrix::editLink(const Individual* indiv1, const Individual* indiv2, const float newWeight, bool accepted, bool forced) {
    if (indiv1 == indiv2 || indiv1 == nullptr || indiv2 == nullptr)
        throw std::invalid_argument("Attempting to edit a non-consistent link");
    
    for (std::size_t link_i(0); link_i<links.size(); link_i++){
        if ((SocialMatrix::links[link_i].first == indiv1 && SocialMatrix::links[link_i].second == indiv2)
            || (SocialMatrix::links[link_i].second == indiv1 && SocialMatrix::links[link_i].first == indiv2)){
            SocialMatrix::editLink(&SocialMatrix::links[link_i], newWeight, accepted, forced);
            break;
        }
    }
}

void SocialMatrix::editLink(SocialMatrix::Link* link, const float newWeight, bool accepted, bool forced) {
    if (link == nullptr)
        throw std::invalid_argument("Attempting to edit a non-consistent link");

    // We truly forbid insignificant moves that are inflating the log and keeping the simulation running forever
    if (!forced && accepted && newWeight > 0 && newWeight < MIN_LINK_WEIGHT){
        if (newWeight < link->weight){
            edgeTrackersManager.addSave(SocialMatrix::currentRound, link->first->agentid, link->second->agentid, link->weight, 0, accepted, forced);
            link->weight = newWeight;
        }
        else{
            std::cerr << "Forbidden link: " << link->weight << " - " << newWeight <<" \n";
        }
        return;
    }
        
    edgeTrackersManager.addSave(SocialMatrix::currentRound, link->first->agentid, link->second->agentid, link->weight, newWeight, accepted, forced);
    if (accepted || forced)
        link->weight = newWeight;
}


unsigned int SocialMatrix::processARound(std::size_t totalrounds) {
    if (SocialMatrix::COMPUTE_CLEARING){
        // Test implementation of a usure of weights and a clearing of small weights
        if (SocialMatrix::currentRound != 0 && SocialMatrix::currentRound != 1 && SocialMatrix::currentRound != totalrounds && SocialMatrix::currentRound % 10 == 0){
            for (std::size_t link_id(0); link_id < SocialMatrix::links.size(); link_id++){
                if (links[link_id].weight > 0 && links[link_id].weight < (MIN_LINK_WEIGHT + DEPRECIATION_RATE))
                    SocialMatrix::editLink(&links[link_id], 0, true, true);
                else if (links[link_id].weight > 0)
                    SocialMatrix::editLink(&links[link_id], links[link_id].weight-DEPRECIATION_RATE, true, true);
            }
        }
    }
        
    std::uniform_int_distribution<std::size_t> distribution(0,GRAPH_SIZE-1);
    unsigned int inactions(0);
    
    std::size_t start_indiv = distribution(gen);
    
    for (std::size_t i(start_indiv); i < GRAPH_SIZE; i++){
        Individual* nodei = SocialMatrix::getIndividual(i);
        inactions += nodei->takeAction() ? 0 : 1;
    }
    for (std::size_t i(0); i < start_indiv; i++){
        Individual* nodei = SocialMatrix::getIndividual(i);
        inactions += nodei->takeAction() ? 0 : 1;
    }
        if (
            (MODE_ECO_LOG && ((totalrounds != 0 && totalrounds > 100 && SocialMatrix::currentRound % totalrounds/10 == 0) || SocialMatrix::currentRound == 1))
            || (!MODE_ECO_LOG && ((totalrounds < 50 && (totalrounds%5==0 || totalrounds < 5)) || (totalrounds > 50 && SocialMatrix::currentRound % UTILITY_COMPUTATION_INTERVAL == 0) ) )
            ) {
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
        }
    
    SocialMatrix::currentRound++;
    return inactions;
}

akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> SocialMatrix::asAdjacencyMatrix() const{
    akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> output;
    for (std::size_t l(0); l<links.size(); l++){
        output(SocialMatrix::links[l].first->agentid+1, SocialMatrix::links[l].second->agentid+1) = SocialMatrix::links[l].weight;
        output(SocialMatrix::links[l].second->agentid+1, SocialMatrix::links[l].first->agentid+1) = SocialMatrix::links[l].weight;
    }
    return output;
}

akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> SocialMatrix::asBinaryAdjacencyMatrix(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE>* adjacencymatrix) const{
    akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> localmatrix;
    if (adjacencymatrix != nullptr){
        for (std::size_t i (0); localmatrix.getStorage() != localmatrix.getStorageEnd(); i++){
            *(localmatrix.getStorage() + i) = (*(adjacencymatrix->getStorage() + i) > 0.005) ? 1 : 0;
        }
    }else {
        for (std::size_t l(0); l<links.size(); l++){
            localmatrix(SocialMatrix::links[l].first->agentid+1, SocialMatrix::links[l].second->agentid+1) = (SocialMatrix::links[l].weight > 0.005) ? 1 : 0;
            localmatrix(SocialMatrix::links[l].second->agentid+1, SocialMatrix::links[l].first->agentid+1) = (SocialMatrix::links[l].weight > 0.005) ? 1 : 0;
        }
    }
    return localmatrix;
}

akml::Matrix<std::size_t, GRAPH_SIZE, GRAPH_SIZE> SocialMatrix::computeDegreesOfSeparation(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix) const{
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

akml::Matrix<float, GRAPH_SIZE+1, 1> SocialMatrix::computeClusteringCoefficients(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix) const{
    bool bam_td = false;
    if (binaryadjacencymatrix == nullptr){
        binaryadjacencymatrix = new akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>;
        *binaryadjacencymatrix = SocialMatrix::asBinaryAdjacencyMatrix();
        bam_td = true;
    }
    akml::Matrix<float, GRAPH_SIZE+1, 1> result;
    std::size_t pairs(0), possiblepairs(0), global_closed_triples(0), global_dim_scopes(0);
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


std::pair<std::string, std::string> SocialMatrix::whereWillYouLog() const {
    return std::make_pair(logPath, logID);
}

float SocialMatrix::getCompatibilityBtwnIndividuals(const Individual* indiv1, const Individual* indiv2) const {
    const std::size_t id1 = indiv1->agentid;
    const std::size_t id2 = indiv2->agentid;
    
    const std::size_t j = std::min(id1, id2);
    const std::size_t i = std::max(id1, id2);
    
    std::size_t cols = 0;
    for (std::size_t k(1); k < i; k++)
        cols += k;

    
    std::size_t pos = cols+j;
    assert((links[pos].first == indiv1 && links[pos].second == indiv2) || (links[pos].first == indiv2 && links[pos].second == indiv1));
    return links[pos].compatibility;
}


bool SocialMatrix::checkLoveTriangleCondition() const {
    for (std::size_t i(0); i < GRAPH_SIZE; i++){
        std::size_t checks = 0;
        for (std::size_t m(0); m < GRAPH_SIZE; m++){
            for (std::size_t k(0); k < GRAPH_SIZE; k++){
                if(m!=k && k!= i && m != i && getCompatibilityBtwnIndividuals(individuals[{m, 0}], individuals[{k, 0}]) < getCompatibilityBtwnIndividuals(individuals[{i, 0}], individuals[{k, 0}]) + getCompatibilityBtwnIndividuals(individuals[{m, 0}], individuals[{i, 0}])){
                    ++checks;
                    goto endloop;
                }
            }
        }
        endloop:
        if (checks < 1)
            return false;
    }
    return true;
}
