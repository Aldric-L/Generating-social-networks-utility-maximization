//
//  Individual.cpp
//  Social Matrix Simulation
//
//  Created by SMS Associates on 21/03/2023.
//

#include "Individual.hpp"
#include "SocialMatrix.hpp"


/*
 * Each individual is constructed with a binary matrix P of dim P_DIMENSION
 * The parameters gamma, is_greedy are random (N(9, 0.35) and U([1, 100-GREEDY_SHARE])
 * Delta is fixed to 2
 */
Individual::Individual(SocialMatrix& world, unsigned long int agentid) : kappa(DEFAULT_KAPPA), delta(DEFAULT_DELTA), world(&world), agentid(agentid), P(Individual::P_DIMENSION, 1)/*, gen(std::random_device{}())*/{
    std::uniform_int_distribution<unsigned short int> distribution(0,1);
    std::size_t sum = 0;
    std::size_t true_sum = 0;
    // We want to prevent the P vector from being the 0 vector
    while (true_sum == 0){
        if (Individual::HETEROGENEOUS_P){
            std::uniform_int_distribution<unsigned short int> inner_distribution(0,10);
            for (std::size_t i(0); i < Individual::P_DIMENSION; i++){
                if (distribution(this->world->getRandomGen())==1)
                    Individual::P(i+1,1) = (inner_distribution(this->world->getRandomGen())>3) ? 1 : 0;
                else
                    Individual::P(i+1,1) = (inner_distribution(this->world->getRandomGen())>6) ? 1 : 0;
                sum += (distribution(this->world->getRandomGen())==1) ? 1 : 0;
                true_sum += Individual::P(i+1,1);
            }
        }else {
            for (std::size_t i(0); i < Individual::P_DIMENSION; i++){
                Individual::P(i+1,1) = (distribution(this->world->getRandomGen())==1) ? 1 : 0;
                sum += (distribution(this->world->getRandomGen())==1) ? 1 : 0;
                true_sum += Individual::P(i+1,1);
            }
        }
        if (sum >= (100 - Individual::GREEDY_SHARE))
            Individual::is_greedy = true;
        else
            Individual::is_greedy = false;
    }

    if (GAMMA_DISP > 0){
        float g = 0;
        std::normal_distribution<float> norm(1,GAMMA_DISP);
        while((g = norm(this->world->getRandomGen())) > 1.8 || g < 0.1 || g==1){ g = norm(this->world->getRandomGen()); }
        Individual::gamma = g+Individual::GAMMA_MEAN-1;
    }else {
        Individual::gamma = Individual::GAMMA_MEAN;
    }
    
    utilityFunc = new ALKYUtility(Individual::gamma, Individual::DEFAULT_KAPPA, Individual::DEFAULT_DELTA);
}


Individual::~Individual(){
    delete utilityFunc;
}

/*
 * Relations and scope are evaluated in the SocialMatrix class, the individual does not compute them directly
 */
akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> Individual::getRelations(){
    return this->world->getIndividualRelations(this);
}

std::vector<SocialMatrix::Link> Individual::getScope() {    
    return this->world->getIndividualScope(this);
}

float Individual::computeUtility(akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>* relations) {
    bool rel_td = false;
    if (relations == nullptr){
        relations = new akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>;
        *relations = this->getRelations();
        rel_td = true;
    }
    std::size_t effectiveRelationCount = 0;
    for (std::size_t i(0); i < GRAPH_SIZE-1; i++){
        if ((*relations)[{i,0}]->weight > 0){
            // We clean relations that are far too weak
            if ((*relations)[{i,0}]->weight < MIN_LINK_WEIGHT)
                (*relations)[{i,0}]->weight = 0;
            else
                effectiveRelationCount++;
            
        }
    }
    
    if (effectiveRelationCount == 0){
        if (rel_td)
            delete relations;
        return 0.f;
    }
    
    akml::DynamicMatrix<float> alpha (effectiveRelationCount, 1);
    akml::DynamicMatrix<float> P_prod(effectiveRelationCount, 1);

    std::size_t internIncrement (0);
    for (std::size_t i(0); i < GRAPH_SIZE-1; i++){
        if ((*relations)[{i,0}]->weight > 0){
            P_prod(internIncrement+1, 1) = (*relations)[{i,0}]->compatibility;
            alpha(internIncrement+1, 1) = (*relations)[{i,0}]->weight;
            internIncrement++;
        }
    }
        
    if (rel_td)
        delete relations;
    
    return utilityFunc->function(P_prod, alpha);
}

akml::DynamicMatrix<float> Individual::computeUtilityGrad(const akml::DynamicMatrix<float>& P_prod, const akml::DynamicMatrix<float>& alpha) {
    akml::DynamicMatrix<float> grad = utilityFunc->derivative(P_prod, alpha);
    
    float maxgrad = akml::max(grad);
    if (maxgrad > 1){
        grad = 0.01 * grad;
    }else if (maxgrad > 0.1){
        grad = 0.1 * grad;
    }
    return grad;
}

/*
 * preprocessTakeAction is a method that chooses what action to do in two cases: when it is turn for the individual to play,
 * or when he is asked to answer to a call (target is here a pointer to the asker).
 *
 * In both cases, this method computes a fake relation matrix with all individuals in scope, in order to give it to
 * the computeUtilityGrad method.
 *
 * If it is his turn: we process an argmax in abs on the grad to moove in the direction of the higher coefficient.
 * If he is asked to answer a call, we only consider the relevant coefficient.
 *
 * If the action is negative, the function returns the desired action but with the 2nd component in nullptr
 * If we need to ask someone, we put the target in the 2nd component.
 */
std::tuple<SocialMatrix::Link*, Individual*, SocialMatrix::Link, bool> Individual::preprocessTakeAction(Individual* target){
    std::vector<SocialMatrix::Link> scope = this->getScope();
    akml::DynamicMatrix<SocialMatrix::Link*> scope_relations (scope.size(), 1);
    akml::DynamicMatrix<SocialMatrix::Link*> true_relations (scope.size(), 1);
    akml::DynamicMatrix<float> alpha (scope.size(), 1);
    akml::DynamicMatrix<float> P_Prod (scope.size(), 1);

    
    for (std::size_t rel(0); rel < scope.size(); rel++){
        scope_relations[{rel, 0}] = &(scope.at(rel));
        true_relations[{rel, 0}] = this->world->findRelation(this, scope.at(rel).second);
        alpha[{rel, 0}] = scope.at(rel).weight;
        P_Prod[{rel, 0}] = scope.at(rel).compatibility;
        
        // Because we want to distinguish individuals that are in the scope but with whom there is no link and those  no link
        // with whom there is no link at all, we create a fake little link of 0.00001 for individuals in scope.
        if (scope.at(rel).weight == 0)
            scope.at(rel).weight = 0.0001;
    }
    
    
    // If greedy we select create a new scope entry
    long long int greedy_target (-1);
    if (Individual::is_greedy){
        std::uniform_int_distribution<unsigned short int> distribution(0,GRAPH_SIZE-2);
        greedy_target = distribution(this->world->getRandomGen());
        if (Individual::GREEDY_FREQ < 100)
            greedy_target = ((greedy_target+1)*100 < Individual::GREEDY_FREQ*(GRAPH_SIZE-1)) ? distribution(this->world->getRandomGen()) : -1;
        
        if (greedy_target != -1){
            bool found = false;
            for (auto& rel : scope){
                if (rel.second->agentid == greedy_target){
                    found = true;
                    break;
                }
            }
            if (!found){
                auto target = this->world->getIndividual(greedy_target);
                scope.emplace_back(this,target, 0.f, this->world->getCompatibilityBtwnIndividuals(this, target));
            }
        }
    }
        
    akml::DynamicMatrix<float> grad (Individual::computeUtilityGrad(P_Prod, alpha));
    
    if (grad.getNRows() == 0){
        SocialMatrix::Link newlinkwanted (nullptr, nullptr, 0);
        return std::make_tuple(nullptr, nullptr, newlinkwanted, false);
    }
    
    for (std::size_t i(0); i < grad.getNRows(); i++){
        // You can't reduce a friendship that does not exist
        if (grad(i+1, 1) < 0 && true_relations(i+1, 1)->weight == 0 ){
            grad(i+1, 1) = 0;
        }
    }

    if (target == nullptr){
        std::size_t max_i = 0;
        if (MEMORY_SIZE > 0 && memoryBuffer.size() > 0){
            float lmax = -MAXFLOAT;
            for (std::size_t i(0); i < grad.getNRows(); i++){
                if (std::abs(grad[{i,0}]) > lmax) {
                    if (std::find(memoryBuffer.begin(), memoryBuffer.end(), scope_relations[{i, 0}]->second) != memoryBuffer.end()){
                        grad[{i,0}] = 0.f;
                    }else{
                        max_i = i;
                        lmax = std::abs(grad[{i,0}]);
                    }
                }
            }
        }else {
            max_i = akml::arg_max(grad, true);
        }
        
        if (std::abs(grad(max_i+1, 1)) < MIN_LINK_WEIGHT){
            return std::make_tuple(nullptr, nullptr, SocialMatrix::Link (nullptr, nullptr, 0), false);
        }
    
        float step = grad(max_i+1, 1);
        // We do not allow to move with a step that is wider that 0.2
        step = std::min(step, (float)MAX_LINK_CHANGE);
        step = std::max(step, (float)-MAX_LINK_CHANGE);
        
        float mov = std::max(0.f, ((alpha(max_i+1, 1) < 0.01) ? 0.f : alpha(max_i+1, 1)) + step);
        
        if (MEMORY_SIZE > 0){
            memoryBuffer.push_front(scope_relations[{max_i, 0}]->second);
            if (memoryBuffer.size() > MEMORY_SIZE)
                memoryBuffer.pop_back();
        }
        SocialMatrix::Link newlinkwanted (true_relations(max_i+1, 1)->first, true_relations(max_i+1, 1)->second, mov);
        
        return std::make_tuple(true_relations(max_i+1, 1), scope_relations[{max_i, 0}]->second, newlinkwanted, (grad(max_i+1, 1) <= 0));
    }else {
        long int target_i = -1;
        for (std::size_t i(0); i < scope_relations.getNRows(); i++){
            if (scope_relations(i+1, 1)->second == target){
                target_i = i;
                break;
            }
        }
        
        // target not found
        if (target_i == -1)
            return std::make_tuple(nullptr, nullptr, SocialMatrix::Link (nullptr, nullptr, 0), false);
        
        // too little desire
        if (std::abs(grad(target_i+1, 1)) < MIN_LINK_WEIGHT)
            return std::make_tuple(true_relations(target_i+1, 1), scope_relations(target_i+1, 1)->second, SocialMatrix::Link (true_relations(target_i+1, 1)->first, true_relations(target_i+1, 1)->second, 0.f), false);
        
        float mov = std::max(0.f, ((alpha(target_i+1, 1) <= 0.0001) ? 0.f : alpha(target_i+1, 1)) + grad(target_i+1, 1));
        
        return std::make_tuple(true_relations(target_i+1, 1), scope_relations(target_i+1, 1)->second, SocialMatrix::Link (true_relations(target_i+1, 1)->first, true_relations(target_i+1, 1)->second, mov), (grad(target_i+1, 1) <= 0));
        
    }
    
}

bool Individual::takeAction(){
    std::tuple<SocialMatrix::Link*, Individual*, SocialMatrix::Link, bool> prefAction = Individual::preprocessTakeAction();
    if (std::get<3>(prefAction) && std::get<0>(prefAction) != nullptr){
        this->world->editLink(std::get<0>(prefAction), std::get<2>(prefAction).weight, true, false);
        return true;
    }else if (std::get<1>(prefAction) != nullptr) {
        return std::get<1>(prefAction)->responseToAction(this, std::get<2>(prefAction).weight);
    }
    return false;
}

bool Individual::responseToAction(Individual* from, float new_weight){
    std::tuple<SocialMatrix::Link*, Individual*, SocialMatrix::Link, bool> prefAction = Individual::preprocessTakeAction(from);
    if (std::get<0>(prefAction) != nullptr && std::get<2>(prefAction).weight > 0){
        if (std::get<2>(prefAction).weight == 0 || std::abs(std::get<2>(prefAction).weight - new_weight) > DEPRECIATION_RATE){
            this->world->editLink(std::get<0>(prefAction), std::min(new_weight, std::get<2>(prefAction).weight), true, false);
            return true;
        }
    }else if (std::get<0>(prefAction) != nullptr ) {
        if (std::get<2>(prefAction).weight <= 0){
            this->world->editLink(std::get<0>(prefAction), new_weight, false, false);
        }else{
            std::cout << "\nAsk to respond to action of " << from << " but he is not found" << std::endl;
        }
    }
    return false;
}

