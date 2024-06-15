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
    if (Individual::HETEROGENEOUS_P){
        std::uniform_int_distribution<unsigned short int> inner_distribution(0,10);
        for (std::size_t i(0); i < Individual::P_DIMENSION; i++){
            if (distribution(this->world->getRandomGen())==1)
                Individual::P(i+1,1) = (inner_distribution(this->world->getRandomGen())>3) ? 1 : 0;
            else
                Individual::P(i+1,1) = (inner_distribution(this->world->getRandomGen())>6) ? 1 : 0;
            sum += (distribution(this->world->getRandomGen())==1) ? 1 : 0;
        }
    }else {
        for (std::size_t i(0); i < Individual::P_DIMENSION; i++){
            Individual::P(i+1,1) = (distribution(this->world->getRandomGen())==1) ? 1 : 0;
            sum += (distribution(this->world->getRandomGen())==1) ? 1 : 0;
        }
    }
    if (sum >= (100 - Individual::GREEDY_SHARE))
        Individual::is_greedy = true;
    else
        Individual::is_greedy = false;
    
    
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

/*
 * For utility computing and take action processing, we build 4 matrices:
 * - P_S: the union of the matrices P of individuals in scope
 * - alpha: a column vector of relation weights fot each individual in scope (same order of P_S)
 * - beta: a column vector of pointers to the individuals in the scope
 * - eta: a column vector of pointers to the relations in the scope
 */
Individual::PSAndAlphaTuple Individual::buildPSAndAlpha(const akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>& relations){
    // We count the number of relations that truly exist
    std::size_t effectiveRelationCount = 0;
    for (std::size_t i(0); i < GRAPH_SIZE-1; i++)
        if ((relations)[{i,0}]->weight > 0)
            effectiveRelationCount++;
    
    //std::vector<akml::Matrix<float, P_DIMENSION, 1>> P_S_temp(effectiveRelationCount);
    // Compatibility matrix
    akml::DynamicMatrix<float> P_prod (effectiveRelationCount, 1);
    // Weights matrix
    akml::DynamicMatrix<float> alpha (effectiveRelationCount, 1);
    // Pointers to individual matrix
    akml::DynamicMatrix<Individual*> beta (effectiveRelationCount, 1);
    // Pointers to relation matrix
    akml::DynamicMatrix<SocialMatrix::Link*> eta (effectiveRelationCount, 1);
    
    std::size_t internIncrement (0);
    for (std::size_t i(0); i < GRAPH_SIZE-1; i++){
        if (relations[{i,0}]->weight > 0){
            if (relations[{i,0}]->first == this){
                //P_S_temp[internIncrement] = relations[{i,0}]->second->getP();
                beta(internIncrement+1, 1) = relations[{i,0}]->second;
            }else {
                //P_S_temp[internIncrement] = relations[{i,0}]->first->getP();
                beta(internIncrement+1, 1) = relations[{i,0}]->first;
            }
            eta(internIncrement+1, 1) = relations[{i,0}];
            alpha(internIncrement+1, 1) = relations[{i,0}]->weight;
            P_prod(internIncrement+1, 1) = relations[{i,0}]->compatibility;
            if (relations[{i,0}]->compatibility == -1)
                throw std::runtime_error("Compatibility not properly stored in the relation.");
            internIncrement++;
        }
    }
    //akml::DynamicMatrix<float> P_S (P_S_temp);
    return { /*.P_S = P_S,*/ .P_prod = P_prod, .alpha = alpha, .beta = beta, .eta = eta};
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
        #if GRAPH_SIZE < 100
            std::cout << "\nUtility (" << agentid << " / " << GRAPH_SIZE << ") : No relation U=0 (Gamma=" << Individual::gamma << ")";
        #endif
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
        
    #if GRAPH_SIZE < 100
    std::cout << "\nUtility (" << agentid << " / " << GRAPH_SIZE << ") : LHS=" << LHS << " RHS=" << RHS1+RHS2 << " Total=" << LHS-(RHS1+RHS2) << "(Gamma=" << Individual::gamma << ")";
    #endif
    
    if (rel_td)
        delete relations;
    
    return /*LHS-(RHS1+RHS2)*/ utilityFunc->function(P_prod, alpha);
}

akml::DynamicMatrix<float> Individual::computeUtilityGrad(akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>& relations, Individual::PSAndAlphaTuple& PS_Alpha) {
    //akml::DynamicMatrix<float> P_prod (akml::matrix_product(akml::transpose(PS_Alpha.P_S), this->P));
    //P_prod = P_prod * (1/(float)P_DIMENSION);
    akml::DynamicMatrix<float> grad = utilityFunc->derivative(PS_Alpha.P_prod, PS_Alpha.alpha);
    
    float maxgrad = akml::max(grad);
    if (maxgrad > 1){
        grad = 0.01 * grad;
    }else if (maxgrad > 0.1){
        grad = 0.1 * grad;
    }else if (maxgrad > 0.01){
        
    }else if (maxgrad > 0.001){
        //grad = 10 * grad;
    }
    //grad.transform([](float val) { return val/10; });
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
    akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> relations = this->getRelations();
    akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> rel_temp;
    std::vector<SocialMatrix::Link> scope = this->getScope();
    
    // Because we want to distinguish individuals that are in the scope but with whom there is no link and those  no link
    // with whom there is no link at all, we create a fake little link of 0.00001 for individuals in scope.
    // If greedy we select create a new scope entry
    long long int greedy_target (-1);
    if (Individual::is_greedy){
        std::uniform_int_distribution<unsigned short int> distribution(0,GRAPH_SIZE-2);
        greedy_target = distribution(this->world->getRandomGen());
        if (Individual::GREEDY_FREQ < 100)
            greedy_target = (((greedy_target+1)/(GRAPH_SIZE-1))*100 < Individual::GREEDY_FREQ) ? distribution(this->world->getRandomGen()) : -1;
    }
    
    
    // Pointers to relation matrix
    std::size_t internIncrement (0);
    akml::DynamicMatrix<SocialMatrix::Link*> true_eta (GRAPH_SIZE-1, 1);
    for (std::size_t rel(0); rel < GRAPH_SIZE-1; rel++){
        rel_temp[{rel, 0}] = new SocialMatrix::Link;
        *rel_temp[{rel, 0}] = *relations[{rel, 0}];
        if (relations[{rel, 0}]->weight == 0){
            for (std::size_t i(0); i < scope.size(); i++){
                if ((scope[i].first == this && relations[{rel, 0}]->first == this && scope[i].second == relations[{rel, 0}]->second)
                    || (scope[i].first == this && relations[{rel, 0}]->second == this && scope[i].second == relations[{rel, 0}]->first)
                    || (scope[i].second == this && relations[{rel, 0}]->second == this && scope[i].first == relations[{rel, 0}]->first)
                    || (scope[i].second == this && relations[{rel, 0}]->first == this && scope[i].first == relations[{rel, 0}]->second)){
                    rel_temp[{rel, 0}]->weight = 0.00001;
                    break;
                }
            }
            
            if (greedy_target != -1 && rel == greedy_target)
                rel_temp[{static_cast<std::size_t>(greedy_target), 0}]->weight = 0.00002;
        }
        if (rel_temp[{rel, 0}]->weight > 0){
            true_eta(internIncrement+1, 1) = relations[{rel, 0}];
            internIncrement++;
        }
    }
    
    // No relation found
    if (internIncrement == 0){
        #if GRAPH_SIZE < 100
        std::cout << "\nWant to take action but no individuals are in scope. Aborted.";
        #endif
        for (std::size_t rel(0); rel < GRAPH_SIZE-1; rel++)
            delete rel_temp[{rel, 0}];
        
        SocialMatrix::Link newlinkwanted (nullptr, nullptr, 0);
        return std::make_tuple(nullptr, nullptr, newlinkwanted, false);
    }
    
    true_eta.resize(internIncrement, 1);

    Individual::PSAndAlphaTuple PS_Alpha = this->buildPSAndAlpha(rel_temp);
    
    //We delete temp pointers
    for (std::size_t rel(0); rel < GRAPH_SIZE-1; rel++)
        delete rel_temp[{rel, 0}];
    
    #if GRAPH_SIZE < 100
    std::cout << "Computing grad" << std::endl;
    #endif
    akml::DynamicMatrix<float> grad (Individual::computeUtilityGrad(relations, PS_Alpha));
    //akml::cout_matrix(grad);
    //std::cout << "Now let's reduce the grad to accessible individuals (scope=" << scope.size() << ")" << std::endl;
    for (std::size_t i(0); i < grad.getNRows(); i++){
        // You can't select an individual that you don't see (no freehate)
        // Since the implementation of dynamic matrices this extra-care is not relevant
        if (PS_Alpha.beta(i+1, 1) == nullptr){
            grad(i+1, 1) = 0;
        // You can't reduce a friendship that does not exist
        }else if (grad(i+1, 1) < 0 && PS_Alpha.alpha(i+1, 1) <= 0.0001 ){
            grad(i+1, 1) = 0;
        }
    }
    #if GRAPH_SIZE < 100
        akml::cout_matrix(grad);
    #endif
    if (target == nullptr){
        std::size_t max_i = 0;
        if (MEMORY_SIZE > 0 && memoryBuffer.size() > 0){
            float lmax = -MAXFLOAT;
            for (std::size_t i(0); i < grad.getNRows(); i++){
                if (std::abs(grad[{i,0}]) > lmax) {
                    if (std::find(memoryBuffer.begin(), memoryBuffer.end(), PS_Alpha.beta[{i, 0}]) != memoryBuffer.end()){
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
            #if GRAPH_SIZE < 100
            std::cout << "\nWant to do an unsignifiant action. Aborted.";
            #endif
            return std::make_tuple(nullptr, nullptr, SocialMatrix::Link (nullptr, nullptr, 0), false);
        }
    
        float step = grad(max_i+1, 1);
        // We do not allow to move with a step that is wider that 0.2
        step = std::min(step, (float)MAX_LINK_CHANGE);
        step = std::max(step, (float)-MAX_LINK_CHANGE);
        
        float mov = std::max(0.f, ((PS_Alpha.alpha(max_i+1, 1) < 0.01) ? 0.f : PS_Alpha.alpha(max_i+1, 1)) + step);
        
        
        #if GRAPH_SIZE < 100
        std::cout << "\nWant to take action on the coef " << max_i << " in the direction of " << grad(max_i+1, 1) << " to reach " << mov << " which corresponds currently to a link of " << ((std::get<1>(PS_Alpha)(max_i+1, 1) <= 0.0001) ? 0 : std::get<1>(PS_Alpha)(max_i+1, 1)) << " which links to the individual " << std::get<2>(PS_Alpha)(max_i+1, 1);
        //std::cout << "We verifiy true_eta : " << true_eta(max_i+1, 1)->first << " - " << true_eta(max_i+1, 1)->second << " - " << true_eta(max_i+1, 1)->weight << "\n";
        #endif
        
        if (MEMORY_SIZE > 0){
            memoryBuffer.push_front(PS_Alpha.beta[{max_i, 0}]);
            if (memoryBuffer.size() > MEMORY_SIZE)
                memoryBuffer.pop_back();
        }
        SocialMatrix::Link newlinkwanted (true_eta(max_i+1, 1)->first, true_eta(max_i+1, 1)->second, mov);
        
        return std::make_tuple(true_eta(max_i+1, 1), PS_Alpha.beta(max_i+1, 1), newlinkwanted, (grad(max_i+1, 1) <= 0));
    }else {
        long int target_i = -1;
        for (std::size_t i(0); i < PS_Alpha.beta.getNRows(); i++){
            if (PS_Alpha.beta(i+1, 1) == target){
                target_i = i;
                break;
            }
        }
        
        // target not found
        if (target_i == -1)
            return std::make_tuple(nullptr, nullptr, SocialMatrix::Link (nullptr, nullptr, 0), false);
        
        // too little desire
        if (std::abs(grad(target_i+1, 1)) < MIN_LINK_WEIGHT)
            return std::make_tuple(true_eta(target_i+1, 1), PS_Alpha.beta(target_i+1, 1), SocialMatrix::Link (true_eta(target_i+1, 1)->first, true_eta(target_i+1, 1)->second, 0.f), false);
        
        float mov = std::max(0.f, ((PS_Alpha.alpha(target_i+1, 1) <= 0.0001) ? 0.f : PS_Alpha.alpha(target_i+1, 1)) + grad(target_i+1, 1));
        
        return std::make_tuple(true_eta(target_i+1, 1), PS_Alpha.beta(target_i+1, 1), SocialMatrix::Link (true_eta(target_i+1, 1)->first, true_eta(target_i+1, 1)->second, mov), (grad(target_i+1, 1) <= 0));
        
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
        #if GRAPH_SIZE < 100
        std::cout << "\nAsk to respond to action of " << from << " but we move to " << std::min(new_weight, std::get<2>(prefAction).weight) << std::endl;
        #endif
        if (std::get<2>(prefAction).weight == 0 || std::abs(std::get<2>(prefAction).weight - new_weight) > DEPRECIATION_RATE){
            this->world->editLink(std::get<0>(prefAction), std::min(new_weight, std::get<2>(prefAction).weight), true, false);
            return true;
        }
    }else if (std::get<0>(prefAction) != nullptr ) {
        if (std::get<2>(prefAction).weight <= 0){
            #if GRAPH_SIZE < 100
            std::cout << "\nAsk to respond to action of " << from << " but our desire is negative or null (" << std::get<2>(prefAction).weight << ")" << std::endl;
            #endif
            this->world->editLink(std::get<0>(prefAction), new_weight, false, false);
        }else{
            std::cout << "\nAsk to respond to action of " << from << " but he is not found" << std::endl;
        }
    }
    return false;
}

