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
 * The parameters gamma, is_greedy are random (N(GAMMA_MEAN, GAMMA_DISP) and U([1, 100-GREEDY_SHARE])
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
        if (SocialMatrix::sig((*relations)[{i,0}]->weight) > 0){
            // We clean relations that are far too weak
            //if ((*relations)[{i,0}]->weight < MIN_LINK_WEIGHT)
            //    (*relations)[{i,0}]->weight = SocialMatrix::sigInverse(0);
            //else
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
        if (SocialMatrix::sig((*relations)[{i,0}]->weight) > 0){
            P_prod(internIncrement+1, 1) = (*relations)[{i,0}]->compatibility;
            alpha(internIncrement+1, 1) = SocialMatrix::sig((*relations)[{i,0}]->weight);
            internIncrement++;
        }
    }
        
    if (rel_td)
        delete relations;
    
    return utilityFunc->function(P_prod, alpha);
}

akml::DynamicMatrix<float> Individual::computeUtilityGrad(const akml::DynamicMatrix<float>& P_prod, const akml::DynamicMatrix<float>& alpha) {
    return akml::hadamard_product(utilityFunc->derivative(P_prod, akml::transform(alpha, SocialMatrix::sig)), akml::transform(alpha, (std::function<float(float)>)[&](float x){return (std::isinf(x)) ? SocialMatrix::sigDerivative(SocialMatrix::sigInverse(MIN_LINK_WEIGHT2/2)) : SocialMatrix::sigDerivative(x); }));
}

/*
 * preprocessTakeAction is a method that chooses what action to do in two cases: when it is turn for the individual to play,
 * or when he is asked to answer to a call (target is here a pointer to the asker).
 *
 * If it is his turn: we process a softmax in abs on the grad to move in the direction of the higher coefficient.
 * If he is asked to answer a call, we only consider the relevant coefficient.
 *
 * If the action is negative, the function returns the desired action but with the 2nd component in nullptr
 * If we need to ask someone, we put the target in the 2nd component.
 */
std::tuple<SocialMatrix::Link*, Individual*, SocialMatrix::Link, bool> Individual::preprocessTakeAction(Individual* target){
    std::vector<SocialMatrix::Link> scope = this->getScope();
    // target_i is the position in scope of the target. If target is not nullptr, it is the position of target. Else, it is the position of the selected individual.
    std::size_t target_i;
    
    // If greedy we select create a new scope entry
    std::size_t greedy_target (std::numeric_limits<std::size_t>::infinity());
    if (Individual::is_greedy && SocialMatrix::SCOPE_DEPTH < GRAPH_SIZE){
        std::uniform_int_distribution<unsigned short int> distribution(0,GRAPH_SIZE-1);
        greedy_target = distribution(this->world->getRandomGen());
        if (Individual::GREEDY_FREQ < 100)
            greedy_target = ((greedy_target)*100 < Individual::GREEDY_FREQ*(GRAPH_SIZE-1)) ? distribution(this->world->getRandomGen()) : std::numeric_limits<std::size_t>::infinity();
        
        if (greedy_target != std::numeric_limits<std::size_t>::infinity()){
            while (greedy_target == this->agentid)
                greedy_target = distribution(this->world->getRandomGen());
            
            bool found = false;
            for (auto& rel : scope){
                if (rel.second->agentid == greedy_target){
                    found = true;
                    break;
                }
            }
            if (!found){
                auto greedy_target_indiv = this->world->getIndividual(greedy_target);
                scope.emplace_back(this,greedy_target_indiv, SocialMatrix::sigInverse(0.f), this->world->getCompatibilityBtwnIndividuals(this, greedy_target_indiv));
            }
        }
    }
    
    // If we have been asked to answer to someone that we don't have in scope (it must be a greedy individual) we add it in the scope
    if (target != nullptr){
        auto pos_i = std::find_if(scope.begin(), scope.end(), [target](const SocialMatrix::Link& element){ return element.second == target;});
        if (pos_i == scope.end()){
            scope.emplace_back(this,target, SocialMatrix::sigInverse(0.f), this->world->getCompatibilityBtwnIndividuals(this, target));
            target_i = scope.size()-1;
        }else {
            target_i = pos_i - scope.begin();
        }
    }
    
    akml::DynamicMatrix<SocialMatrix::Link*> scope_relations (scope.size(), 1);
    akml::DynamicMatrix<SocialMatrix::Link*> true_relations (scope.size(), 1);
    akml::DynamicMatrix<float> alpha (scope.size(), 1);
    akml::DynamicMatrix<float> P_Prod (scope.size(), 1);
    for (std::size_t rel(0); rel < scope.size(); rel++){
        // We build scope relations because it has links with the target in .second contrary to true relations...
        scope_relations[{rel, 0}] = &(scope.at(rel));
        true_relations[{rel, 0}] = this->world->findRelation(this, scope.at(rel).second);
        alpha[{rel, 0}] = scope.at(rel).weight;
        P_Prod[{rel, 0}] = scope.at(rel).compatibility;
    }
        
    akml::DynamicMatrix<float> grad (Individual::computeUtilityGrad(P_Prod, alpha));
    
    if (grad.getNRows() == 0)
        return std::make_tuple(nullptr, nullptr, SocialMatrix::Link (nullptr, nullptr, 0), false);
    
    for (std::size_t i(0); i < grad.getNRows(); i++){
        // You can't reduce a friendship that does not exist
        if (grad[{i, 0}] < 0 && true_relations[{i, 0}]->weight <= SocialMatrix::sigInverse(0.f) )
            grad[{i, 0}] = 0;
        
        if (std::isnan(grad[{i, 0}]))
            grad[{i, 0}] = 0.f;
    }

    if (target == nullptr){
        akml::DynamicMatrix<float> softMaxScores = grad;
        akml::DynamicMatrix<float> softMaxScoresAdjusted(grad.getNRows(), 1);
        akml::DynamicMatrix<float> softMaxScoresIds(grad.getNRows(), 1);
        std::size_t nulls = 0;
        for (std::size_t grad_i(0); grad_i < softMaxScores.getNRows(); grad_i++){
            if (std::abs(softMaxScores[{grad_i, 0}]) < MIN_GRADIENT_MOVE ||
                (std::abs(softMaxScores[{grad_i, 0}]) >= MIN_GRADIENT_MOVE && MEMORY_SIZE > 0 && memoryBuffer.size() > 0 &&
                 std::find_if(memoryBuffer.begin(), memoryBuffer.end(), [&scope_relations, grad_i](const std::pair<Individual*, bool>& element){ return element.first == scope_relations[{grad_i, 0}]->second;}) != memoryBuffer.end())){
                softMaxScores[{grad_i, 0}] = 0;
                nulls++;
            }else {
                softMaxScoresAdjusted[{grad_i-nulls, 0}] = std::abs(softMaxScores[{grad_i, 0}]);
                softMaxScoresIds[{grad_i-nulls, 0}] = grad_i;
            }
        }
        softMaxScoresAdjusted.resize(softMaxScores.getNRows()-nulls, 1);
        softMaxScoresIds.resize(softMaxScores.getNRows()-nulls, 1);

        if (nulls == softMaxScores.getNRows()){
            if (MEMORY_SIZE > 0 && memoryBuffer.size() > 0){
                for (std::size_t mem_i(memoryBuffer.size()); mem_i > 0; mem_i--){
                    if (memoryBuffer.at(mem_i-1).second == true){
                        memoryBuffer.erase(memoryBuffer.begin() + mem_i-1);
                        break;
                    }
                }
            }
            return std::make_tuple(nullptr, nullptr, SocialMatrix::Link (nullptr, nullptr, 0), false);
        }
        target_i = 0;
        if (softMaxScoresAdjusted.getNRows() > 1)
            target_i = sampleFromSoftmax(softmax(softMaxScoresAdjusted), this->world->getRandomGen());
            
        target_i = softMaxScoresIds[{target_i, 0}];
    }
    
    //std::cout << this << " wants to ask " << true_relations[{target_i, 0}]->second << " because its current link " << SocialMatrix::sig(alpha[{target_i, 0}]) << " makes him a grad of " << grad[{target_i, 0}] << "\n";
    float optiWeight = -MAXFLOAT;
    if (SocialMatrix::sig(alpha[{target_i, 0}]) <= MIN_LINK_WEIGHT2 && grad[{target_i, 0}] <= 0)
        goto returnsection;
    
    if (SocialMatrix::sig(alpha[{target_i, 0}]) == 0){
        if (SocialMatrix::sig(grad[{target_i, 0}]) < MAX_GRADIENT_MOVE){
            optiWeight = grad[{target_i, 0}];
        }else {
            optiWeight = SocialMatrix::sigInverse(MAX_GRADIENT_MOVE);
        }
        goto returnsection;
    }
    
    optiWeight = alpha[{target_i, 0}] + 0.1 * grad[{target_i, 0}];
    if (SocialMatrix::sig(optiWeight) <= MIN_LINK_WEIGHT2){
        optiWeight = SocialMatrix::sigInverse(0.f);
        goto returnsection;
    }
    
    if (std::abs(SocialMatrix::sig(optiWeight) - SocialMatrix::sig(alpha[{target_i, 0}])) <= MIN_GRADIENT_MOVE){
        optiWeight = alpha[{target_i, 0}];
        goto returnsection;
    }
    
    if (std::abs(SocialMatrix::sig(optiWeight) - SocialMatrix::sig(alpha[{target_i, 0}])) > MAX_GRADIENT_MOVE){
        optiWeight = SocialMatrix::sigInverse(SocialMatrix::sig(alpha[{target_i, 0}]) + ((grad[{target_i, 0}] > 0) ? 1.f : -1.f) * MAX_GRADIENT_MOVE);
    }
        
    returnsection:
    if (optiWeight == alpha[{target_i, 0}] && target == nullptr)
        return std::make_tuple(nullptr, nullptr, SocialMatrix::Link (nullptr, nullptr, 0), false);
    
    
    //std::cout << this << " wants to go to " << SocialMatrix::sig(optiWeight) << "\n";
    
    return std::make_tuple(true_relations[{target_i, 0}], scope_relations[{target_i, 0}]->second, SocialMatrix::Link (true_relations[{target_i, 0}]->first, true_relations[{target_i, 0}]->second, optiWeight), ((target == nullptr && grad[{target_i, 0}] <= 0) || (target != nullptr && optiWeight > alpha[{target_i, 0}])));
    
}

bool Individual::takeAction(){
    std::tuple<SocialMatrix::Link*, Individual*, SocialMatrix::Link, bool> prefAction = Individual::preprocessTakeAction();
    if (std::get<3>(prefAction) && std::get<0>(prefAction) != nullptr){
        this->world->editLink(std::get<0>(prefAction), std::get<2>(prefAction).weight, true, false);
        this->remember(std::get<1>(prefAction), true);
        return true;
    }else if (std::get<1>(prefAction) != nullptr) {
        bool attempt = std::get<1>(prefAction)->responseToAction(this, std::get<2>(prefAction).weight);
        this->remember(std::get<1>(prefAction), attempt);
        return attempt;
    }
    return false;
}

bool Individual::responseToAction(Individual* from, float new_weight){
    std::tuple<SocialMatrix::Link*, Individual*, SocialMatrix::Link, bool> prefAction = Individual::preprocessTakeAction(from);
    if (std::get<0>(prefAction) != nullptr){
        if (SocialMatrix::sig(std::get<2>(prefAction).weight) > 0 && std::get<3>(prefAction)){
            this->world->editLink(std::get<0>(prefAction), std::min(new_weight, std::get<2>(prefAction).weight), true, false);
            return true;
        }
        // We allow to reduce an unwanted friendship
        if (!std::get<3>(prefAction) && SocialMatrix::sig(std::get<2>(prefAction).weight) < std::get<0>(prefAction)->weight){
            this->world->editLink(std::get<0>(prefAction), std::get<2>(prefAction).weight, true, false);
            return true;
        }
        this->world->editLink(std::get<0>(prefAction), new_weight, false, false);
        return false;
    }
    throw std::runtime_error("Unable to find the individual that has asked for the relation");
    return false;
}

void Individual::remember(Individual* target, bool accepted) {
    if (MEMORY_SIZE > 0){
        memoryBuffer.push_front(std::make_pair(target, accepted));
        if (memoryBuffer.size() > MEMORY_SIZE)
            memoryBuffer.pop_back();
    }
}

std::vector<double> Individual::softmax(const akml::DynamicMatrix<float>& logits) {
    std::vector<double> softmax_probs(logits.getNRows());
    double max_logit = *std::max_element(logits.getStorage(), logits.getStorageEnd());

    // Numerator of softmax
    double sum_exp = 0.0;
    for (std::size_t i(0); i < logits.getNRows(); ++i) {
        softmax_probs[i] = std::exp(logits[{i,0}] - max_logit);
        sum_exp += softmax_probs[i];
    }

    // Denumerator
    for (std::size_t i(0); i < softmax_probs.size(); ++i)
        softmax_probs[i] /= sum_exp;

    return softmax_probs;
}

std::size_t Individual::sampleFromSoftmax(const std::vector<double>& softmaxProbs, std::mt19937& gen) {
    // Create a random number generator
    std::uniform_real_distribution<> dis(0.f, 1.f);

    double random_value = dis(gen);
    double cumulative_sum = 0.f;

    for (std::size_t i(0); i < softmaxProbs.size(); ++i) {
        cumulative_sum += softmaxProbs[i];
        if (random_value < cumulative_sum)
            return i;
    }
    // Return the last index if none is selected (due to floating-point issues)
    return softmaxProbs.size() - 1;
}
