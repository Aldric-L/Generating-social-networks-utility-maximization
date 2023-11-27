//
//  Individual.cpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//

#include "Individual.hpp"
#include "SexualMarket.hpp"

Individual::Individual(SexualMarket& world, int agentid){
	this->world = &world;
    Individual::agentid = agentid;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned short int> distribution(0,1);
    for (int i(0); i < P_DIMENSION; i++){
        Individual::P(i+1,1) = (distribution(gen)==1) ? 1 : 0;
    }
    std::normal_distribution<float> norm(0.5,0.15);
    float g = 0;
    while((g = norm(gen)) > 0.9 || g < 0.1){ g = norm(gen); }
    Individual::gamma = g+1;
}

akml::Matrix<float, P_DIMENSION, 1> Individual::getP(){
	return Individual::P;
}

std::array<SexualMarket::Link*, GRAPH_SIZE-1> Individual::getRelations(){
    return this->world->getIndividualRelations(this);
}

std::vector<SexualMarket::Link> Individual::getScope() {
    return this->world->getIndividualScope(this);
}

Individual::PSAndAlphaTuple Individual::buildPSAndAlpha(const std::array<SexualMarket::Link*, GRAPH_SIZE-1>& relations){
    std::array<akml::Matrix<float, P_DIMENSION, 1>, GRAPH_SIZE-1> P_S_temp;
    // Weights matrix
    akml::Matrix<float, GRAPH_SIZE-1, 1> alpha;
    // Pointers to individual matrix
    akml::Matrix<Individual*, GRAPH_SIZE-1, 1> beta;
    // Pointers to relation matrix
    akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1> eta;

    // In order to try to avoid using dynamic matrices, we will keep fixed sized matrices with 0 where we should not have a column
    for (int i(0); i < GRAPH_SIZE-1; i++){
        if (relations[i]->weight > 0){
            if (relations[i]->first == this){
                P_S_temp[i] = relations[i]->second->getP();
                beta(i+1, 1) = relations[i]->second;
            }else {
                P_S_temp[i] = relations[i]->first->getP();
                beta(i+1, 1) = relations[i]->first;
            }
            eta(i+1, 1) = relations[i];
        }else{
            beta(i+1, 1) = nullptr;
            eta(i+1, 1) = nullptr;
        }
        alpha(i+1, 1) = relations[i]->weight;
    }
    
    akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1> P_S (P_S_temp);
    return std::make_tuple(P_S, alpha, beta, eta);
}

float Individual::computeUtility(std::array<SexualMarket::Link*, GRAPH_SIZE-1>* relations) {
    bool rel_td = false;
    if (relations == nullptr){
        relations = new std::array<SexualMarket::Link*, GRAPH_SIZE-1>;
        *relations = this->getRelations();
        rel_td = true;
    }
    unsigned short int effectiveRelationCount = 0;
    std::array<akml::Matrix<float, P_DIMENSION, 1>, GRAPH_SIZE-1> P_S_temp;
    akml::Matrix<float, GRAPH_SIZE-1, 1> alpha;
    float RHS;

    
    // In order to try to avoid using dynamic matrices, we will keep fixed sized matrices with 0 where we should not have a column
    for (int i(0); i < GRAPH_SIZE-1; i++){
        if ((*relations)[i]->weight > 0){
            effectiveRelationCount++;
            if ((*relations)[i]->first == this){
                P_S_temp[i] = (*relations)[i]->second->getP();
            }else {
                P_S_temp[i] = (*relations)[i]->first->getP();
            }
            RHS += std::pow((*relations)[i]->weight, Individual::gamma);
        }
        alpha(i+1, 1) = (*relations)[i]->weight;
    }
    RHS = std::pow(RHS, 1/(Individual::gamma))/(GRAPH_SIZE);
    
    akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1> P_S (P_S_temp);
    akml::Matrix<float, GRAPH_SIZE-1, P_DIMENSION> P_S_transpose = akml::transpose(P_S);
    akml::Matrix<float, GRAPH_SIZE-1, 1> P_prod = akml::matrix_product(P_S_transpose, Individual::P);
    
    P_prod.transform([](float val) { return val/P_DIMENSION; });
    
    float LHS = akml::inner_product(alpha, P_prod) * 4;
        
    std::cout << "Utility : LHS=" << LHS << " RHS=" << RHS << " Total=" << LHS-RHS << "(Gamma=" << Individual::gamma << ")" << std::endl;
    
    if (rel_td)
        delete relations;
    
    return LHS-RHS;
}

akml::Matrix<float, GRAPH_SIZE-1, 1> Individual::computeUtilityGrad(std::array<SexualMarket::Link*, GRAPH_SIZE-1>* relations, Individual::PSAndAlphaTuple* PS_Alpha) {
    // The following lines are implemented to avoid memory leaks : pointers allow us to use this function without providing
    // any argument. In fact, it has been done that way to allow us to not recalculate variables that have been already used
    // in a previous calculation. By pointers, we can take back our calculations where we were.
    bool rel_td = false;
    bool ps_td = false;
    if (relations == nullptr){
        relations = new std::array<SexualMarket::Link*, GRAPH_SIZE-1>;
        *relations = this->getRelations();
        rel_td = true;
    }
    if (PS_Alpha == nullptr){
        PS_Alpha = new Individual::PSAndAlphaTuple;
        *PS_Alpha = this->buildPSAndAlpha(*relations);
        ps_td = true;
    }
    akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1> PS_temp = std::get<0>(*PS_Alpha);
    akml::Matrix<float, GRAPH_SIZE-1, 1> Alpha_temp = std::get<1>(*PS_Alpha);
    akml::Matrix<float, GRAPH_SIZE-1, 1> grad;
    akml::Matrix<float, GRAPH_SIZE-1, P_DIMENSION> P_S_transpose = akml::transpose(PS_temp);
    akml::Matrix<float, GRAPH_SIZE-1, 1> P_prod = akml::matrix_product(P_S_transpose, this->P);
    P_prod.transform([](float val) { return (val*4)/P_DIMENSION; });
    
    float scalaralpha = 0;
    for (unsigned short int line(1); line <= GRAPH_SIZE-1; line++){
        scalaralpha += (Alpha_temp(line, 1) != 0) ? std::pow(Alpha_temp(line, 1),this->gamma)/GRAPH_SIZE : 0;
    }
    
    Alpha_temp.transform([this, &scalaralpha](float val) { return (val != 0) ? (-1/(Individual::gamma))*(this->gamma)*std::pow(val, ((this->gamma)-1)/(Individual::gamma))*scalaralpha: 0; });
    grad = P_prod + Alpha_temp;
    //akml::cout_matrix(P_prod);
    //akml::cout_matrix(std::get<1>(PS_ALPHA));
    akml::cout_matrix(grad);
    
    if (rel_td)
        delete relations;
    if (ps_td)
        delete PS_Alpha;
    return grad;
}

std::tuple<SexualMarket::Link*, Individual*, SexualMarket::Link, bool> Individual::preprocessTakeAction(Individual* target){
    std::array<SexualMarket::Link*, GRAPH_SIZE-1> relations = this->getRelations();
    std::array<SexualMarket::Link*, GRAPH_SIZE-1> rel_temp;
    std::vector<SexualMarket::Link> scope = this->getScope();
    
    // Because we want to distinguish individuals that are in the scope but with whom there is no link and those  no link
    // with whom there is no link at all, we create a fake little link of 0.00001 for individuals in scope.
    // Pointers to relation matrix
    akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1> true_eta;
    for (std::size_t rel(0); rel < GRAPH_SIZE-1; rel++){
        rel_temp[rel] = new SexualMarket::Link;
        *rel_temp[rel] = *relations[rel];
        if (relations[rel]->weight == 0){
            for (std::size_t i(0); i < scope.size(); i++){
                if ((scope[i].first == this && relations[rel]->first == this && scope[i].second == relations[rel]->second)
                    || (scope[i].first == this && relations[rel]->second == this && scope[i].second == relations[rel]->first)
                    || (scope[i].second == this && relations[rel]->second == this && scope[i].first == relations[rel]->first)
                    || (scope[i].second == this && relations[rel]->first == this && scope[i].first == relations[rel]->second)){
                    rel_temp[rel]->weight = 0.00001;
                    break;
                }
            }
        }
        true_eta(rel+1, 1) = relations[rel];
    }

    Individual::PSAndAlphaTuple PS_Alpha = this->buildPSAndAlpha(std::ref(rel_temp));
    
    //We delete temp pointers
    for (std::size_t rel(0); rel < GRAPH_SIZE-1; rel++){
        delete rel_temp[rel];
    }
    
    std::cout << "Computing grad" << std::endl;
    akml::Matrix<float, GRAPH_SIZE-1, 1> grad = Individual::computeUtilityGrad(&relations, &PS_Alpha);
    std::cout << "Now lets reduce the grad to accessible individuals (scope=" << scope.size() << ")" << std::endl;
    for (std::size_t i(0); i < GRAPH_SIZE-1; i++){
        // You can't select an individual that you don't see (no freehate)
        if (std::get<2>(PS_Alpha)(i+1, 1) == nullptr){
            grad(i+1, 1) = 0;
        // You can't reduce a friendship that does not exist
        }else if (grad(i+1, 1) < 0 && std::get<1>(PS_Alpha)(i+1, 1) <= 0.0001 ){
            grad(i+1, 1) = 0;
        }
    }
    akml::cout_matrix(grad);
    //akml::cout_matrix(std::get<2>(PS_Alpha));
    if (target == nullptr){
        unsigned short int max_i = akml::arg_max(grad, true);
        if (std::abs(grad(max_i+1, 1)) < 0.001){
            std::cout << "Want to do an unsignifiant action. Aborted.";
            SexualMarket::Link newlinkwanted (nullptr, nullptr, 0);
            return std::make_tuple(nullptr, nullptr, newlinkwanted, false);
        }
        float mov = std::max(0.f, ((std::get<1>(PS_Alpha)(max_i+1, 1) <= 0.0001) ? 0.f : std::get<1>(PS_Alpha)(max_i+1, 1)) + grad(max_i+1, 1));
        std::cout << "Want to take action on the coef " << max_i << " in the direction of " << grad(max_i+1, 1) << " to reach " << mov << " which corresponds currently to a link of " << ((std::get<1>(PS_Alpha)(max_i+1, 1) <= 0.0001) ? 0 : std::get<1>(PS_Alpha)(max_i+1, 1)) << " which links to the individual " << std::get<2>(PS_Alpha)(max_i+1, 1);
        //std::cout << "We verifiy true_eta : " << true_eta(max_i+1, 1)->first << " - " << true_eta(max_i+1, 1)->second << " - " << true_eta(max_i+1, 1)->weight << "\n";
        
        SexualMarket::Link newlinkwanted (true_eta(max_i+1, 1)->first, true_eta(max_i+1, 1)->second, mov);
        
        return std::make_tuple(true_eta(max_i+1, 1), std::get<2>(PS_Alpha)(max_i+1, 1), newlinkwanted, (grad(max_i+1, 1) <= 0));
    }else {
        short int target_i = -1;
        for (unsigned short int i(0); i < GRAPH_SIZE-1; i++){
            if (std::get<2>(PS_Alpha)(i+1, 1) == target){
                target_i = i;
                break;
            }
        }
        
        // target not found
        if (target_i == -1){
            SexualMarket::Link newlinkwanted (nullptr, nullptr, 0);
            return std::make_tuple(nullptr, nullptr, newlinkwanted, false);
        }
        
        float mov = std::max(0.f, ((std::get<1>(PS_Alpha)(target_i+1, 1) <= 0.0001) ? 0.f : std::get<1>(PS_Alpha)(target_i+1, 1)) + grad(target_i+1, 1));
        SexualMarket::Link newlinkwanted (true_eta(target_i+1, 1)->first, true_eta(target_i+1, 1)->second, mov);
        
        return std::make_tuple(true_eta(target_i+1, 1), std::get<2>(PS_Alpha)(target_i+1, 1), newlinkwanted, (grad(target_i+1, 1) <= 0));
        
    }
    
}

void Individual::takeAction(){
    std::tuple<SexualMarket::Link*, Individual*, SexualMarket::Link, bool> prefAction = Individual::preprocessTakeAction();
    if (std::get<3>(prefAction) && std::get<0>(prefAction) != nullptr){
        this->world->editLink(std::get<0>(prefAction), std::get<2>(prefAction).weight, true);
    }else if (std::get<1>(prefAction) != nullptr) {
        std::get<1>(prefAction)->responseToAction(this, std::get<2>(prefAction).weight);
    }
}

void Individual::responseToAction(Individual* from, float new_weight){
    std::tuple<SexualMarket::Link*, Individual*, SexualMarket::Link, bool> prefAction = Individual::preprocessTakeAction(from);
    if (std::get<0>(prefAction) != nullptr && std::get<2>(prefAction).weight > 0){
        std::cout << "Ask to respond to action of " << from << " but we move to " << std::min(new_weight, std::get<2>(prefAction).weight) << std::endl;
        this->world->editLink(std::get<0>(prefAction), std::min(new_weight, std::get<2>(prefAction).weight), true);
    }else if (std::get<0>(prefAction) != nullptr) {
        if (std::get<2>(prefAction).weight <= 0){
            std::cout << "Ask to respond to action of " << from << " but our desire is negative or null (" << std::get<2>(prefAction).weight << ")" << std::endl;
            this->world->editLink(std::get<0>(prefAction), std::min(new_weight, std::get<2>(prefAction).weight), false);
        }else{
            std::cout << "Ask to respond to action of " << from << " but he is not found" << std::endl;
        }
    }
}
