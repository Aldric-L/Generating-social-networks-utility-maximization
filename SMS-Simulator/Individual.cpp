//
//  Individual.cpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//

#include "Individual.hpp"
#include "SexualMarket.hpp"

Individual::Individual(SexualMarket& world){
	this->world = &world;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned short int> distribution(0,1);
    for (int i(0); i < P_DIMENSION; i++){
        //Individual::P[i] = (distribution(gen)==1);
        Individual::P(i+1,1) = (distribution(gen)==1) ? 1 : 0;
    }
    std::normal_distribution<float> norm(0,1);
    Individual::gamma = std::abs(norm(gen));
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

std::tuple<akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1>,akml::Matrix<float, GRAPH_SIZE-1, 1>, akml::Matrix<Individual*, GRAPH_SIZE-1, 1>> Individual::buildPSAndAlpha(const std::array<SexualMarket::Link*, GRAPH_SIZE-1>& relations){
    std::array<akml::Matrix<float, P_DIMENSION, 1>, GRAPH_SIZE-1> P_S_temp;
    akml::Matrix<float, GRAPH_SIZE-1, 1> alpha;
    akml::Matrix<Individual*, GRAPH_SIZE-1, 1> beta;

    // In order to try to avoid using dynamic matrices, we will keep fixed sized matrices with 0 where we should not have a column
    for (int i(0); i < GRAPH_SIZE-1; i++){
            if (relations[i]->first == this){
                if (relations[i]->weight > 0)
                    P_S_temp[i] = relations[i]->second->getP();
                if (relations[i]->weight > 0)
                    beta(i+1, 1) = relations[i]->second;
                else
                    beta(i+1, 1) = nullptr;
            }else {
                if (relations[i]->weight > 0)
                    P_S_temp[i] = relations[i]->first->getP();
                if (relations[i]->weight > 0)
                    beta(i+1, 1) = relations[i]->first;
                else
                    beta(i+1, 1) = nullptr;
            }
        alpha(i+1, 1) = relations[i]->weight;
    }
    
    akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1> P_S (P_S_temp);
    return std::make_tuple(P_S, alpha, beta);
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
    RHS = std::pow(RHS, 2)/P_DIMENSION;
    
    akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1> P_S (P_S_temp);
    akml::Matrix<float, GRAPH_SIZE-1, P_DIMENSION> P_S_transpose = akml::transpose(P_S);
    akml::Matrix<float, GRAPH_SIZE-1, 1> P_prod = akml::matrix_product(P_S_transpose, Individual::P);
    
    //akml::cout_matrix(P_S);
    //akml::cout_matrix(alpha);
    //akml::cout_matrix(P_prod);
    
    P_prod.transform([](float val) { return val/P_DIMENSION; });
    
    //akml::cout_matrix(P_prod);
    float LHS = akml::inner_product(alpha, P_prod);
        
    std::cout << "Utility : LHS=" << LHS << " RHS=" << RHS << " Total=" << LHS-RHS << "(Gamma=" << Individual::gamma << ")" << std::endl;
    
    if (rel_td)
        delete relations;
    
    return LHS-RHS;
}

akml::Matrix<float, GRAPH_SIZE-1, 1> Individual::computeUtilityGrad(std::array<SexualMarket::Link*, GRAPH_SIZE-1>* relations, std::tuple<akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1>,akml::Matrix<float, GRAPH_SIZE-1, 1>, akml::Matrix<Individual*, GRAPH_SIZE-1, 1>>* PS_Alpha) {
    bool rel_td = false;
    bool ps_td = false;
    if (relations == nullptr){
        relations = new std::array<SexualMarket::Link*, GRAPH_SIZE-1>;
        *relations = this->getRelations();
        rel_td = true;
    }
    if (PS_Alpha == nullptr){
        PS_Alpha = new std::tuple<akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1>,akml::Matrix<float, GRAPH_SIZE-1, 1>, akml::Matrix<Individual*, GRAPH_SIZE-1, 1>>;
        *PS_Alpha = this->buildPSAndAlpha(*relations);
        ps_td = true;
    }
    akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1> PS_temp = std::get<0>(*PS_Alpha);
    akml::Matrix<float, GRAPH_SIZE-1, 1> Alpha_temp = std::get<1>(*PS_Alpha);
    akml::Matrix<float, GRAPH_SIZE-1, 1> grad;
    akml::Matrix<float, GRAPH_SIZE-1, P_DIMENSION> P_S_transpose = akml::transpose(PS_temp);
    akml::Matrix<float, GRAPH_SIZE-1, 1> P_prod = akml::matrix_product(P_S_transpose, this->P);
    P_prod.transform([](float val) { return val/P_DIMENSION; });
    
    float scalaralpha = 0;
    for (unsigned short int line(1); line <= GRAPH_SIZE-1; line++){
        scalaralpha += (Alpha_temp(line, 1) != 0) ? std::pow(Alpha_temp(line, 1),this->gamma)/P_DIMENSION : 0;
    }
    
    Alpha_temp.transform([this, &scalaralpha](float val) { return (val != 0) ? (-2)*Individual::gamma*std::pow(val, Individual::gamma-1)*scalaralpha : 0; });
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


void Individual::takeAction(){
    std::array<SexualMarket::Link*, GRAPH_SIZE-1> relations = this->getRelations();
    std::array<SexualMarket::Link*, GRAPH_SIZE-1> rel_temp;
    std::vector<SexualMarket::Link> scope = this->getScope();
    
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
    }

    std::tuple<akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1>,akml::Matrix<float, GRAPH_SIZE-1, 1>, akml::Matrix<Individual*, GRAPH_SIZE-1, 1>> PS_Alpha = this->buildPSAndAlpha(std::ref(rel_temp));
    
    std::cout << "Computing grad" << std::endl;
    akml::Matrix<float, GRAPH_SIZE-1, 1> grad = Individual::computeUtilityGrad(&relations, &PS_Alpha);
    std::cout << "Now lets reduce the grad to accessible individuals (scope=" << scope.size() << ")" << std::endl;
    for (std::size_t i(0); i < GRAPH_SIZE-1; i++){
        if (std::get<2>(PS_Alpha)(i+1, 1) == nullptr){
            grad(i+1, 1) = 0;
        }
    }
    akml::cout_matrix(grad);
    akml::cout_matrix(std::get<2>(PS_Alpha));
    unsigned short int max_i = akml::arg_max(grad, true);
    std::cout << "Want to take action on the coef " << max_i << " in the direction of " << grad(max_i+1, 1) << " which corresponds currently to a link of " << ((std::get<1>(PS_Alpha)(max_i+1, 1) <= 0.0001) ? 0 : std::get<1>(PS_Alpha)(max_i+1, 1)) << " which links to the individual " << std::get<2>(PS_Alpha)(max_i+1, 1);
}

bool Individual::responseToAction(Individual* from, float new_weight){
}
