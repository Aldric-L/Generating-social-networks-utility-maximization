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

float Individual::computeUtility() {
    std::array<SexualMarket::Link*, GRAPH_SIZE-1> relations = Individual::getRelations();
    unsigned short int effectiveRelationCount = 0;
    std::array<akml::Matrix<float, P_DIMENSION, 1>, GRAPH_SIZE-1> P_S_temp;
    akml::Matrix<float, GRAPH_SIZE-1, 1> alpha;
    float RHS;

    
    // In order to try to avoid using dynamic matrices, we will keep fixed sized matrices with 0 where we should not have a column
    for (int i(0); i < GRAPH_SIZE-1; i++){
        if (relations[i]->weight > 0){
            effectiveRelationCount++;
            if (relations[i]->first == this){
                P_S_temp[i] = relations[i]->second->getP();
            }else {
                P_S_temp[i] = relations[i]->first->getP();
            }
            RHS += std::pow(relations[i]->weight, Individual::gamma);
        }
        alpha(i+1, 1) = relations[i]->weight;
    }
    RHS = std::pow(RHS, 2);
    
    akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1> P_S (P_S_temp);
    akml::Matrix<float, GRAPH_SIZE-1, P_DIMENSION> P_S_transpose = akml::transpose(P_S);
    akml::Matrix<float, GRAPH_SIZE-1, 1> P_prod = akml::matrix_product(P_S_transpose, Individual::P);
    
    akml::cout_matrix(P_S);
    akml::cout_matrix(alpha);
    akml::cout_matrix(P_prod);
    
    P_prod.transform([](float val) { return val/P_DIMENSION; });
    
    akml::cout_matrix(P_prod);
    float LHS = akml::inner_product(alpha, P_prod);
        
    std::cout << "Utility : LHS=" << LHS << " RHS=" << RHS << " Total=" << LHS-RHS << "(Gamma=" << Individual::gamma << ")" << std::endl;
    return LHS-RHS;
}
