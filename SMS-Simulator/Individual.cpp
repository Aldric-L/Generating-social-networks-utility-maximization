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
        Individual::P[i] = (distribution(gen)==1);
    }
}

std::array<bool, P_DIMENSION> Individual::getP(){
	return Individual::P;
}

std::array<SexualMarket::Link*, GRAPH_SIZE-1> Individual::getRelations(){
    return this->world->getIndividualRelations(this);
}

std::vector<SexualMarket::Link> Individual::getScope() {
    return this->world->getIndividualScope(this);
}
