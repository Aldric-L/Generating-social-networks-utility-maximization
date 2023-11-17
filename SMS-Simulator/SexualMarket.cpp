//
//  SexualMarket.cpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//
#include "SexualMarket.hpp"
#include "Individual.hpp"

SexualMarket::SexualMarket(){
    for (int indiv(0); indiv<GRAPH_SIZE; indiv++){
        individuals[indiv] = new Individual(*this);
    }
    int link_i(0);
    for (int indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (int indiv_t(0); indiv_t<indiv; indiv_t++){
            Link l;
            l.first = individuals[indiv];
            l.second = individuals[indiv_t];
            l.weight = 0.f;
            SexualMarket::links[link_i] = l;
            link_i++;
        }
    }
}

SexualMarket::~SexualMarket(){
    for (int indiv(0); indiv<GRAPH_SIZE; indiv++){
        delete individuals[indiv];
    }
}

void SexualMarket::initializeLinks(){
    int link_i(0);
    for (int indiv(0); indiv<GRAPH_SIZE; indiv++){
        for (int indiv_t(0); indiv_t<indiv; indiv_t++){
            if (indiv == indiv_t+1)
                SexualMarket::links[link_i].weight = 0.1;
            else if (indiv == indiv_t+7)
                SexualMarket::links[link_i].weight = 0.5;
            link_i++;
        }
    }
    
}

std::array<SexualMarket::Link*, GRAPH_SIZE-1> SexualMarket::getIndividualRelations(Individual* indiv) {
    std::array<Link*, GRAPH_SIZE-1> linksForIndividuals;
    unsigned short int incr = 0;
    for (int i(0); i < LINKS_NB; i++){
        if (SexualMarket::links[i].first == indiv || SexualMarket::links[i].second == indiv){
            incr++;
            linksForIndividuals[incr] = &SexualMarket::links[i];
        }
    }
    if (incr != GRAPH_SIZE-1)
        throw std::invalid_argument("Something went really wrong");
    return linksForIndividuals;
}

std::array<Individual*, GRAPH_SIZE> SexualMarket::getIndividuals() {
    return SexualMarket::individuals;
}
