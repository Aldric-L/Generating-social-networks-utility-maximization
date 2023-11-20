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
    std::array<Link*, GRAPH_SIZE-1> linksForIndividual;
    unsigned short int incr = 0;
    for (int i(0); i < LINKS_NB; i++){
        if (SexualMarket::links[i].first == indiv || SexualMarket::links[i].second == indiv){
            linksForIndividual[incr] = &SexualMarket::links[i];
            incr++;
        }
    }
    if (incr != GRAPH_SIZE-1)
        throw std::invalid_argument("Something went really wrong");
    return linksForIndividual;
}

std::array<Individual*, GRAPH_SIZE> SexualMarket::getIndividuals() {
    return SexualMarket::individuals;
}

std::vector<SexualMarket::Link> SexualMarket::getIndividualScope(Individual* indiv) {
    std::vector<SexualMarket::Link> scope;
    std::array<Link*, GRAPH_SIZE-1> linksofIndividual = SexualMarket::getIndividualRelations(indiv);
    for (int level0(0); level0 < GRAPH_SIZE-1; level0++){
        if (linksofIndividual[level0]->weight > 0){
            Individual* target0 = nullptr;
            Individual* target1 = nullptr;
            
            if (linksofIndividual[level0]->first == indiv)
                target0 = linksofIndividual[level0]->second;
            else
                target0 = linksofIndividual[level0]->first;
            
            std::array<Link*, GRAPH_SIZE-1> linksofTarget = SexualMarket::getIndividualRelations(target0);
            for (int level1(0); level1 < GRAPH_SIZE-1; level1++){
                if (linksofTarget[level1]->weight > 0){
                    if (linksofTarget[level1]->first == target0)
                        target1 = linksofTarget[level1]->second;
                    else
                        target1 = linksofTarget[level1]->first;
                    
                    if (target0 != target1 && target1 != indiv && target0 != indiv){
                        Link l;
                        l.first = indiv;
                        l.second = target1;
                        l.weight = linksofTarget[level1]->weight * linksofIndividual[level0]->weight;
                        scope.push_back(l);
                    }
                }
            }
            
        }
    }
    return scope;
}
