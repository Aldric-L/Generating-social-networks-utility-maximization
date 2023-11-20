//
//  SexualMarket.hpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//


#ifndef SexualMarket_hpp
#define SexualMarket_hpp

#include "Constants.hpp"
class Individual;

class SexualMarket {
    public:
        struct Link {
            Individual* first;
            Individual* second;
            float weight;
        };
        
    protected:
        std::array<Individual*, GRAPH_SIZE> individuals;
        std::array<Link, LINKS_NB> links;

    public:
        SexualMarket();
        ~SexualMarket();    
    
        //temporaire, pour l'instant l'initialisation est hardcoded
        void initializeLinks();
        std::array<Link*, GRAPH_SIZE-1> getIndividualRelations(Individual* indiv);
        std::array<Individual*, GRAPH_SIZE> getIndividuals();
        std::vector<SexualMarket::Link> getIndividualScope(Individual* indiv);
};

 /* SexualMarket_hpp */
#endif
