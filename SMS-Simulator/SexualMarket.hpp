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
            
            inline Link() {};
            inline Link(Individual* first, Individual* second, float weight) : first(first), second(second), weight(weight) {};
        };
        
    protected:
        std::array<Individual*, GRAPH_SIZE> individuals;
        std::array<Link, LINKS_NB> links;

    public:
        SexualMarket();
        ~SexualMarket();    
    
        // WIP - for now initialization is hardcoded
        void initializeLinks();
        std::array<Link*, GRAPH_SIZE-1> getIndividualRelations(Individual* indiv);
        std::array<Individual*, GRAPH_SIZE> getIndividuals();
        std::vector<SexualMarket::Link> getIndividualScope(Individual* indiv);
        void editLink(Individual* indiv1, Individual* indiv2, float newWeight);
};

 /* SexualMarket_hpp */
#endif
