//
//  SexualMarket.hpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//


#ifndef SexualMarket_hpp
#define SexualMarket_hpp

#include "Constants.hpp"
#include "AKML-lib/AgentBasedUtilities/Save.hpp"
#include "AKML-lib/AgentBasedUtilities/CSV_Saver.hpp"

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
        int currentRound = 0;
        std::array<Individual*, GRAPH_SIZE> individuals;
        std::array<Link, LINKS_NB> links;
        typedef akml::Save<5, int, unsigned short int, unsigned short int, float, float> EdgeSaveTrackerType;
        typedef akml::Save<3, int, unsigned short int, float> UtilitySaveTrackerType;
        akml::CSV_Saver<EdgeSaveTrackerType> edgeTrackersManager;
        akml::CSV_Saver<UtilitySaveTrackerType> utilityTrackersManager;

    public:
        SexualMarket();
        ~SexualMarket();    
    
        // WIP - for now initialization is hardcoded
        void initializeLinks();
        std::array<Link*, GRAPH_SIZE-1> getIndividualRelations(Individual* indiv);
        std::array<Individual*, GRAPH_SIZE> getIndividuals();
        std::vector<SexualMarket::Link> getIndividualScope(Individual* indiv);
        void editLink(Individual* indiv1, Individual* indiv2, float newWeight);
        void editLink(SexualMarket::Link* link, float newWeight);

};

 /* SexualMarket_hpp */
#endif
