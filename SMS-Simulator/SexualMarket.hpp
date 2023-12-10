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
        // We replace std::array<Link, LINKS_NB> links; by std::vector in order to register the Links in the heap and not in the stack
        // This is, in fact, very sad. I love static things. I grief the static array (no way to use my static Matrices with objects).
        std::vector<Link> links;
        std::size_t currentRound = 0;
        akml::Matrix<Individual*, GRAPH_SIZE, 1> individuals;
        typedef akml::Save<6, std::size_t, unsigned long int, unsigned long int, float, float, bool> EdgeSaveTrackerType;
        typedef akml::Save<3, std::size_t, unsigned long int, float> UtilitySaveTrackerType;
        typedef akml::Save<4, std::size_t, unsigned long int, float, std::string> VerticesSaveTrackerType;
        akml::CSV_Saver<EdgeSaveTrackerType> edgeTrackersManager;
        akml::CSV_Saver<UtilitySaveTrackerType> utilityTrackersManager;
        akml::CSV_Saver<VerticesSaveTrackerType> verticesTrackersManager;

    public:
        static inline bool SHOULD_I_LOG = true;
        SexualMarket();
        ~SexualMarket();    
    
        // WIP - for now initialization is hardcoded
        void initializeLinks();
        akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1> getIndividualRelations(Individual* indiv);
        akml::Matrix<Individual*, GRAPH_SIZE, 1> getIndividuals();
        Individual* getIndividual(const std::size_t indiv_id);
        std::vector<SexualMarket::Link> getIndividualScope(Individual* indiv);
        void editLink(Individual* indiv1, Individual* indiv2, float newWeight, bool accepted=true);
        void editLink(SexualMarket::Link* link, float newWeight, bool accepted=true);
        unsigned int processARound();
        akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> asAdjacencyMatrix();


};

 /* SexualMarket_hpp */
#endif
