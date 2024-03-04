//
//  SocialMatrix.hpp
//  Social Matrix Simulation
//
//  Created by SMS Associates on 21/03/2023.
//


#ifndef SocialMatrix_hpp
#define SocialMatrix_hpp

#include "Constants.hpp"
#include "AKML-lib/AgentBasedUtilities/Save.hpp"
#include "AKML-lib/AgentBasedUtilities/CSV_Saver.hpp"

class Individual;

class SocialMatrix {
    public:
        struct Link {
            Individual* first;
            Individual* second;
            float weight;
            
            inline Link() {};
            inline Link(Individual* first, Individual* second, float weight) : first(first), second(second), weight(weight) {};
        };
        
    protected:
        // We replace std::array<Link, LINKS_NB> links; by std::vector in order to register the Links in the heap and not in the stack. This is, in fact, very sad. I love static things. I grief the static array.
        std::vector<Link> links;
        akml::Matrix<Individual*, GRAPH_SIZE, 1> individuals;
        typedef akml::Save<6, std::size_t, unsigned long int, unsigned long int, float, float, bool> EdgeSaveTrackerType;
        typedef akml::Save<3, std::size_t, unsigned long int, float> UtilitySaveTrackerType;
        typedef akml::Save<8, std::size_t, unsigned long int, float, bool, float, float, std::size_t, std::string> VerticesSaveTrackerType;
        typedef akml::MatrixSave<GRAPH_SIZE+1, float> ClusteringSaveTrackerType;
        akml::CSV_Saver<EdgeSaveTrackerType> edgeTrackersManager;
        akml::CSV_Saver<UtilitySaveTrackerType> utilityTrackersManager;
        akml::CSV_Saver<VerticesSaveTrackerType> verticesTrackersManager;
        akml::CSV_Saver<ClusteringSaveTrackerType> clusteringTrackersManager;

    public:
        static inline bool SHOULD_I_LOG = true;
        std::size_t currentRound = 0;
        SocialMatrix();
        ~SocialMatrix();    
    
        void initializeLinks();
        akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> getIndividualRelations(Individual* indiv);
        akml::Matrix<Individual*, GRAPH_SIZE, 1> getIndividuals();
        Individual* getIndividual(const std::size_t indiv_id);
        std::vector<SocialMatrix::Link> getIndividualScope(Individual* indiv);
        void editLink(Individual* indiv1, Individual* indiv2, float newWeight, bool accepted=true);
        void editLink(SocialMatrix::Link* link, float newWeight, bool accepted=true);
        unsigned int processARound(std::size_t totalrounds=0);
        akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> asAdjacencyMatrix();
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> asBinaryAdjacencyMatrix(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE>* adjacencymatrix = nullptr);
        akml::Matrix<std::size_t, GRAPH_SIZE, GRAPH_SIZE> computeDegreesOfSeparation(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix = nullptr);
        akml::Matrix<float, GRAPH_SIZE+1, 1> computeClusteringCoefficients(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix = nullptr);


};

 /* SocialMatrix_hpp */
#endif
