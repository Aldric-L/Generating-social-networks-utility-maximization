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
            float compatibility = -MAXFLOAT;
            
            inline Link() {};
            inline Link(Individual* first, Individual* second, float weight, float compatibility) : first(first), second(second), weight(weight), compatibility(compatibility) {};
            Link(Individual* first, Individual* second, float weight);
        };
        
    protected:
        // We replace std::array<Link, LINKS_NB> links; by std::vector in order to register the Links in the heap and not in the stack. This is, in fact, very sad. I love static things. I grieve the static array.
        std::vector<Link> links;
        akml::Matrix<Individual*, GRAPH_SIZE, 1> individuals;
        typedef akml::Save<7, std::size_t, unsigned long int, unsigned long int, float, float, bool, bool> EdgeSaveTrackerType;
        typedef akml::Save<3, std::size_t, unsigned long int, float> UtilitySaveTrackerType;
        typedef akml::Save<8, std::size_t, unsigned long int, float, bool, float, float, std::size_t, std::string> VerticesSaveTrackerType;
        typedef akml::MatrixSave<GRAPH_SIZE+1, float> ClusteringSaveTrackerType;
        typedef akml::FullMatrixSave<akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE>> FinalSaveTrackerType;
        akml::CSV_Saver<EdgeSaveTrackerType> edgeTrackersManager;
        akml::CSV_Saver<UtilitySaveTrackerType> utilityTrackersManager;
        akml::CSV_Saver<VerticesSaveTrackerType> verticesTrackersManager;
        akml::CSV_Saver<ClusteringSaveTrackerType> clusteringTrackersManager;
        akml::CSV_Saver<FinalSaveTrackerType> finalAdjacencyMatrixTrackersManager;
        std::mt19937 gen;
        std::string logPath = "";
        std::string logID = "";

    public:
        static inline bool SHOULD_I_LOG = true;
        static inline bool COMPUTE_CLUSTERING = true;
        static inline bool COMPUTE_CLEARING = true;
        static inline std::string GLOBAL_LOG_PREFIX = "";
        static inline bool MODE_FOLDER_LOG = true;
        static inline unsigned short int UTILITY_COMPUTATION_INTERVAL = 10;
        static inline std::size_t SCOPE_DEPTH = 1;
    
        std::size_t currentRound = 0;
        SocialMatrix();
        SocialMatrix(const akml::DynamicMatrix<float>& compatibilityMatrix);
        ~SocialMatrix();
    
        void initializeLinks();
        void initializeLinks(const akml::DynamicMatrix<float>& adjacencyMatrix);
        void forceEditCompatibilityMatrix(const akml::DynamicMatrix<float>& compatibilityMatrix);
        bool checkLoveTriangleCondition();
        akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> getIndividualRelations(Individual* indiv);
        akml::Matrix<Individual*, GRAPH_SIZE, 1> getIndividuals();
        Individual* getIndividual(const std::size_t indiv_id);
        std::vector<SocialMatrix::Link> getIndividualScope(Individual* indiv, Individual* original=nullptr, const std::size_t scopeDepth=0);
        void editLink(const Individual* indiv1, const Individual* indiv2, const float newWeight, bool accepted=true, bool forced=false);
        void editLink(SocialMatrix::Link* link, const float newWeight, bool accepted=true, bool forced=false);
        unsigned int processARound(std::size_t totalrounds=0);
        akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> asAdjacencyMatrix() const;
        akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE> asBinaryAdjacencyMatrix(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE>* adjacencymatrix = nullptr) const;
        akml::Matrix<std::size_t, GRAPH_SIZE, GRAPH_SIZE> computeDegreesOfSeparation(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix = nullptr) const;
        akml::Matrix<float, GRAPH_SIZE+1, 1> computeClusteringCoefficients(akml::Matrix<bool, GRAPH_SIZE, GRAPH_SIZE>* binaryadjacencymatrix = nullptr) const;
        std::pair<std::string, std::string> whereWillYouLog() const;
        inline std::mt19937& getRandomGen() { return gen; }
        float getCompatibilityBtwnIndividuals(const Individual* indiv1, const Individual* indiv2);
        Link* findRelation(const Individual* indiv1, const Individual* indiv2);
};

 /* SocialMatrix_hpp */
#endif
