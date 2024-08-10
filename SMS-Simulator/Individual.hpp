//
//  Individual.hpp
//  Social MatrixSocialMatrix Simulation
//
//  Created by SMS Associates on 21/03/2023.
//
#ifndef Individual_hpp
#define Individual_hpp

#include "Constants.hpp"
#include "SocialMatrix.hpp"
#include "UtilityFunction.hpp"

#include <cmath>    // for std::exp
#include <numeric>  // for std::accumulate

class SocialMatrix;

class Individual {
	protected:
        akml::DynamicMatrix<float> P;
        SocialMatrix *const world;
        UtilityFunction *utilityFunc;
        std::deque<std::pair<Individual*, bool>> memoryBuffer;

	public:
        static inline unsigned short int GREEDY_SHARE = 0;
        static inline unsigned short int GREEDY_FREQ = 10;
        static inline float DEFAULT_DELTA = 2;
        static inline float DEFAULT_KAPPA = 10;
        static inline float GAMMA_MEAN = 9;
        static inline float GAMMA_DISP = 0.35;
        static inline bool HETEROGENEOUS_P = false;
        static inline unsigned short int MEMORY_SIZE = 0;
        static inline unsigned short int P_DIMENSION = 100;
        
        bool is_greedy;
        float gamma;
        float delta;
        float kappa;
        const unsigned long int agentid;
    
		Individual(SocialMatrix& world, unsigned long int agentid);
        ~Individual();
        bool takeAction();
        bool responseToAction(Individual* from, float new_weight);
        akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> getRelations();
        std::vector<SocialMatrix::Link> getScope();
        float computeUtility(akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>* relations);
        void remember(Individual* target, bool accepted);
    
        inline const UtilityFunction* getUtilityFunction() const{ return utilityFunc; }
        inline akml::DynamicMatrix<float>& getP() { return Individual::P; }
    
    protected:        
        akml::DynamicMatrix<float> computeUtilityGrad(const akml::DynamicMatrix<float>& P_prod, const akml::DynamicMatrix<float>& alpha);
        std::tuple<SocialMatrix::Link*, Individual*, SocialMatrix::Link, bool> preprocessTakeAction(Individual* target=nullptr);
        static std::vector<double> softmax(const akml::DynamicMatrix<float>& logits);
        static std::size_t sampleFromSoftmax(const std::vector<double>& softmaxProbs, std::mt19937& gen);
};
#endif /* Individual_hpp */
