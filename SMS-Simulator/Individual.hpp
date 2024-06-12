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

class SocialMatrix;

class Individual {
	protected:
        struct PSAndAlphaTuple {
            akml::DynamicMatrix<float> P_S; // the union of P vectors of individuals in scope
            akml::DynamicMatrix<float> alpha; // a column vector of weights for each individual in scope (same order of P_S)
            akml::DynamicMatrix<Individual*> beta; // a column vector of pointers to the individuals in the scope
            akml::DynamicMatrix<SocialMatrix::Link*> eta; // a column vector of pointers to the relations in the scope
        };
    
        akml::Matrix<float, P_DIMENSION, 1> P;
        SocialMatrix *const world;
        UtilityFunction *utilityFunc;
        std::deque<Individual*> memoryBuffer;
        std::mt19937 gen;
    
        PSAndAlphaTuple buildPSAndAlpha (const akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>& relations);
        akml::DynamicMatrix<float> computeUtilityGrad(akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>& relations, PSAndAlphaTuple& PS_Alpha);
        std::tuple<SocialMatrix::Link*, Individual*, SocialMatrix::Link, bool> preprocessTakeAction(Individual* target=nullptr);
        

	public:
        static inline unsigned short int GREEDY_SHARE = 0;
        static inline unsigned short int GREEDY_FREQ = 10;
        static inline float DEFAULT_DELTA = 2;
        static inline float DEFAULT_KAPPA = 10;
        static inline float GAMMA_MEAN = 9;
        static inline bool HETEROGENEOUS_P = true;
        static inline unsigned short int MEMORY_SIZE = 0;
        
        bool is_greedy;
        float gamma;
        float delta;
        float kappa;
        const unsigned long int agentid;
    
		Individual(SocialMatrix& world, unsigned long int agentid);
        ~Individual();
        akml::Matrix<float, P_DIMENSION, 1>& getP();
        bool takeAction();
        bool responseToAction(Individual* from, float new_weight);
        akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> getRelations();
        std::vector<SocialMatrix::Link> getScope();
        float computeUtility(akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>* relations);
        inline const UtilityFunction* getUtilityFunction() const{ return utilityFunc; }
};
#endif /* Individual_hpp */
