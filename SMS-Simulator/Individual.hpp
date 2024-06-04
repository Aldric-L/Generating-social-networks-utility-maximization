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

class SocialMatrix;

class Individual {
	protected:
        typedef std::tuple<akml::DynamicMatrix<float>, akml::DynamicMatrix<float>, akml::DynamicMatrix<Individual*>, akml::DynamicMatrix<SocialMatrix::Link*>> PSAndAlphaTuple;
    
        akml::Matrix<float, P_DIMENSION, 1> P;
        SocialMatrix *world;
        //std::mt19937 gen;
    
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
        
        bool is_greedy;
        float gamma;
        float delta;
        float kappa;
        unsigned long int agentid;
    
		Individual(SocialMatrix& world, unsigned long int agentid);
        akml::Matrix<float, P_DIMENSION, 1>& getP();
        bool takeAction();
        bool responseToAction(Individual* from, float new_weight);
        akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1> getRelations();
        std::vector<SocialMatrix::Link> getScope();
        float computeUtility(akml::Matrix<SocialMatrix::Link*, GRAPH_SIZE-1, 1>* relations);
};
#endif /* Individual_hpp */
