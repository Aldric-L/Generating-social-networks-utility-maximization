//
//  Individual.hpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//
#ifndef Individual_hpp
#define Individual_hpp

#include "Constants.hpp"
#include "SexualMarket.hpp"

class SexualMarket;

class Individual {
	protected:
        typedef std::tuple<akml::DynamicMatrix<float>, akml::DynamicMatrix<float>, akml::DynamicMatrix<Individual*>, akml::DynamicMatrix<SexualMarket::Link*>> PSAndAlphaTuple;
    
        akml::Matrix<float, P_DIMENSION, 1> P;
        SexualMarket *world;
        float utility = 0;
    
        PSAndAlphaTuple buildPSAndAlpha (const akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1>& relations);
        akml::DynamicMatrix<float> computeUtilityGrad(akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1>& relations, PSAndAlphaTuple& PS_Alpha);
        std::tuple<SexualMarket::Link*, Individual*, SexualMarket::Link, bool> preprocessTakeAction(Individual* target=nullptr);
        

	public:
        static inline unsigned short int GREEDY_SHARE = 0;
        
        bool is_greedy;
        float gamma;
        float delta;
        unsigned long int agentid;
    
		Individual(SexualMarket& world, unsigned long int agentid);
        akml::Matrix<float, P_DIMENSION, 1>& getP();
        bool takeAction();
        bool responseToAction(Individual* from, float new_weight);
        akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1> getRelations();
        std::vector<SexualMarket::Link> getScope();
        float computeUtility(akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1>* relations);
};
#endif /* Individual_hpp */
