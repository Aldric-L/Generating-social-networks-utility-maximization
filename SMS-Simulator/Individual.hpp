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
        typedef std::tuple<akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1>, akml::Matrix<float, GRAPH_SIZE-1, 1>, akml::Matrix<Individual*, GRAPH_SIZE-1, 1>, akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1>> PSAndAlphaTuple;
    
        akml::Matrix<float, P_DIMENSION, 1> P;
        SexualMarket *world;
        float utility = 0;
    
        PSAndAlphaTuple buildPSAndAlpha (const std::array<SexualMarket::Link*, GRAPH_SIZE-1>& relations);
        akml::Matrix<float, GRAPH_SIZE-1, 1> computeUtilityGrad(std::array<SexualMarket::Link*, GRAPH_SIZE-1>* relations=nullptr, PSAndAlphaTuple* PS_Alpha=nullptr);
        std::tuple<SexualMarket::Link*, Individual*, SexualMarket::Link, bool> preprocessTakeAction(Individual* target=nullptr);
        

	public:
        float gamma;
        float delta;
        int agentid;
		Individual(SexualMarket& world, int agentid);
        akml::Matrix<float, P_DIMENSION, 1> getP();
        void takeAction();
        void responseToAction(Individual* from, float new_weight);
        std::array<SexualMarket::Link*, GRAPH_SIZE-1> getRelations();
        std::vector<SexualMarket::Link> getScope();
        float computeUtility(std::array<SexualMarket::Link*, GRAPH_SIZE-1>* relations);
};
#endif /* Individual_hpp */
