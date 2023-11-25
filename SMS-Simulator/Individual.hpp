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
    using PSAndAlphaType = std::tuple<akml::Matrix<float, P_DIMENSION, GRAPH_SIZE-1>, akml::Matrix<float, GRAPH_SIZE-1, 1>, akml::Matrix<Individual*, GRAPH_SIZE-1, 1>, akml::Matrix<SexualMarket::Link*, GRAPH_SIZE-1, 1>>;
    
	protected:
        akml::Matrix<float, P_DIMENSION, 1> P;
        //std::array <bool, P_DIMENSION> P;
        SexualMarket *world;
        float utility = 0;
        float gamma;
    
        PSAndAlphaType buildPSAndAlpha (const std::array<SexualMarket::Link*, GRAPH_SIZE-1>& relations);
        akml::Matrix<float, GRAPH_SIZE-1, 1> computeUtilityGrad(std::array<SexualMarket::Link*, GRAPH_SIZE-1>* relations=nullptr, PSAndAlphaType* PS_Alpha=nullptr);
        

	public:
		Individual(SexualMarket& world);
        akml::Matrix<float, P_DIMENSION, 1> getP();
        void takeAction();
        void responseToAction(Individual* from, float new_weight);
        std::array<SexualMarket::Link*, GRAPH_SIZE-1> getRelations();
        std::vector<SexualMarket::Link> getScope();
        float computeUtility(std::array<SexualMarket::Link*, GRAPH_SIZE-1>* relations);
};
#endif /* Individual_hpp */
