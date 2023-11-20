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
        akml::Matrix<float, P_DIMENSION, 1> P;
        //std::array <bool, P_DIMENSION> P;
        SexualMarket *world;
        float utility = 0;
        float gamma;

	public:
		Individual(SexualMarket& world);
        akml::Matrix<float, P_DIMENSION, 1> getP();
        void takeAction();
        bool responseToAction();
        std::array<SexualMarket::Link*, GRAPH_SIZE-1> getRelations();
        std::vector<SexualMarket::Link> getScope();
        float computeUtility();
		
};
#endif /* Individual_hpp */
