//
//  Individual.hpp
//  Sexual Market Simulation
//
//  Created by Yann Kerzreho on 21/03/2023.
//
#include "SexualMarket.hpp"
#ifndef Individual_hpp
#define Individual_hpp

#define M_DIMENSIONS 3

#include <stdio.h>
#include <array>

class Individual {
	private:
		bool available_on_SMS = true;
		const std::array <int, M_DIMENSIONS> M;
		const std::array <float, M_DIMENSIONS> W;
		std::array <int, 2> coord;
		Individual* target;
		SexualMarket* world;
    
	protected:
		int const vision_radius = 12;
		int const action_radius = vision_radius / 2;

	public:
		// Getters
		bool is_available();
		std::array <int, M_DIMENSIONS> getM();
		std::array <float, M_DIMENSIONS> getW();
		std::array <int, 2> getCoord();
		Individual* getTarget();

};

#endif /* Individual_hpp */
