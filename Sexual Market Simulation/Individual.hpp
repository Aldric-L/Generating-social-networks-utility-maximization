//
//  Individual.hpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//
#pragma once

#include <stdio.h>
#include <iostream>
#include <array>
#include <vector>

#ifndef Individual_hpp
#define Individual_hpp

#define M_DIMENSIONS 3
#define BREAKAVERSION 30


class SexualMarket;

class Individual {
	private:
		bool available_on_SMS = true;
		std::array <int, M_DIMENSIONS> M;
		std::array <float, M_DIMENSIONS> W;
		std::array <int, 2> coord;
		int delta_breaklost = BREAKAVERSION;
		Individual* target = nullptr;
		SexualMarket* world = nullptr;
        std::vector <Individual*> visibleIndividuals;
    
	protected:
		int const vision_radius = 12;
		int const action_radius = vision_radius / 2;

	public:
		// Getters
		Individual(SexualMarket& world, std::array<int, 2> coord);
		bool is_available();
		std::array <int, M_DIMENSIONS> getM();
		std::array <float, M_DIMENSIONS> getW();
		std::array <int, 2> getCoord();
		Individual* getTarget();
		std::vector <Individual*> getVisibleIndividuals();

        // Setters
        void setVisibleIndividuals();
    
        //Other
        float score(Individual& scoredIndividual);
};

#endif /* Individual_hpp */
