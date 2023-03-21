//
//  Individual.cpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//

#include "Individual.hpp"
#include "SexualMarket.hpp"

Individual::Individual(SexualMarket& world, std::array<int, 2> coord){
	Individual::world = &world;
	Individual::coord = coord;

	for (int i(0); i < M_DIMENSIONS; i++) {
		Individual::M[i] = rand();
		Individual::W[i] = (float)rand();
	}
}

bool Individual::is_available()
{
	return Individual::available_on_SMS;
}

std::array<int, M_DIMENSIONS> Individual::getM()
{
	return Individual::M;
}

std::array<float, M_DIMENSIONS> Individual::getW()
{
	return Individual::W;
}

std::array<int, 2> Individual::getCoord()
{
	return Individual::coord;
}

Individual* Individual::getTarget()
{
	return Individual::target;
}

