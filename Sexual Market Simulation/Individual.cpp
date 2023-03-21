//
//  Individual.cpp
//  Sexual Market Simulation
//
//  Created by SMS and Associates on 21/03/2023.
//

#include "Individual.hpp"

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

