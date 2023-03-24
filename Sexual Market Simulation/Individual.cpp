//
//  Individual.cpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//

#include "Individual.hpp"
#include "SexualMarket.hpp"

Individual::Individual(SexualMarket& world, std::array<int, 2> coord){
	this->world = &world;
	this->coord = coord;

	for (int i(0); i < M_DIMENSIONS; i++) {
		this->M[i] = rand();
		this->W[i] = (float)rand();
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

std::vector<Individual*> Individual::getVisibleIndividuals()
{
	std::vector<Individual*> visibles = this->world->getVisibleIndividuals(this->vision_radius, this->getCoord());
	std::vector<Individual*> to_return;
	for (int i(0); i < visibles.size(); i++) {
		if (visibles[i] != this)
			to_return.push_back(visibles[i]);
	}
	return to_return;
}

void Individual::setVisibleIndividuals()
{
    visibleIndividuals = getVisibleIndividuals();
}

float Individual::score(Individual &scoredIndividual)
{
    float score = 0
    for(int i = 0; i < M_DIMENSIONS; i++){
        score += getM()[i] * scoredIndividual.getM()[i];
    }
    return score;
}
