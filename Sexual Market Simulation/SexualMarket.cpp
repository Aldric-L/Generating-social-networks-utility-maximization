//
//  SexualMarket.cpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//

#include "SexualMarket.hpp"
#include "Individual.hpp"

SexualMarket::SexualMarket(std::array<int, 2> dimension, unsigned int nb_individuals) :
	dimension(dimension), nb_individuals(nb_individuals) {
	
	std::vector<std::array<int, 2>> used_coord;

	for (int i(0); i < nb_individuals; i++) {
		std::array<int, 2> cord = { 1 + rand()%255, 1 + rand()% 255 };

		while (
			std::find(std::begin(used_coord), std::end(used_coord), cord) != std::end(used_coord) || cord[0] > dimension[0] || cord[1] > dimension[1]) {
			cord = { 1 + rand() % 255, 1 + rand() % 255 };
		}
		used_coord.push_back(cord);
		individuals.push_back(new Individual(*this, cord));
	}

}


std::array<int, 2> SexualMarket::getDimension()
{
	return SexualMarket::dimension;
}

std::vector<Individual*> SexualMarket::getAliveIndividuals()
{
	return SexualMarket::individuals;
}

// A continuer de terminer
std::vector<Individual*> SexualMarket::getVisibleIndividuals(int const vision_radius, std::array<int, 2> pos)
{	
	std::vector<Individual*> visible_individuals;

	std::array<std::array<int, 2>, 2> maxcoord;
	maxcoord[0] = { pos[0] - vision_radius, pos[0] + vision_radius };
	maxcoord[1] = { pos[1] - vision_radius, pos[1] + vision_radius };

	if (maxcoord[0][0] < 0)
		maxcoord[0][0] = 0;
	if (maxcoord[1][0] < 0)
		maxcoord[1][0] = 0;
	if (maxcoord[0][1] > this->getDimension()[0])
		maxcoord[0][1] = this->getDimension()[0];
	if (maxcoord[1][1] > this->getDimension()[1])
		maxcoord[1][1] = this->getDimension()[1];

	std::cout << "Pos: " << pos[0] << ":" << pos[1] << std::endl;
	std::cout << "Max " << maxcoord[0][0] << ":" << maxcoord[0][1] << " - " << maxcoord[1][0] << ":" << maxcoord[1][1] << std::endl;

	for (int i(0); i < SexualMarket::individuals.size(); i++) {
		std::array<int, 2> coord = individuals[i]->getCoord();

		if (maxcoord[0][0] <= coord[0] &&
			coord[0] <= maxcoord[0][1] && maxcoord[1][0] <= coord[1] 
			&& coord[1] <= maxcoord[1][1])
			visible_individuals.push_back(individuals[i]);

		return visible_individuals;
	}

}
