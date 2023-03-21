//
//  SexualMarket.hpp
//  Sexual Market Simulation
//
//  Created by SMS Associates on 21/03/2023.
//
#pragma once

#include <stdio.h>
#include <array>
#include <vector>

#ifndef SexualMarket_hpp
#define SexualMarket_hpp

class Individual;

class SexualMarket {

private:
    std::array<int, 2> dimension;
    unsigned int nb_individuals;
    std::vector<Individual*> individuals;

public:
    // Constructor
    SexualMarket(std::array<int, 2> dimension, unsigned int nb_individuals);

    // Getter
    std::array<int, 2> getDimension();
    std::vector<Individual*> getAliveIndividuals();

    std::vector<Individual*> getVisibleIndividuals(int const vision, std::array<int, 2> pos);



};

 /* SexualMarket_hpp */
#endif