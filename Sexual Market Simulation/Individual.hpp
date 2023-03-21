//
//  Individual.hpp
//  Sexual Market Simulation
//
//  Created by Yann Kerzreho on 21/03/2023.
//

#ifndef Individual_hpp
#define Individual_hpp

#define M_DIMENSIONS 3

#include <stdio.h>

class Individual {
	bool available_on_SMS = true;
	int M[M_DIMENSIONS];
	float W[M_DIMENSIONS];
	int coord[2];
	Individual* target;
};

#endif /* Individual_hpp */
