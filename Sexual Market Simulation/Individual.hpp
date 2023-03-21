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
private:
	bool available_on_SMS = true;
	int M[M_DIMENSIONS] const;
	float W[M_DIMENSIONS] const;
	int coord[2];
	Individual* target;
    
public:
};

#endif /* Individual_hpp */
