//
//  SexualMarket.hpp
//  Sexual Market Simulation
//
//  Created by Yann Kerzreho on 21/03/2023.
//

#ifndef SexualMarket_hpp
#define SexualMarket_hpp

#include <stdio.h>

class SexualMarket{

    
public:
    // Constructor
    SexualMarket(int[2] dimension) : dimension(dimension);
    
    // Getter
    int getDimension() const { return dimension; }

private:
    int[2] const dimension;
}

#endif /* SexualMarket_hpp */
