//
//  Constants.hpp
//  Social Matrix Simulation
//
//  Created by Aldric Labarthe on 14/11/2023.
//

#include "CMakeConsts.h"

#include <stdio.h>
#include <iostream>
#include <array>
#include <vector>
#include <utility>
#include <random>
#include <cassert>
#include <ctime>
#include <chrono>
#include <memory>
#include <thread>
#include <random>

#ifndef GRAPH_SIZE
    #define GRAPH_SIZE 20
#endif

#define P_DIMENSION 100
#define COMPUTE_CLUSTERING 1
#define COMPUTE_CLEARING 1
#define MODE_ECO_LOG 0
#define MODE_FOLDER_LOG 1

#if MODE_FOLDER_LOG
    #include <filesystem>
#endif

#define LINKS_NB (GRAPH_SIZE*GRAPH_SIZE-GRAPH_SIZE)/2

#ifndef Constants_h
#define Constants_h

//#include "SocialMatrix.hpp"
//#include "Individual.hpp"

#endif /* Constants_h */
