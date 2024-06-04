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
#include <algorithm>
#include <filesystem>

#ifndef GRAPH_SIZE
    #define GRAPH_SIZE 20
#endif

#ifndef MAX_THREADS_USAGE
    #define MAX_THREADS_USAGE 85
#endif

#define P_DIMENSION 100
#define MODE_ECO_LOG 0
#define MODE_FOLDER_LOG 1

#define LINKS_NB ((GRAPH_SIZE*GRAPH_SIZE-GRAPH_SIZE)/2)

#ifndef Constants_h
#define Constants_h

//#include "SocialMatrix.hpp"
//#include "Individual.hpp"

#endif /* Constants_h */
