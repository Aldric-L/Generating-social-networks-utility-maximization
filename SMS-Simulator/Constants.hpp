//
//  Constants.hpp
//  Social Matrix Simulation
//
//  Created by Aldric Labarthe on 14/11/2023.
//

#include "CMakeConsts.h"

#define SMS_VERSION "0.71 20/08/24"

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
#include <deque>

#ifndef GRAPH_SIZE
    #define GRAPH_SIZE (20)
#endif

#ifndef MAX_THREADS_USAGE
    #define MAX_THREADS_USAGE (85)
#endif

#define MODE_ECO_LOG 0

#define DEPRECIATION_RATE (0)
//#define MIN_LINK_WEIGHT (0.005)
#define MIN_LINK_WEIGHT (-MAXFLOAT)
#define MAX_LINK_CHANGE (0.1)


#define MIN_LINK_WEIGHT2 0.05
#define MIN_GRADIENT_MOVE (0.0001)
#define MAX_GRADIENT_MOVE (0.1)

#define MAX_INACTIVE_ROUNDS (50)


#ifndef Constants_h
#define Constants_h

//#include "SocialMatrix.hpp"
//#include "Individual.hpp"

#endif /* Constants_h */
