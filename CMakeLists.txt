cmake_minimum_required(VERSION 3.27)
project(SocialMatrixSimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(SOURCES SMS-Simulator/Constants.hpp
SMS-Simulator/Individual.hpp
SMS-Simulator/Individual.cpp
SMS-Simulator/main.cpp
SMS-Simulator/UtilityFunction.hpp
SMS-Simulator/UtilityFunction.cpp
SMS-Simulator/OptimalMatrix.cpp
SMS-Simulator/OptimalMatrix.hpp
SMS-Simulator/SocialMatrix.hpp
SMS-Simulator/SocialMatrix.cpp)
add_executable(SocialMatrixSimulation ${SOURCES})

if(DEFINED GRAPH_SIZE)
    message("Setting GRAPH SIZE in CMake to ${GRAPH_SIZE}")
    set(CMAKE_GRAPH_SIZE ${GRAPH_SIZE})
endif(DEFINED GRAPH_SIZE)

if(DEFINED MAX_THREADS_USAGE)
    message("Setting MAX_THREADS_USAGE in CMake to ${MAX_THREADS_USAGE}")
    set(CMAKE_MAX_THREADS_USAGE ${MAX_THREADS_USAGE})
endif(DEFINED MAX_THREADS_USAGE)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("CMake Debug mode enabled")
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)
    
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/SMS-Simulator/CMakeConsts.h.in ${CMAKE_CURRENT_SOURCE_DIR}/SMS-Simulator/CMakeConsts.h @ONLY)

include_directories (SMS-Simulator/)
add_subdirectory(SMS-Simulator/AKML-lib/)
target_link_libraries(SocialMatrixSimulation PRIVATE AKML)
target_include_directories(SocialMatrixSimulation PRIVATE src/AKML-lib/)
set_property(TARGET AKML PROPERTY CXX_STANDARD 20)
