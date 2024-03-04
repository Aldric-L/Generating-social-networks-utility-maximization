cmake_minimum_required(VERSION 3.27)
project(SocialMatrixSimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

set(SOURCES SMS-Simulator/Constants.hpp
SMS-Simulator/Individual.hpp
SMS-Simulator/Individual.cpp
SMS-Simulator/main.cpp
SMS-Simulator/SocialMatrix.hpp
SMS-Simulator/SocialMatrix.cpp)
add_executable(SocialMatrixSimulation ${SOURCES})

if(DEFINED GRAPH_SIZE)
    message("Setting GRAPH SIZE in CMake to ${GRAPH_SIZE}")
    set(CMAKE_GRAPH_SIZE ${GRAPH_SIZE})
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/SMS-Simulator/CMakeConsts.h.in ${CMAKE_CURRENT_SOURCE_DIR}/SMS-Simulator/CMakeConsts.h @ONLY)
endif(DEFINED GRAPH_SIZE)


include_directories (SMS-Simulator/)
add_subdirectory(SMS-Simulator/AKML-lib/)
target_link_libraries(SocialMatrixSimulation PRIVATE AKML)
target_include_directories(SocialMatrixSimulation PRIVATE src/AKML-lib/)
set_property(TARGET AKML PROPERTY CXX_STANDARD 20)
