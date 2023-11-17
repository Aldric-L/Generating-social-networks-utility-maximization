project(SexualMarketSimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

set(SOURCES SMS-Simulator/Constants.hpp
SMS-Simulator/Individual.hpp
SMS-Simulator/Individual.cpp
SMS-Simulator/main.cpp
SMS-Simulator/SexualMarket.hpp
SMS-Simulator/SexualMarket.cpp)
add_executable(SexualMarketSimulation ${SOURCES})

