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

add_subdirectory(SMS-Simulator/AKML-lib/)
target_link_libraries(SexualMarketSimulation PRIVATE AKML)
target_include_directories(SexualMarketSimulation PRIVATE src/AKML-lib/)
set_property(TARGET AKML PROPERTY CXX_STANDARD 20)
