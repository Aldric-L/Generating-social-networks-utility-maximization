#!/bin/bash -e

GH="${1:-200}"
echo "GRAPH_SIZE is set to: $GH"
cmake -B build -D CMAKE_BUILD_TYPE=Debug -D GRAPH_SIZE=$GH || exit 1
cmake --build build --config Debug || exit 1
./build/SocialMatrixSimulation
