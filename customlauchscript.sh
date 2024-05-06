#!/bin/bash -e

GH="${1:-200}"
BUILD_TYPE="Debug"


while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -R)
        BUILD_TYPE="Release"
        shift # past argument
        ;;
        *)    # unknown option
        ;;
    esac
    shift # past argument or value
done

echo "GRAPH_SIZE is set to: $GH"
echo "Build type is set to: $BUILD_TYPE"
cmake -B build -D CMAKE_BUILD_TYPE=$BUILD_TYPE -D GRAPH_SIZE=$GH || exit 1
cmake --build build --config $BUILD_TYPE || exit 1
./build/SocialMatrixSimulation
