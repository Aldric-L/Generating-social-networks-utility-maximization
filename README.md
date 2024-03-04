#  Social Matrix Simulation
### SMS Associates - 2023

## Build and installation 

You need to use CMake. We provide you a script to launch the simulator named `customlaunchscript.sh`. You should use as first argument the number of individuals you want to simulate. Example: `./customlaunchscript.sh 500` will compile the project for 500 individuals. The default value is set to 200. 

Under the hood, we use the following instructions, to compile in Debug mode for now.

```
cmake -B build -DCMAKE_BUILD_TYPE=Debug 
cmake --build build --config Debug
```
