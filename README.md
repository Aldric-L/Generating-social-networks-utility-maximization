#  Social Matrix Simulation
### SMS Associates Yann Kerzreho & Aldric Labarthe @ ENS Paris-Saclay - 2023-2024

Our model is a bilateral network of agents where each of them play iteratively to optimize their own utility. In this purpose they optimize the weights of their relations within a given scope. Each node/agent has its own carateristics. The model is provided in two forms: a global optimization problem and a round by round simulation with unperfect information.

## Usage

### Build and installation

The simulator is implemented in C++ as a CLI application. It does not rely on any external material (except the library AKML that is statically linked in the final executable). Despite being highly versatile the **simulator is built to generate fixed-sized graphs**. Hence, you need to recompile it each time you want to change the size of the graph. 

To make it simple, we use CMake. The project is developed using CLang and LLVM but is tested with gcc for compatibility. We provide you a script to launch the simulator named `customlaunchscript.sh`. 

**How to use the launch script "customlaunchscript.sh"?**
The script features three main options:

- the graphsize (unsigned int) : the first parameter without any flag (default: 200)

- the max CPU share allowed (unsigned int, (0,100)) : the second parameter without any flag (default: 85)

- the compiler mode (Debug/Release) : use -R for release, debug is default

Example: `./customlaunchscript.sh 500 90 -R` will compile the project for 500 individuals in Release mode with a CPU max usage of 90%. 


If you want to build without the script, you can use:

```
cmake -B build -DCMAKE_BUILD_TYPE=Debug 
cmake --build build --config Debug
```

### CLI Usage
Use the flag -h to ask for help!
```
Usage: ./SocialMatrixSimulation [options]
Options:
  -R --rounds                       How many rounds?
  -S --simuls                       How many simulations?
  -O --optiGraph                    Should I compute optimal graph?
  -op --optiPrecision                Precision of the optimal graph 1e{-x}?
  -P --dimP                         Dimension of the vector of personnality
  -m --memLength                    Length of the memory of past requests
  -s --greedyS                      Share of greedy individuals [0,100]
  -f --greedyF                      Frequency of the greedy bonus [0,100]
  -D --delta                        Utility parameter
  -K --kappa                        Utility parameter
  -G --gamma                        Utility parameter
  -g --gammaDisp                    Dispersion of the utility parameter gamma
  -p --htroP                        Enable/Disable the two groups of P
  -c --clearing                     Enable/Disable the clearing and decaying mechanism
  -C --clustering                   Enable/Disable the computing of clustering coefficients
  -t --checkLTCond                  Should we check the Love Triangle Condition
  -u --utFreq                       Frequency of the utility log
  -l --log                          Should we log results?
  -L --logPrefix                    Select a subfolder for registering logs
  -a --adjacencyFName               Import adjacency matrix from csv file
  -c --compatibilityFName           Import compatibility matrix from csv file
```

