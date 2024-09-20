# SOCIAL NETWORKS AS A UTILITY MAXIMIZATION DECISION
### SMS Associates Yann Kerzreho & Aldric Labarthe @ ENS Paris-Saclay - 2023-2024 - Social Matrix Simulation

Our model is a bilateral network of agents where each of them play iteratively to optimize their own utility. In this purpose they optimize the weights of their relations within a given scope. Each node/agent has its own carateristics. The model is provided in two forms: a global optimization problem and a round by round simulation with imperfect information.

## Paper and reproducibility
All data used in the paper and all the statistical analysis performed is available in the folder Paper. The file PaperMain.R can regerate all tables, where all the other files are R functions for reproducibility. The data is available in two compressed folder at this link ([https://1drv.ms/f/s!An5zxDZ6MkIwo4JOOjl38hN-FeKA-A?e=pnTXqT](https://1drv.ms/f/s!An5zxDZ6MkIwo4JOOjl38hN-FeKA-A?e=pnTXqT)) with both the simulations outputs and inputs. Beware that paths may not be properly configured as these documents are not meant to be executed but are provided for full transparency. 

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

### How to use the Simulator

#### General principle
The Simulator is designed to conduct **round by round** simulations. Optimal graphs (static optimization) are possible but you should keep in mind that they are, at begining, only a side tool of the round by round simulation. 

Every simulation is **single-threaded and full CPU**. Expect the CPU usage of the core to be 100%. Beware that the max CPU share allowed is used to set an upper limit for the number of simulations that can be conducted in parallel (on separated threads). The failure of one simulation do not affect the behavior of the others of the same batch. The system treats simulations by batches and waits for the end of the batch to jump to the next batch.

#### CLI Usage
Use the flag -h to ask for help!
```
Usage: ./SocialMatrixSimulation [options]
Options:
  -R --rounds                                How many rounds?
  -S --simuls                                How many simulations?
  -O --optiGraph                             Should I compute optimal graph?
  -op --optiPrecision                        Precision of the optimal graph 1e{-x}?
  -oi --optiIt                               Max. Iterations for finding the optimal graph
  -P --dimP                                  Dimension of the vector of personnality
  -m --memLength                             Length of the memory of past requests
  -s --greedyS                               Share of greedy individuals [0,100]
  -f --greedyF                               Frequency of the greedy bonus [0,100]
  -d --scopeDepth                            Individuals' scope depth in the network
  -D --delta                                 Utility parameter
  -K --kappa                                 Utility parameter
  -G --gamma                                 Utility parameter
  -g --gammaDisp                             Dispersion of the utility parameter gamma
  -p --htroP                                 Enable/Disable the two groups of P
  -c --clearing                              Enable/Disable the clearing and decaying mechanism
  -C --clustering                            Enable/Disable the computing of clustering coefficients
  -t --checkLTCond                           Should we check the Love Triangle Condition
  -u --utFreq                                Frequency of the utility log
  -l --log                                   Select the log policy (NONE, NO_EDGES, ALL)
  -L --logPrefix                             Select a subfolder for registering logs
  -lF --logSubFolder                         Create a subfolder for each simulation
  -i --initDensity                           Initialization density factor [1,64]
  -A --adjacencyFName                        Import adjacency matrix from csv file
  -Co --compatibilityFName                   Import compatibility matrix from csv file
```

Note:
- If you only want to compute the optimal matrix use `-O=True` and `-R=0`

- `memLength`, `greedyS`, `greedyF`, `clearing`, `utFreq` are only useful in round by round simulations

- `optiPrecision` cannot go beyond 1e-6 in almost every modern system due to the use of floats

- There is no gatekeeper for the CLI instructions, be careful with domains according to the paper.

#### Outputs

At each begining of a simulation, a specific folder is created with the unique-id of the simulation following the pattern `sim_<unique-id>`. In this folder several log files and output files can be found depending on the CLI arguments:

- **Edges log file** (csv): `SMS-Save-Edges-<unique-id>.csv` for round by round simulations. It contains all demands asked by agents. Columns are `"round", "vertex1", "vertex2", "old_weight", "new_weight", "accepted", "forced"`.

- **Utility log file** (csv): `SMS-Save-Utility-<unique-id>.csv` for round by round simulations. It contains the utilities of each agent evaluated at a frequency defined in the CLI. Columns are `"round", "agentid", "utility"`.

- **Vertices log file** (csv): `SMS-Save-Vertices-<unique-id>.csv`. It describes all agents in the simulation (with the value of their individual parameters). Columns are `"round", "agentid", "gamma", "isgreedy", "meandist", "vardist", "maxdist", "P"`. Note that `"meandist", "vardist", "maxdist"` are not computed for pure optimal graph simulations (-1 is written as undefined).

- **Clustering log file** (csv): `SMS-Save-Clustering-<unique-id>.csv` for round by round simulations. It computes the local clustering coefficient of each agent. Each row is the vector of all clustering coefficients for a given round. Beware that the row can contain another undefined value (abandonned global clustering coefficient).

- **AdjacencyMatrix log file** (csv): `SMS-Save-AdjacencyMatrix-<unique-id>.csv` for round by round simulations. Every 250 rounds, the adjacency matrix is logged. Matrices are concatenated by rows. Please note that this log is fully redundant and only here to make data processing easier (the edges log file is often too heavy to be easily processed by R). 

- **OptimalGraph log file** (csv): `SMS-Save-OptimalGraph-<unique-id>.csv` for optimal graph computation. The adjacency matrix is logged, followed by the compatibility matrix. Matrices are concatenated by rows. Please note that the ADAM optimization process stops when it reaches 1 million epochs, whether it has reached convergence or not. 

- **OptimalUtility log file** (csv): `SMS-Save-OptimalUtility-<unique-id>.csv` for optimal graph computation. Outputs in a single column the individual utility that agents are experiencing at optimum. 

- **SimulInfos** (txt): `SMS-SimulInfos-<unique-id>.txt` this files contains all CLI parameters and keeps track of the simulation. Errors and enlapsed time are to be found here.

## Contact
Please contact Aldric L. if you find any bug in this software by the GitHub issues. Remember that this work is licensed under the GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007.

