# PHGS-CVRP: A Parallel Hybrid Genetic Search for the CVRP

This is a parallel implementation of the Hybrid Genetic Search Algorithm for the CVRP [1], [2]. It uses OpenMP to parallelize
the whole algorithm. It uses as many threads as OpenMP offers. Every thread manages its own randomly generated population.
Every process uses different seed for the random number generator, thus increasing randomness in each process' starting
position and direction in the search space.

Some parts of this README are inherited from [2].

## References

[1] Vidal, T. (2022). Hybrid genetic search for the CVRP: Open-source implementation and SWAP* neighborhood. Computers & Operations Research, 140, 105643.
https://doi.org/10.1016/j.cor.2021.105643 (Available [HERE](https://arxiv.org/abs/2012.10384) in technical report form).

[2] https://github.com/vidalt/HGS-CVRP/releases/tag/v1.0.0


## Other programming languages

There exist wrappers for this code in the following languages:
* **C**: The **C_Interface** file contains a simple C API
* **Python**: The [PyHygese](https://github.com/chkwon/PyHygese) package is maintained to interact with the latest release of this algorithm
* **Julia**: The [Hygese.jl](https://github.com/chkwon/Hygese.jl) package is maintained to interact with the latest release of this algorithm

We encourage you to consider using these wrappers in your different projects.
Please contact me if you wish to list other wrappers and interfaces in this section.

## Compiling the executable 

You need [`CMake`](https://cmake.org) to compile.

Build with:
```console
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles"
make bin
```
This will generate the executable file `hgs` in the `build` directory.

## Running the algorithm

After building the executable, try an example: 
```console
./hgs ../Instances/CVRP/X-n157-k13.vrp mySolution.sol -seed 1 -t 30
```

The following options are supported:
```
Call with: ./hgs instancePath solPath [-it nbIter] [-t myCPUtime] [-bks bksPath] [-seed mySeed] [-veh nbVehicles] [-log verbose]
[-it <int>] sets a maximum number of iterations without improvement. Defaults to 20,000                                     
[-t <double>] sets a time limit in seconds. If this parameter is set, the code will be run iteratively until the time limit           
[-seed <int>] [Deprecated, since every process uses its own seed] sets a fixed seed. Defaults to 0                                                                                    
[-veh <int>] sets a prescribed fleet size. Otherwise a reasonable UB on the fleet size is calculated                      
[-round <bool>] rounding the distance to the nearest integer or not. It can be 0 (not rounding) or 1 (rounding). Defaults to 1. 
[-log <bool>] sets the verbose level of the algorithm log. It can be 0 or 1. Defaults to 1.                                       

Additional Arguments:
[-nbIterTraces <int>] Number of iterations between traces display during HGS execution. Defaults to 500
[-nbGranular <int>] Granular search parameter, limits the number of moves in the RI local search. Defaults to 20               
[-mu <int>] Minimum population size. Defaults to 25                                                                            
[-lambda <int>] Number of solutions created before reaching the maximum population size (i.e., generation size). Defaults to 40
[-nbElite <int>] Number of elite individuals. Defaults to 5                                                                    
[-nbClose <int>] Number of closest solutions/individuals considered when calculating diversity contribution. Defaults to 4     
[-nbIterPenaltyManagement <int>] Number of iterations between penalty updates. Defaults to 100
[-targetFeasible <double>] target ratio of feasible individuals between penalty updates. Defaults to 0.2
[-penaltyIncrease <double>] penalty increase if insufficient feasible individuals between penalty updates. Defaults to 1.2
[-penaltyDecrease <double>] penalty decrease if sufficient feasible individuals between penalty updates. Defaults to 0.85
```

There exist different conventions regarding distance calculations in the academic literature.
The default code behavior is to apply integer rounding, as it should be done on the X instances of Uchoa et al. (2017).
To change this behavior (e.g., when testing on the CMT or Golden instances), give a flag `-round 0`, when you run the executable.

The progress of the algorithm in the standard output will be displayed as:

``
It [N1] | T(s) [T] | Feas [NF] [BestF] [AvgF] | Inf [NI] [AvgBestI] [AvgI] | Div [DivF] [DivI] 

``
```
[N1]: Total number of iterations
[T]: Wall time spent until now
[NF] and [NI]: Average number of feasible and infeasible solutions in the subpopulations 
[BestF]: Value of the best feasible solution in the subpopulations 
[AvgBestI]: Average best infeasible solution value
[AvgF] and [AvgI]: Average value of the solutions in the feasible and infeasible subpopulations 
[DivF] and [DivI]: Diversity of the feasible and infeasible subpopulations
```

## Code structure

The main classes containing the logic of the algorithm are the following:
* **Params**: Stores the main data structures for the method
* **Individual**: Represents an individual solution in the genetic algorithm, also provides I/O functions to read and write individual solutions in CVRPLib format.
* **Population**: Stores the solutions of the genetic algorithm into two different groups according to their feasibility. Also includes the functions in charge of diversity management.
* **Genetic**: Contains the main procedures of the genetic algorithm as well as the crossover
* **LocalSearch**: Includes the local search functions, including the SWAP* neighborhood
* **Split**: Algorithms designed to decode solutions represented as giant tours into complete CVRP solutions
* **CircleSector**: Small code used to represent and manage arc sectors (to efficiently restrict the SWAP* neighborhood)

In addition, additional classes have been created to facilitate interfacing:
* **AlgorithmParameters**: Stores the parameters of the algorithm
* **CVRPLIB** Contains the instance data and functions designed to read input data as text files according to the CVRPLIB conventions
* **commandline**: Reads the line of command
* **main**: Main code to start the algorithm
* **C_Interface**: Provides a C interface for the method

## Compiling the shared library

The original code at [2] included a possibility to compile the code as a shared library. This repository does not aim 
to provide the same.

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**
- Copyright(c) 2020 Thibaut Vidal
- Copyright(c) 2023 Roman Reimche
