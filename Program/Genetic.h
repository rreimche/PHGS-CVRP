/*MIT License

Copyright(c) 2020 Thibaut Vidal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#ifndef GENETIC_H
#define GENETIC_H

#include "Population.h"
#include "Individual.h"

// Averages of state/performance data over all populations at a point in time
struct StateAvg {
    double avgFeasibleSubpopSize = 0;
    double avgInfeasibleSubpopSize = 0;
    double avgFeasibleBestCosts = 0;
    double avgInfeasibleBestCosts = 0;
    double avgFeasibleAvgCosts = 0;
    double avgInfeasibleAvgCosts = 0;
    double avgFeasibleDiversity = 0;
    double avgInfeasibleDiversity = 0;
};


class Genetic
{
public:
    const int nMaxThreads;                // Max number of threads
    Params & paramsGlobal;                   // global copy of problem parameters
	std::vector<Params> paramsPerThread;	// Problem parameters
	//Split split;					// Split algorithm
    std::vector<Split> splits;       // Split algorithms one per thread
	//LocalSearch localSearch;		// Local Search structure
    std::vector<LocalSearch> localSearchPerThread; // Local Search instances per thread
	//Population population;			// Population (public for now to give access to the solutions, but should be be improved later on)
    std::vector<Population> populations; // Populations, one per thread
	//Individual offspring;			// First individual to be used as input for the crossover
    std::vector<Individual> offsprings;
    // Keeps tracks of the time stamps of successive best solutions: {time, generation, value}
    std::vector<std::tuple<double, int, double>> searchProgress;

    Individual* bestOfTheBest;            // Best Solution

	// OX Crossover
    void crossoverOX(Individual &result, const Individual &parent1, const Individual &parent2, Split &split);

    // Running the genetic algorithm until maxIterNonProd consecutive iterations or a time limit
    void run() ;

    // Get best individual of all populations
    Individual* getBestOfTheBest();

	// Constructor
	Genetic(Params & params);

    StateAvg getState();
    void printState(int nbIter, StateAvg avg) const;

    void exportSearchProgress(std::string fileName, std::string instanceName);

};

#endif
