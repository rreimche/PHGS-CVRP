#include "Genetic.h"
#include <omp.h>

void Genetic::run()
{	
	/* INITIAL POPULATION */
/*#pragma omp parallel for num_threads(nMaxThreads)
    for(int i=0; i < nMaxThreads; i++){
        populations[i].generatePopulation();
        if (paramsPerThread[thread_num].verbose) std::cout << "----- THREAD " << omp_get_thread_num() <<" DONE WITH POPULATION INIT" << std::endl;
    }*/

    int nbIterGlobal = 0;


#pragma omp parallel num_threads(nMaxThreads)
    {

        int nbIterNonProd = 1;
        int thread_num = omp_get_thread_num();

        /* INITIAL POPULATION */
        populations[thread_num].generatePopulation();
        //if (paramsPerThread[thread_num].verbose) std::cout << "----- THREAD " << omp_get_thread_num() <<" DONE WITH POPULATION INIT" << std::endl;

        Population& population = populations[thread_num];
        LocalSearch& localSearch = localSearchPerThread[thread_num];
        Individual& offspring = offsprings[thread_num];
        Split& split = splits[thread_num];

#pragma omp barrier

        //if (paramsPerThread[thread_num].verbose) std::cout << "----- THREAD " << thread_num <<" STARTING GENETIC ALGORITHM" << std::endl;
        if (paramsPerThread[thread_num].verbose) std::printf("----- THREAD %d STARTING GENETIC ALGORITHM\n", thread_num);

#pragma omp barrier
        for (int nbIterThread = 0 ; nbIterNonProd <= paramsPerThread[thread_num].ap.nbIter && (paramsPerThread[thread_num].ap.timeLimit == 0 || (omp_get_wtime()-paramsGlobal.startTime) < paramsGlobal.ap.timeLimit) ; nbIterThread++)
        {
            /* SELECTION AND CROSSOVER */
            crossoverOX(offspring, population.getBinaryTournament(), population.getBinaryTournament(), split);

            //if (paramsPerThread[thread_num].verbose) std::cout << "----- THREAD " << thread_num <<" FINISHED CROSSOVER IN GENERATION " << nbIterThread << std::endl;

            /* LOCAL SEARCH */
            localSearch.run(offspring, paramsPerThread[thread_num].penaltyCapacity, paramsPerThread[thread_num].penaltyDuration);

            //if (paramsPerThread[thread_num].verbose) std::cout << "----- THREAD " << thread_num <<" FINISHED LOCAL SEARCH IN GENERATION " << nbIterThread << std::endl;

            /* POPULATION MANAGEMENT AND REPAIR PROCEDURE */
            bool isNewBest = population.addIndividual(offspring,true);
            if (!offspring.eval.isFeasible && paramsPerThread[thread_num].ran()%2 == 0) // Repair half of the solutions in case of infeasibility
            {
                localSearch.run(offspring, paramsPerThread[thread_num].penaltyCapacity*10., paramsPerThread[thread_num].penaltyDuration*10.);
                if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false) || isNewBest);
            }

            //if (paramsPerThread[thread_num].verbose) std::cout << "----- THREAD " << thread_num <<" FINISHED LOCAL SEARCH IN GENERATION " << nbIterThread << std::endl;

            /* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
            if (isNewBest) nbIterNonProd = 1;
            else nbIterNonProd ++ ;

            /* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
            if (nbIterThread % paramsPerThread[thread_num].ap.nbIterPenaltyManagement == 0) population.managePenalties();
            //if (nbIterThread % paramsPerThread[thread_num].ap.nbIterTraces == 0) population.printState(nbIterThread, nbIterNonProd, thread_num);

            /* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/
            if (paramsPerThread[thread_num].ap.timeLimit != 0 && nbIterNonProd == paramsPerThread[thread_num].ap.nbIter)
            {
                population.restart();
                nbIterNonProd = 1;

                //if (paramsPerThread[thread_num].verbose) std::cout << "----- THREAD " << thread_num <<" HAS RESTARTED THE POPULATION IN GENERATION " << nbIterThread << std::endl;
            }

            if((nbIterThread + 1) % paramsPerThread[thread_num].ap.exchangeRate == 0) {

                #pragma omp barrier
                // At this point, all the threads have finished the iteration

                #pragma omp single
                {
                    bestOfTheBest = getBestOfTheBest();
                }

                // TODO can we have duplicates?
                populations[thread_num].addIndividual(*bestOfTheBest, true);


                #pragma omp barrier
                // At this point, exchange between the threads is finished

            }

            #pragma omp single
            {
                if (paramsGlobal.verbose && (nbIterThread + 1) % paramsGlobal.ap.nbIterTraces == 0) {
                    searchProgress.push_back(std::make_tuple(omp_get_wtime() - paramsGlobal.startTime , nbIterThread + 1, bestOfTheBest->eval.penalizedCost));
                    StateAvg state = getState();
                    printState(nbIterThread, state);
                }

                nbIterGlobal++;
            }

        }
    }

    searchProgress.push_back(std::make_tuple(omp_get_wtime() - paramsGlobal.startTime, nbIterGlobal + 1, bestOfTheBest->eval.penalizedCost));
	if (paramsGlobal.verbose) std::cout << "----- PARALLEL GENETIC ALGORITHM FINISHED. TIME SPENT: " << (omp_get_wtime() - paramsGlobal.startTime) << std::endl;

}

void Genetic::crossoverOX(Individual &result, const Individual &parent1, const Individual &parent2, Split &split)
{

    int thread_num = omp_get_thread_num();

	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (paramsPerThread[thread_num].nbClients + 1, false);

	// Picking the beginning and end of the crossover zone
	std::uniform_int_distribution<> distr(0, paramsPerThread[thread_num].nbClients-1);
	int start = distr(paramsPerThread[thread_num].ran);
	int end = distr(paramsPerThread[thread_num].ran);

	// Avoid that start and end coincide by accident
	while (end == start) end = distr(paramsPerThread[thread_num].ran);

	// Copy from start to end
	int j = start;
	while (j % paramsPerThread[thread_num].nbClients != (end + 1) % paramsPerThread[thread_num].nbClients)
	{
		result.chromT[j % paramsPerThread[thread_num].nbClients] = parent1.chromT[j % paramsPerThread[thread_num].nbClients];
		freqClient[result.chromT[j % paramsPerThread[thread_num].nbClients]] = true;
		j++;
	}

	// Fill the remaining elements in the order given by the second parent
	for (int i = 1; i <= paramsPerThread[thread_num].nbClients; i++)
	{
		int temp = parent2.chromT[(end + i) % paramsPerThread[thread_num].nbClients];
		if (!freqClient[temp])
		{
			result.chromT[j % paramsPerThread[thread_num].nbClients] = temp;
			j++;
		}
	}

	// Complete the individual with the Split algorithm
	split.generalSplit(result, parent1.eval.nbRoutes);
}

// NOTE optimisation potential: all the threads do the same, how can we reduce this? is it worth working on?
Individual* Genetic::getBestOfTheBest(){
    Individual* newBestOfTheBest = nullptr;
    bool updated = false;
    for (int i = 0; i < nMaxThreads; i++) {
        if (populations[i].getBestFound() != nullptr) {
            if (newBestOfTheBest == nullptr || populations[i].getBestFound()->biasedFitness > newBestOfTheBest->biasedFitness) {
                newBestOfTheBest = new Individual(*populations[i].getBestFound());
                updated = true;
            }
        }
    }

    if (!updated) newBestOfTheBest = new Individual(*bestOfTheBest);

    return newBestOfTheBest;
}


StateAvg Genetic::getState()
{
    StateAvg state;

    int sumFeasibleSubpopSize = 0;
    int sumInfeasibleSubpopSize = 0;
    double sumFeasibleBestCosts = 0;
    double sumInfeasibleBestCosts = 0;
    double sumFeasibleAvgCosts = 0;
    double sumInfeasibleAvgCosts = 0;
    double sumFeasibleDiversity = 0;
    double sumInfeasibleDiversity = 0;

    for(int i = 0; i < nMaxThreads; i++){
        sumFeasibleSubpopSize += populations[i].getFeasibleSubpopSize();
        sumInfeasibleSubpopSize += populations[i].getInfeasibleSubpopSize();
        sumFeasibleBestCosts += populations[i].getBestFeasible()->eval.penalizedCost;
        sumInfeasibleBestCosts += populations[i].getBestInfeasible()->eval.penalizedCost;
        sumFeasibleAvgCosts += populations[i].getAverageFeasibleCost();
        sumInfeasibleAvgCosts += populations[i].getAverageInfeasibleCost();
        sumFeasibleDiversity += populations[i].getFeasibleDiversity();
        sumInfeasibleDiversity += populations[i].getInfeasibleDiversity();
    }


    state.avgFeasibleSubpopSize = sumFeasibleSubpopSize / nMaxThreads;
    state.avgInfeasibleSubpopSize = sumInfeasibleSubpopSize / nMaxThreads;
    state.avgFeasibleBestCosts = sumFeasibleBestCosts / nMaxThreads;
    state.avgInfeasibleBestCosts = sumInfeasibleBestCosts / nMaxThreads;
    state.avgFeasibleAvgCosts = sumFeasibleAvgCosts / nMaxThreads;
    state.avgInfeasibleAvgCosts = sumInfeasibleAvgCosts / nMaxThreads;
    state.avgFeasibleDiversity = sumFeasibleDiversity / nMaxThreads;
    state.avgInfeasibleDiversity = sumInfeasibleDiversity / nMaxThreads;

    return state;
}

void Genetic::printState(int nbIter, StateAvg avg) const {
    std::printf("It %6d | T(s) %.2f", nbIter + 1,  omp_get_wtime()-paramsGlobal.startTime);

    if(avg.avgFeasibleSubpopSize != 0) std::printf(" | Feas %.2f %.2f %.2f", avg.avgFeasibleSubpopSize, bestOfTheBest->eval.penalizedCost, avg.avgFeasibleAvgCosts);
    else std::printf(" | NO-FEASIBLE");

    if (avg.avgInfeasibleSubpopSize != 0) std::printf(" | Inf %.2f %.2f %.2f", avg.avgInfeasibleSubpopSize, avg.avgInfeasibleBestCosts, avg.avgInfeasibleAvgCosts);
    else std::printf(" | NO-INFEASIBLE");

    std::printf(" | Div %.2f %.2f", avg.avgFeasibleDiversity, avg.avgInfeasibleDiversity);

    std::cout << std::endl;
}

void Genetic::exportSearchProgress(std::string fileName, std::string instanceName)
{
    std::ofstream myfile(fileName);
    for (std::tuple<double, int, double> state : searchProgress)
    myfile << instanceName << ";" << paramsGlobal.ap.seed << ";" << std::get<2>(state) << ";" << std::get<1>(state) << ";" << std::get<0>(state) << std::endl;
}

//TODO add parameter for number of threads/demes
Genetic::Genetic(Params & params) : 
    //nMaxThreads(omp_get_max_threads()),
    nMaxThreads(omp_get_max_threads()),
    bestOfTheBest(nullptr),
    paramsPerThread(std::vector<Params>(nMaxThreads, params)),
    paramsGlobal(params)
	{

        // Avoid dangling references by reserving the memory in advance

        splits.reserve(nMaxThreads);
        populations.reserve(nMaxThreads);

        // Populate vectors without dependencies on each other
        for(int i = 0; i < nMaxThreads; i++) {
            paramsPerThread[i].ran = std::minstd_rand(i);
            splits.emplace_back(paramsPerThread[i]);
            localSearchPerThread.emplace_back(paramsPerThread[i]);
            offsprings.emplace_back(paramsPerThread[i]);
        }

        // Populate vectors with dependencies on prior results
        for (int i = 0; i < nMaxThreads; i++) {
            populations.emplace_back(paramsPerThread[i], splits[i], localSearchPerThread[i]);
        }
}



