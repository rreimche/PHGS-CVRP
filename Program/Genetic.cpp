#include "Genetic.h"
#include <omp.h>

void Genetic::run()
{	
	/* INITIAL POPULATION */
#pragma omp parallel for num_threads(nMaxThreads)
    for(int i=0; i < nMaxThreads; i++){
        populations[i].generatePopulation();
        if (params.verbose) std::cout << "----- THREAD " << omp_get_thread_num() <<" DONE WITH POPULATION INIT" << std::endl;
    }


#pragma omp parallel num_threads(nMaxThreads)
    {
        int nbIterNonProd = 1;
        int thread_num = omp_get_thread_num();

        Population population = populations[thread_num];
        LocalSearch localSearch = localSearchPerThread[thread_num];
        Individual offspring = offsprings[thread_num];
        Split split = splits[thread_num];

        if (params.verbose) std::cout << "----- THREAD " << thread_num <<" STARTING GENETIC ALGORITHM" << std::endl;

        for (int nbIterThread = 0 ; nbIterNonProd <= params.ap.nbIter && (params.ap.timeLimit == 0 || (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC < params.ap.timeLimit) ; nbIterThread++)
        {
            /* SELECTION AND CROSSOVER */
            crossoverOX(offspring, population.getBinaryTournament(), population.getBinaryTournament(), split);

            /* LOCAL SEARCH */
            localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
            bool isNewBest = population.addIndividual(offspring,true);
            if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
            {
                localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
                if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false) || isNewBest);
            }

            /* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
            if (isNewBest) nbIterNonProd = 1;
            else nbIterNonProd ++ ;

            /* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
            if (nbIterThread % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
            if (nbIterThread % params.ap.nbIterTraces == 0) population.printState(nbIterThread, nbIterNonProd, thread_num);

            /* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/
            if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
            {
                population.restart();
                nbIterNonProd = 1;
            }

            if((nbIterThread + 1) % exchangeRate == 0){

                #pragma omp barrier
                // At this point, all the threads have finished the iteration

                //TODO exchange individuals
                bestOfTheBest = getBestOfTheBest();
                std::cout << "THREAD " << thread_num << " AT BARRIER " << std::endl;

                #pragma omp barrier
                // At this point, exchange between the threads is finished

            }

        }
    }

	if (params.verbose) std::cout << "----- PARALLEL GENETIC ALGORITHM FINISHED. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;
}

void Genetic::crossoverOX(Individual &result, const Individual &parent1, const Individual &parent2, Split &split)
{
	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

	// Picking the beginning and end of the crossover zone
	std::uniform_int_distribution<> distr(0, params.nbClients-1);
	int start = distr(params.ran);
	int end = distr(params.ran);

	// Avoid that start and end coincide by accident
	while (end == start) end = distr(params.ran);

	// Copy from start to end
	int j = start;
	while (j % params.nbClients != (end + 1) % params.nbClients)
	{
		result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
		freqClient[result.chromT[j % params.nbClients]] = true;
		j++;
	}

	// Fill the remaining elements in the order given by the second parent
	for (int i = 1; i <= params.nbClients; i++)
	{
		int temp = parent2.chromT[(end + i) % params.nbClients];
		if (freqClient[temp] == false)
		{
			result.chromT[j % params.nbClients] = temp;
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

Genetic::Genetic(Params & params) : 
	params(params), 
    nMaxThreads(omp_get_max_threads()),
    exchangeRate(params.ap.exchangeRate),
    bestOfTheBest(nullptr)
	{
        // Initialize local search instances per thread
        for (int i = 0; i < nMaxThreads; i++) {
            splits.emplace_back(params);
            localSearchPerThread.emplace_back(params);
            populations.emplace_back(params, splits[i], localSearchPerThread[i]);
            offsprings.emplace_back(params);
        }
}

