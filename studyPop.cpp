#include "studyPop.h"



// Constructor: When an object of type "studyPop" is created, it automatically creates its first generation
studyPop::studyPop(int currentRound) {

	if(PARASITES_EXIST) {

		/*** Initialize parasiteInfectionTables ***/

		// perform this part if parasite host model being used is GFG
		if(strcmp("GFG", MODEL_OF_PARASITE_INFECTION) == 0) {
			if(PARASITE_PLOIDY == 1) {
                // Fitness tables based on Agrawal 2009
				for(int i=0; i<9 ; i++) {
					hostFitnessTable[i][2] = 0;
					parasiteFitnessTable[i][2] = 0;
					hostFitnessTable[i][5] = 0;
					parasiteFitnessTable[i][5] = 0;
					hostFitnessTable[i][6] = 0;
					parasiteFitnessTable[i][6] = 0;
					hostFitnessTable[i][7] = 0;
					parasiteFitnessTable[i][7] = 0;
					hostFitnessTable[i][8] = 0;
					parasiteFitnessTable[i][8] = 0;
				}

				// Host Fitness based on Agrawal 2009
				hostFitnessTable[0][0]=hostFitnessTable[0][1]=hostFitnessTable[0][3] = pow((1-PLEIOTROPIC_COST_FOR_HOST), 2.0);
				hostFitnessTable[0][4] = pow((1-PLEIOTROPIC_COST_FOR_HOST), 2.0) * (1-V);
				hostFitnessTable[1][0]=hostFitnessTable[1][1]=hostFitnessTable[1][3] = (1-PLEIOTROPIC_COST_FOR_HOST)* (1-PLEIOTROPIC_COST_FOR_HOST*D1);
				hostFitnessTable[1][4] = (1-PLEIOTROPIC_COST_FOR_HOST)* (1-PLEIOTROPIC_COST_FOR_HOST*D1) * (1-V);
				hostFitnessTable[2][0]=hostFitnessTable[2][1] = (1-PLEIOTROPIC_COST_FOR_HOST);
				hostFitnessTable[2][3]=hostFitnessTable[2][4] = (1-PLEIOTROPIC_COST_FOR_HOST)* (1-V);
				hostFitnessTable[3][0]=hostFitnessTable[3][1]=hostFitnessTable[3][3] = (1-PLEIOTROPIC_COST_FOR_HOST)* (1-PLEIOTROPIC_COST_FOR_HOST*D1);
				hostFitnessTable[3][4] = (1-PLEIOTROPIC_COST_FOR_HOST)* (1-PLEIOTROPIC_COST_FOR_HOST*D1) * (1-V);
				hostFitnessTable[4][0]=hostFitnessTable[4][1]=hostFitnessTable[4][3] = pow((1-PLEIOTROPIC_COST_FOR_HOST*D1), 2.0);
				hostFitnessTable[4][4] = pow((1-PLEIOTROPIC_COST_FOR_HOST*D1), 2.0) * (1-V);
				hostFitnessTable[5][0]=hostFitnessTable[5][1] = (1-PLEIOTROPIC_COST_FOR_HOST*D1);
				hostFitnessTable[5][3]=hostFitnessTable[5][4] = (1-PLEIOTROPIC_COST_FOR_HOST*D1)* (1-V);
				hostFitnessTable[6][0]=hostFitnessTable[6][3] = (1-PLEIOTROPIC_COST_FOR_HOST);
				hostFitnessTable[6][1]=hostFitnessTable[6][4] = (1-PLEIOTROPIC_COST_FOR_HOST)* (1-V);
				hostFitnessTable[7][0]=hostFitnessTable[7][3] = (1-PLEIOTROPIC_COST_FOR_HOST*D1);
				hostFitnessTable[7][1]=hostFitnessTable[7][4] = (1-PLEIOTROPIC_COST_FOR_HOST*D1)* (1-V);
				hostFitnessTable[8][0]=hostFitnessTable[8][1]=hostFitnessTable[8][3]=hostFitnessTable[8][4] = 1-V;

				// Parasite Fitness based on Agrawal 2009
				for(int i=0; i<8 ; i++)
					parasiteFitnessTable[i][0] = 1-T_COST;
				parasiteFitnessTable[8][0] = 1;
				for(int i=0; i<6 ; i++)
					parasiteFitnessTable[i][1] = (1-PLEIOTROPIC_COST_FOR_PARASITE)*(1-T_COST);
				parasiteFitnessTable[6][1] = (1-PLEIOTROPIC_COST_FOR_PARASITE);
				parasiteFitnessTable[7][1] = (1-PLEIOTROPIC_COST_FOR_PARASITE);
				parasiteFitnessTable[8][1] = (1-PLEIOTROPIC_COST_FOR_PARASITE);
				for(int i=0; i<8 ; i++) {
					if(i != 2 && i != 5) {
						parasiteFitnessTable[i][3] = (1-PLEIOTROPIC_COST_FOR_PARASITE)*(1-T_COST);
					}
				}
				parasiteFitnessTable[2][3] = (1-PLEIOTROPIC_COST_FOR_PARASITE);
				parasiteFitnessTable[5][3] = (1-PLEIOTROPIC_COST_FOR_PARASITE);
				parasiteFitnessTable[8][3] = (1-PLEIOTROPIC_COST_FOR_PARASITE);
				for(int i=0; i<9 ; i++)
					parasiteFitnessTable[i][4] = pow((1-PLEIOTROPIC_COST_FOR_PARASITE), 2.0);
			}
		}

        // perform this part if parasite host model being used is MA
		else if(strcmp("MA", MODEL_OF_PARASITE_INFECTION) == 0) {
			if(PARASITE_PLOIDY == 1) {

				for(int i=0; i<9 ; i++) {
					hostFitnessTable[i][2] = 0;
					parasiteFitnessTable[i][2] = 0;
					hostFitnessTable[i][5] = 0;
					parasiteFitnessTable[i][5] = 0;
					hostFitnessTable[i][6] = 0;
					parasiteFitnessTable[i][6] = 0;
					hostFitnessTable[i][7] = 0;
					parasiteFitnessTable[i][7] = 0;
					hostFitnessTable[i][8] = 0;
					parasiteFitnessTable[i][8] = 0;
				}

				// find values for host & parasite fitness tables based on Agrawal 2009
				for(int i=0; i<9 ; i++) {
					for(int j=0; j<5 ; j++) {
						hostFitnessTable[i][j] = 1;
						parasiteFitnessTable[i][j] = 1-T_COST;
					}
				}
				hostFitnessTable[0][0]=hostFitnessTable[2][1]=hostFitnessTable[6][3]=hostFitnessTable[8][4] = 1-V;
				hostFitnessTable[1][0]=hostFitnessTable[3][0]=hostFitnessTable[1][1]=hostFitnessTable[5][1]=hostFitnessTable[3][3]=hostFitnessTable[7][3]=hostFitnessTable[5][4]=hostFitnessTable[7][4] = 1 - (D1 + (1-D1)*D2)*V;
				hostFitnessTable[4][0]=hostFitnessTable[4][1]=hostFitnessTable[4][3]=hostFitnessTable[4][4]= 1 - pow(D1 + (1-D1)*D2, 2.0)*V;
				parasiteFitnessTable[0][0]=parasiteFitnessTable[2][1]=parasiteFitnessTable[6][3]=parasiteFitnessTable[8][4] = 1;
				parasiteFitnessTable[1][0]=parasiteFitnessTable[3][0]=parasiteFitnessTable[1][1]=parasiteFitnessTable[5][1]=parasiteFitnessTable[3][3]=parasiteFitnessTable[7][3]=parasiteFitnessTable[5][4]=parasiteFitnessTable[7][4] = 1 - (1 - D1+(1-D1)*D2)*T_COST;
				parasiteFitnessTable[4][0]=parasiteFitnessTable[4][1]=parasiteFitnessTable[4][3]=parasiteFitnessTable[4][4] = 1 - (1 - pow(D1 + (1-D1)*D2, 2.0))*T_COST;
			}
		}

        // perform this part if parasite host model being used is IMA
		else if(strcmp("IMA", MODEL_OF_PARASITE_INFECTION) == 0) {
			if(PARASITE_PLOIDY == 1) {

				for(int i=0; i<9 ; i++) {
					hostFitnessTable[i][2] = 0;
					parasiteFitnessTable[i][2] = 0;
					hostFitnessTable[i][5] = 0;
					parasiteFitnessTable[i][5] = 0;
					hostFitnessTable[i][6] = 0;
					parasiteFitnessTable[i][6] = 0;
					hostFitnessTable[i][7] = 0;
					parasiteFitnessTable[i][7] = 0;
					hostFitnessTable[i][8] = 0;
					parasiteFitnessTable[i][8] = 0;
				}

				// find values for host & parasite fitness tables based on Agrawal 2009
				for(int i=0; i<9 ; i++) {
					for(int j=0; j<5 ; j++) {
						if(j==2)
							continue;
						hostFitnessTable[i][j] = 1-V;
						parasiteFitnessTable[i][j] = 1;
					}
				}
				hostFitnessTable[0][0]=hostFitnessTable[2][1]=hostFitnessTable[6][3]=hostFitnessTable[8][4] = 1;
				hostFitnessTable[1][0]=hostFitnessTable[3][0]=hostFitnessTable[1][1]=hostFitnessTable[5][1]=hostFitnessTable[3][3]=hostFitnessTable[7][3]=hostFitnessTable[5][4]=hostFitnessTable[7][4] = 1 - V - (D1 + (1-D1)*D2)*V;
				hostFitnessTable[4][0]=hostFitnessTable[4][1]=hostFitnessTable[4][3]=hostFitnessTable[4][4]= 1 - V - pow(D1 + (1-D1)*D2, 2.0)*V;
				parasiteFitnessTable[0][0]=parasiteFitnessTable[2][1]=parasiteFitnessTable[6][3]=parasiteFitnessTable[8][4] = 1 - T_COST;
				parasiteFitnessTable[1][0]=parasiteFitnessTable[3][0]=parasiteFitnessTable[1][1]=parasiteFitnessTable[5][1]=parasiteFitnessTable[3][3]=parasiteFitnessTable[7][3]=parasiteFitnessTable[5][4]=parasiteFitnessTable[7][4] = 1 - (D1+(1-D1)*D2)*T_COST;
				parasiteFitnessTable[4][0]=parasiteFitnessTable[4][1]=parasiteFitnessTable[4][3]=parasiteFitnessTable[4][4] = 1 - pow(D1 + (1-D1)*D2, 2.0)*T_COST;
			}
		}


		// In the following parts we initialize Parasite frequencies

		// use GSL library random num geberator to generate random frequencies
		const gsl_rng_type * T;
		gsl_rng * r;
		// Initialize the random number generator
		T = gsl_rng_taus;
		r = gsl_rng_alloc (T);
		gsl_rng_set (r, gsl_rng_default_seed);

		double freqA = gsl_ran_flat(r, 0.1, 0.9);
		double freq_a =  1-freqA;
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
		double freqB = gsl_ran_flat(r, 0.1, 0.9);
		double freq_b = 1-freqB;
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);

        // free up the memory taken by the random num generator
		gsl_rng_free(r);

		for(int i=0; i<9 ; i++) {
			parasites[i] = 0;
		}

		if(PARASITE_PLOIDY == 1) {
			parasites[0] = freqA * freqB;
			parasites[1] = freqA * freq_b;
			parasites[3] = freq_a * freqB;
			parasites[4] = freq_a * freq_b;
		}
		else if(PARASITE_PLOIDY == 2) {
			for(int i=0; i<9 ; i++) {
				parasites[i] = 1;
				for(int j=0; j<(int)(i/3) ; j++) {
					parasites[i] *= freq_a;
				}
				for(int j=0; j<2-(int)(i/3) ; j++) {
					parasites[i] *= freqA;
				}
				for(int j=0; j<i-((int)(i/3)*3) ; j++) {
					parasites[i] *= freq_b;
				}
				for(int j=0; j<2-i+((int)(i/3)*3) ; j++) {
					parasites[i] *= freqB;
				}
			}
		}
	}


	// initialize tempHostfreq. The values for immunityGenotypeFreq is calculated later in the code
	for(int i=0; i<9 ; i++) {
		immunityGenotypeFreq[i] = tempHostFreq[i] = 0;
	}


	for(int i=0; i<POPULATION; i++) {
		gsl_rng_default_seed += time(NULL)^i;
		indiv newIndiv(currentRound);

		immunityGenotypeFreq[(int)((newIndiv.immunityGeneA[0]+newIndiv.immunityGeneA[1])*3 + newIndiv.immunityGeneB[0]+newIndiv.immunityGeneB[1])] += 1;

		thePop.push_back(newIndiv);
	}


	for(int i=0; i<9 ; i++) {
		immunityGenotypeFreq[i] /= POPULATION;
	}

	// initialize other variables
	highestNewW=0;
	highestNewDadW = 0;
	meanRecombination = 1;
	newMeanRecombination = 0;
	selfedChildren[0] = 0.0;
	selfedChildren[1] = 0.0;
	outcrossedChildren[0] = 0.0;
	outcrossedChildren[1] = 0.0;
	population = POPULATION;
}

// Construct a population using a list of known individuals
studyPop::studyPop(vector<indiv> initPop) {

	thePop = initPop;

	if(initPop.size() == 0) {
		highestW = 0;
		highestDadW = 1 - K*INITIALSELFING;
		highestNewW=0;
		highestNewDadW = 0;
		population = 0;
		selfedChildren[0] = 0.0;
		selfedChildren[1] = 0.0;
		outcrossedChildren[0] = 0.0;
		outcrossedChildren[1] = 0.0;
	}

	else {

		highestW = thePop.at(0).get_fitness();
		for(int i=1; i<thePop.size() ; i++) {
			if(thePop.at(i).get_fitness() > highestW)
				highestW = thePop.at(i).get_fitness();
		}


		highestDadW = (1 - K*thePop.at(0).calcSelfingRate())*thePop.at(0).get_fitness();
		for(int i=1; i<thePop.size() ; i++) {
			if((1 - K*thePop.at(i).calcSelfingRate())*thePop.at(i).get_fitness() > highestDadW)
				highestDadW = (1 - K*thePop.at(i).calcSelfingRate())*thePop.at(i).get_fitness();
		}


		highestNewW = 0;
		highestNewDadW = 0;
		selfedChildren[0] = 0.0;
		selfedChildren[1] = 0.0;
		outcrossedChildren[0] = 0.0;
		outcrossedChildren[1] = 0.0;
		population = POPULATION;
		}
}




// This function is called to make next generation from existing indivs in a "studyPop"
void studyPop::makeNewGeneration(int currentGen, int currentRound) {

	if(PARASITES_EXIST) {
		makeNewParasiteGeneration();
		addMutationToImmunityGenes(1);
		addMutationToImmunityGenes(2);
	}

	for(int i=0; i<POPULATION; i++) {
		// The indiv made here is stored in "newGen" vector
		gsl_rng_default_seed += time(NULL);
		makeIndiv(currentGen, currentRound);
	}

	// Erase the last population from memory and copy new generation into "thePop" vector
	thePop = newGen;
	newGen.clear();
	population /= 2;

	// update values corresponding to the new generation & prepare the two temp varialbles for later use (make them zero)
	highestW = highestNewW;
	highestDadW = highestNewDadW;

	meanRecombination = newMeanRecombination/(2*POPULATION);
	for(int i=0; i<9 ; i++) {
		immunityGenotypeFreq[i] = tempHostFreq[i]/POPULATION;
		tempHostFreq[i] = 0;
	}

	highestNewW = 0;
	highestNewDadW = 0;
	newMeanRecombination = 0;

}





// pick an individual for reproduction
int studyPop::pickAMom() {

	// This value will later be returned to determine the index of mom in "thePop"
	int mom;

	// this is used for fitness check in this function
	double requiredFitness;

	// use GSL to generate random numbers
	const gsl_rng_type * T;
	gsl_rng * r;

	// Initialize the random number generator
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed);


	// pick different moms until one it passes fitness check
	do{
		mom = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
		// change seed for better randomization
		gsl_rng_default_seed += time(NULL)^mom;
		gsl_rng_set (r, gsl_rng_default_seed);
		requiredFitness = gsl_ran_flat(r, 0, 1);
	}while(requiredFitness > (thePop.at(mom).get_fitness()/highestW));

	gsl_rng_free(r);

	return mom;
}




// pick an individual for reproduction with an outcrossing indiv
int studyPop::pickADad(int chosenMom) {

	// This value will be returned to determine the index of dad in "thePop"
	int dad;

	// perform fitness check for outcrossing in dad
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_default_seed += time(NULL);
	// Initialize the random number generator
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed);

	double requiredFitnessToOutcross;

	// pick different dad until one it passes fitness check & is also not the same as the mom chosen already
	do{
		do{
			dad = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
		}while(dad == chosenMom);
		// change seed for more randomization
		gsl_rng_default_seed += time(NULL)^dad;
		gsl_rng_set (r, gsl_rng_default_seed);
		requiredFitnessToOutcross = gsl_ran_flat(r, 0, 1);
	}while(requiredFitnessToOutcross > (((1-K*thePop.at(dad).calcSelfingRate())*thePop.at(dad).get_fitness())/highestDadW));

    // free up the memory taken by the random number generator
	gsl_rng_free(r);

	return dad;
}





void studyPop::makeIndiv(int currentGen, int currentRound) {


	int mom, dad;
	vector<int> *momChromosome, *dadChromosome;
	// variables to save R,S,N, immuneA & immuneB alleles inherited from parents
	double momS, momR, dadS, dadR, momN, dadN;
	int momA, momB, dadA, dadB;
	// MH: Variable to determine if gamete is asexually-produced
	int isasex = 0;

	// pick the mom for the new indiv
	mom = pickAMom();

	// check if the mom selfs or outcrosses MH OR ASEXUAL
	const gsl_rng_type * T;
	gsl_rng * r;
	// Initialize the random number generator
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL)^mom;
	gsl_rng_set (r, gsl_rng_default_seed);

	double outCrossCheck = gsl_ran_flat(r, 0, 1);
	double AsexCheck = gsl_ran_flat(r, 0, 1);
	/*
	printf("\n");
	printf("Asex Check is %f\n",AsexCheck);
	*/

	if(AsexCheck < ASEXP){
		// In this case, offspring is clonal
		dad = mom;
		isasex = 1;
	} else {
		if(outCrossCheck > thePop.at(mom).calcSelfingRate()) {
			// in this case mom would outcross
			dad = pickADad(mom);
		}
		else {
			// in this case mom will self
			dad = mom;
		}
	}
	/*
	printf("IsAsex is %d\n",isasex);
	printf("Dad is %d, Mom is %d\n",dad,mom);
	*/

	// Get the chromosomes from mom and dad
	// The last five elements are not mutations (check makeHamete function in individual.cpp to find out)
	/* printf("Creating Mom chromosome\n"); */
	momChromosome = thePop.at(mom).makeGamete(this->meanRecombination,isasex,0);
	/* printf("Creating Dad chromosome\n"); */
	dadChromosome = thePop.at(dad).makeGamete(this->meanRecombination,isasex,1);

	// renew the seed
	gsl_rng_default_seed += time(NULL)*(momChromosome->size() + dadChromosome->size());
	gsl_rng_set (r, gsl_rng_default_seed);


	// A variable to store the mutation amount of Selfing, recomb and neutral genes
	double geneMutation;

	momN = *(thePop.at(mom).neutralGene + momChromosome->at(momChromosome->size()-5));
	// After mutation check, add mutation to inherited neutral gene
	double randomMutationChecker;
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<NEUTRALMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAN);
		momN = momN + geneMutation;
		gsl_rng_default_seed += time(NULL)^(mom);
		gsl_rng_set (r, gsl_rng_default_seed);
	}


	momS = *(thePop.at(mom).selfingGene + momChromosome->at(momChromosome->size()-4));
	// After mutation check, add mutation to inherited selfing gene and keep it in [0,1]
	// MH: I've altered this to remove the dependence on time, so no selfing arises if not defined
	randomMutationChecker = gsl_ran_flat(r, 0, 1);

	/*if(currentGen <= 2000)
		randomMutationChecker = 1;*/
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	
	if(randomMutationChecker<SELFINGMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAS);
		momS = ((momS + geneMutation) > 0 ? ((momS + geneMutation) < 1 ? (momS + geneMutation) : 1) : 0);
		gsl_rng_default_seed += time(NULL)^((int)floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5)+1);
		gsl_rng_set (r, gsl_rng_default_seed);
	}


	momR = *(thePop.at(mom).recombinationGene + momChromosome->at(momChromosome->size()-3));
	// After mutation check, add mutations to inherited recombination gene and keep it >0
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<RECOMBMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAR);
		momR = (momR + geneMutation)>0 ? (momR + geneMutation) : 0;
		gsl_rng_default_seed += time(NULL)^((int)floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5)+1);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	// mutation is not added here, it is added after a generation is made
	momA = *(thePop.at(mom).immunityGeneA + momChromosome->at(momChromosome->size()-2));
	momB = *(thePop.at(mom).immunityGeneB + momChromosome->at(momChromosome->size()-1));


	dadN = *(thePop.at(dad).neutralGene + dadChromosome->at(dadChromosome->size()-5));
	// After mutation check, add mutation to inherited neutral gene
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<NEUTRALMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAN);
		dadN = dadN + geneMutation;
		gsl_rng_default_seed += time(NULL)^(dad);
		gsl_rng_set (r, gsl_rng_default_seed);
	}


	dadS = *(thePop.at(dad).selfingGene + dadChromosome->at(dadChromosome->size()-4));
	// add mutation to inherited selfing gene and keep it in [0,1]
	// MH: I've altered this to remove the dependence on time, so no selfing arises if not defined
	randomMutationChecker = gsl_ran_flat(r, 0, 1);

	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

	if(randomMutationChecker<SELFINGMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r,SIGMAS);
		dadS = ((dadS + geneMutation) > 0 ? ((dadS + geneMutation) < 1 ? (dadS + geneMutation) : 1) : 0);
		gsl_rng_default_seed += time(NULL)^((int)floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5)+1);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	dadR = *(thePop.at(dad).recombinationGene + dadChromosome->at(dadChromosome->size()-3));
	// add mutations to inherited recombination gene and keep it >0
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<RECOMBMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAR);
		dadR = (dadR + geneMutation)>0 ? (dadR + geneMutation) : 0;
		gsl_rng_default_seed += time(NULL)^((int)floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5)+1);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	// mutation for these alleles are added after a generation is made
	dadA = *(thePop.at(dad).immunityGeneA + dadChromosome->at(dadChromosome->size()-2));
	dadB = *(thePop.at(dad).immunityGeneB + dadChromosome->at(dadChromosome->size()-1));

	// removing the last five elements which were just indices of inherited R,S,N,immuneA,immuneB
	momChromosome->pop_back();
	momChromosome->pop_back();
	momChromosome->pop_back();
	momChromosome->pop_back();
	momChromosome->pop_back();
	dadChromosome->pop_back();
	dadChromosome->pop_back();
	dadChromosome->pop_back();
	dadChromosome->pop_back();
	dadChromosome->pop_back();

	// Add deleterious and beneficial mutations to the chromosomes inherited from parents
	addMutation(*momChromosome, *dadChromosome, currentRound);

    // create the new indiv with the given info
	indiv theChild(*momChromosome, *dadChromosome, momS, dadS, momR, dadR, momN, dadN, momA, dadA, momB, dadB);

	tempHostFreq[(int)((theChild.immunityGeneA[0]+theChild.immunityGeneA[1])*3 + theChild.immunityGeneB[0]+theChild.immunityGeneB[1])] += 1;

	// save whether this child is a result of selfing or outcrossing
	// these values will later be used to calculate inbreeding dep if "EXPERIMENTAL_INBREEDING_DEP" is zero
	// NOTE: in our simulations "EXPERIMENTAL_INBREEDING_DEP" was 1
	if(mom == dad && isasex == 0) {
		selfedChildren[0] += 1.0;
		selfedChildren[1] += theChild.get_fitness();
	}
	else {
		outcrossedChildren[0] += 1.0;
		outcrossedChildren[1] += theChild.get_fitness();
	}


	newMeanRecombination += momR + dadR;

	// use delete because we used new in "makeGamete" method
	delete momChromosome;
	delete dadChromosome;

	newGen.push_back(theChild);

    // free up the space taken by the random number generator
	gsl_rng_free(r);

    // one indiv was added to the new population (new generation)
	population++;
}


// make a child with given parents (in case of selfing, both parents are the same indiv)
// MH NOTE: I've yet to update this section to account for asex: but if it is only used 
// when testing ID, then I suspect accounting for asex is not needed...?
indiv* studyPop::makeForcedIndiv(int mom, int dad, int currentRound) {

	vector<int> *momChromosome, *dadChromosome;
	// variables to save R,S,N,immuneA and immuneB inherited from parents
	double momS, momR, dadS, dadR, momN, dadN;
	int momA, momB, dadA, dadB;



	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);



	// Get the chromosomes from mom and dad
	momChromosome = thePop.at(mom).makeGamete(this->meanRecombination,0,0);
	dadChromosome = thePop.at(dad).makeGamete(this->meanRecombination,0,1);

	// renew the seed for more randomization
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

	// A variable to store the mutation amount of Selfing, recomb and neutral genes (R,S,N)
	double geneMutation;

	momN = *(thePop.at(mom).neutralGene + momChromosome->at(momChromosome->size()-5));
	// After mutation check, add mutation to inherited neutral gene
	double randomMutationChecker;
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<NEUTRALMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAN);
		momN = momN + geneMutation;
		gsl_rng_default_seed += time(NULL)^(mom);
		gsl_rng_set (r, gsl_rng_default_seed);
	}




	momS = *(thePop.at(mom).selfingGene + momChromosome->at(momChromosome->size()-4));
	// After mutation check, add mutation to inherited selfing gene and keep it in [0,1]
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<SELFINGMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAS);
		momS = ((momS + geneMutation) > 0 ? ((momS + geneMutation) < 1 ? (momS + geneMutation) : 1) : 0);
		gsl_rng_default_seed += time(NULL)^((int)floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5)+1);
		gsl_rng_set (r, gsl_rng_default_seed);
	}


	momR = *(thePop.at(mom).recombinationGene + momChromosome->at(momChromosome->size()-3));
	// After mutation check, add mutations to inherited recombination gene and keep it >0
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<RECOMBMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAR);
		momR = (momR + geneMutation)>0 ? (momR + geneMutation) : 0;
		gsl_rng_default_seed += time(NULL)^((int)floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5)+1);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	// mutation is not added here, it is added later, after a generation is made
	momA = *(thePop.at(mom).immunityGeneA + momChromosome->at(momChromosome->size()-2));
	momB = *(thePop.at(mom).immunityGeneB + momChromosome->at(momChromosome->size()-1));


	dadN = *(thePop.at(dad).neutralGene + dadChromosome->at(dadChromosome->size()-5));
	// After mutation check, add mutation to inherited neutral gene
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<NEUTRALMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAN);
		dadN = dadN + geneMutation;
		gsl_rng_default_seed += time(NULL)^(dad);
		gsl_rng_set (r, gsl_rng_default_seed);
	}


	dadS = *(thePop.at(dad).selfingGene + dadChromosome->at(dadChromosome->size()-4));
	// add mutation to inherited selfing gene and keep it in [0,1]
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<SELFINGMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAS);
		dadS = ((dadS + geneMutation) > 0 ? ((dadS + geneMutation) < 1 ? (dadS + geneMutation) : 1) : 0);
		gsl_rng_default_seed += time(NULL)^((int)floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5)+1);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	dadR = *(thePop.at(dad).recombinationGene + dadChromosome->at(dadChromosome->size()-3));
	// add mutations to inherited recombination gene and keep it >0
	randomMutationChecker = gsl_ran_flat(r, 0, 1);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);
	if(randomMutationChecker<RECOMBMUTATIONRATE) {
		geneMutation = gsl_ran_gaussian(r, SIGMAR);
		dadR = (dadR + geneMutation)>0 ? (dadR + geneMutation) : 0;
		gsl_rng_default_seed += time(NULL)^((int)floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5)+1);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	// mutation for the immunity alleles are added after a generation is made
	dadA = *(thePop.at(dad).immunityGeneA + dadChromosome->at(dadChromosome->size()-2));
	dadB = *(thePop.at(dad).immunityGeneB + dadChromosome->at(dadChromosome->size()-1));

	// removing the last three elements which were just indices of inherited recomb & selfing & neutral genes
	momChromosome->pop_back();
	momChromosome->pop_back();
	momChromosome->pop_back();
	momChromosome->pop_back();
	momChromosome->pop_back();
	dadChromosome->pop_back();
	dadChromosome->pop_back();
	dadChromosome->pop_back();
	dadChromosome->pop_back();
	dadChromosome->pop_back();

	// Add mutations to the chromosomes inherited from parents
	addMutation(*momChromosome, *dadChromosome, currentRound);

    // create the child with the given info
	indiv* theChild = new indiv(*momChromosome, *dadChromosome, momS, dadS, momR, dadR, momN, dadN, momA, dadA, momB, dadB);

	// use delete because we used new in "makeGamete" method
	delete momChromosome;
	delete dadChromosome;

    // free up the memory taken by the random number generator
	gsl_rng_free(r);

	return theChild;
}




void studyPop::addMutation(vector<int>& parent1, vector<int>& parent2, int currentRound) {

	// These variables save the average number of mutations on each chromosome
	double avgDelMutations;
	double avgBenMutations;
	double P_ben;

    P_ben = PROB_BEN;

	avgDelMutations = (0.5)*MUTATIONRATE*(1-P_ben); // multiplied by half to get # for each chromosome
	avgBenMutations = (0.5)*MUTATIONRATE*P_ben;   // multiplied by half to get # for each chromosome

	// Pick number of mutation loci on each gamete in the next several lines
	int parent1DelMutationsNum, parent1BenMutationsNum, parent2DelMutationsNum, parent2BenMutationsNum;

	// GSL objects to generate random numbers
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_default_seed += time(NULL);
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r, gsl_rng_default_seed);

	parent1DelMutationsNum = gsl_ran_poisson(r,avgDelMutations);
	// Get a new seed for more randomization
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set(r, gsl_rng_default_seed);
	parent1BenMutationsNum = gsl_ran_poisson(r,avgBenMutations);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set(r, gsl_rng_default_seed);
	parent2DelMutationsNum = gsl_ran_poisson(r,avgDelMutations);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set(r, gsl_rng_default_seed);
	parent2BenMutationsNum = gsl_ran_poisson(r,avgBenMutations);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set(r, gsl_rng_default_seed);

	// Pick exact loci of mutation on each chromosome & then put them on the chromosome
	vector<int> allMutationsOnGamete1, benMutationsOnGamete1, delMutationsOnGamete1;
	vector<int> allMutationsOnGamete2, benMutationsOnGamete2, delMutationsOnGamete2;
	// Gamete from parent1
	for(int i=0; i<parent1DelMutationsNum ; i++) {
		gsl_rng_default_seed += time(NULL)^i;
		gsl_rng_set(r, gsl_rng_default_seed);
		delMutationsOnGamete1.push_back(floor(gsl_ran_flat(r, 0.5, L+0.4999) + 0.5));
	}

	// Find beneficial mutations loci
	for(int i=0; i<parent1BenMutationsNum ;) {
		gsl_rng_default_seed += time(NULL)^i;
		gsl_rng_set(r, gsl_rng_default_seed);
		int tempMutation = floor(gsl_ran_flat(r, 0.5, L+0.4999) + 0.5);
		int q=0;
		for(; q<delMutationsOnGamete1.size() ; q++) {
			if(delMutationsOnGamete1.at(q) == tempMutation)
				break;
		}
		if(q==delMutationsOnGamete1.size()) {
			benMutationsOnGamete1.push_back(tempMutation);
			i++;
		}
	}

    // put the beneficial and del mutations in a same vector & sort
	allMutationsOnGamete1.insert(allMutationsOnGamete1.end(), delMutationsOnGamete1.begin(), delMutationsOnGamete1.end());
	allMutationsOnGamete1.insert(allMutationsOnGamete1.end(), benMutationsOnGamete1.begin(), benMutationsOnGamete1.end());

	sort(allMutationsOnGamete1.begin(), allMutationsOnGamete1.end());

    // make the beneficial mutations negative (and leave the the del mutations positive integers)
	for(int p=0; p<allMutationsOnGamete1.size() ; p++) {
		for(int q=0; q<benMutationsOnGamete1.size() ; q++) {
			if(allMutationsOnGamete1.at(p) == benMutationsOnGamete1.at(q)) {
				allMutationsOnGamete1.at(p) *= -1;
			}
		}
	}

	//Merge the sorted mutations with the unmutated chromosome
	int placeOnChromosome = 0;
	if(parent1.size() > 0) {
		for(int i=0; i<allMutationsOnGamete1.size();) {
			if(abs(allMutationsOnGamete1.at(i)) < abs(parent1.at(placeOnChromosome))) {
				parent1.insert(parent1.begin() + placeOnChromosome, allMutationsOnGamete1.at(i));
				i++;
				placeOnChromosome++;
			}
			else if(abs(allMutationsOnGamete1.at(i)) > abs(parent1.at(placeOnChromosome))) {
				if(placeOnChromosome < parent1.size()-1) {
					placeOnChromosome++;
				}
				else {
					parent1.insert(parent1.begin() + placeOnChromosome + 1, allMutationsOnGamete1.begin() + i, allMutationsOnGamete1.end());
					break;
				}
			}
			// This else is for when the mutation already exists on the chromosome.
			else {
				if(allMutationsOnGamete1.at(i) == parent1.at(placeOnChromosome))
					i++;
				else{
					// replace beneficial/deleterious with each other
					allMutationsOnGamete1.at(i) *= -1;
					i++;
				}
			}

		}
	}

	// this else is executed when original parent chromosome has no mutations at all
	else
		parent1 = allMutationsOnGamete1;





	// Gamete from parent2
	for(int i=0; i<parent2DelMutationsNum ; i++) {
		gsl_rng_default_seed += time(NULL)^i;
		gsl_rng_set(r, gsl_rng_default_seed);
		delMutationsOnGamete2.push_back(floor(gsl_ran_flat(r, 0.5, L+0.4999) + 0.5));
	}


	// Find beneficial mutations loci
	for(int i=0; i<parent2BenMutationsNum ;) {
		gsl_rng_default_seed += time(NULL)^i;
		gsl_rng_set(r, gsl_rng_default_seed);
		int tempMutation = floor(gsl_ran_flat(r, 0.5, L+0.4999) + 0.5);
		int q=0;
		for(; q<delMutationsOnGamete2.size() ; q++) {
			if(delMutationsOnGamete2.at(q) == tempMutation)
				break;
		}
		if(q==delMutationsOnGamete2.size()) {
			benMutationsOnGamete2.push_back(tempMutation);
			i++;
		}
	}

    // put all the mutations in a same vector and sort
	allMutationsOnGamete2.insert(allMutationsOnGamete2.end(), delMutationsOnGamete2.begin(), delMutationsOnGamete2.end());
	allMutationsOnGamete2.insert(allMutationsOnGamete2.end(), benMutationsOnGamete2.begin(), benMutationsOnGamete2.end());

	sort(allMutationsOnGamete2.begin(), allMutationsOnGamete2.end());

    // make the ben mutations negative
	for(int p=0; p<allMutationsOnGamete2.size() ; p++) {
		for(int q=0; q<benMutationsOnGamete2.size() ; q++) {
			if(allMutationsOnGamete2.at(p) == benMutationsOnGamete2.at(q)) {
				allMutationsOnGamete2.at(p) *= -1;
			}
		}
	}


	//Merge the sorted mutations with the unmutated chromosome
	placeOnChromosome = 0;
	if(parent2.size() > 0) {
		for(int i=0; i<allMutationsOnGamete2.size();) {
			if(abs(allMutationsOnGamete2.at(i)) < abs(parent2.at(placeOnChromosome))) {
				parent2.insert(parent2.begin() + placeOnChromosome, allMutationsOnGamete2.at(i));
				i++;
				placeOnChromosome++;
			}
			else if(abs(allMutationsOnGamete2.at(i)) > abs(parent2.at(placeOnChromosome))) {
				if(placeOnChromosome < parent2.size()-1) {
					placeOnChromosome++;
				}
				else {
					parent2.insert(parent2.begin() + placeOnChromosome + 1, allMutationsOnGamete2.begin() + i, allMutationsOnGamete2.end());
					break;
				}
			}
			// This else is for when the mutation already exists on the chromosome.
			else {
				if(allMutationsOnGamete2.at(i) == parent2.at(placeOnChromosome))
					i++;
				else{
					// replace beneficial/deleterious with each other
					allMutationsOnGamete2.at(i) *= -1;
					i++;
				}

			}

		}

	}

	// this else is executed when original parent chromosome has no mutations at all
	else {
		parent2 = allMutationsOnGamete2;
	}

    // free up the memory taken by the random number generator
	gsl_rng_free(r);
}

void studyPop::addParasitesEffect(indiv& newIndiv) {
    // get the immunity gene info from the indiv
	int hostGenotypeCode = (newIndiv.immunityGeneA[0] + newIndiv.immunityGeneA[1]) * 3;
	hostGenotypeCode += newIndiv.immunityGeneB[0] + newIndiv.immunityGeneB[1];

	long double fitnessReductionByParasites=0;

    // use the host fitness table (initialized in the constructor function) to find reduction in fitness
	if(PARASITE_PLOIDY == 1) {
		for(int i=0; i<5 ; i++) {
			if(i != 2) {
				fitnessReductionByParasites += hostFitnessTable[hostGenotypeCode][i] * parasites[i];
			}
		}
	}

	newIndiv.fitness *= fitnessReductionByParasites;
}




void studyPop::makeNewParasiteGeneration() {

	// calculate the fitness of each genotype
	double parasiteGenotypeFitness[9];
	for(int i=0; i<9 ; i++) {
		parasiteGenotypeFitness[i] = 0;
	}
	if(PARASITE_PLOIDY == 1) {
		for(int i=0; i<5 ; i++) {
			if(i != 2) {
				for(int j=0; j<9 ; j++) {
					parasiteGenotypeFitness[i] += immunityGenotypeFreq[j]*parasiteFitnessTable[j][i];
				}
			}
		}
	}
	else if(PARASITE_PLOIDY == 2) {
		for(int i=0; i<9 ; i++) {
			for(int j=0; j<9 ; j++) {
				parasiteGenotypeFitness[i] += immunityGenotypeFreq[j]*parasiteFitnessTable[j][i];
			}
		}
	}


	// Calculate the average fitness in the parasite population
	double meanFitness=0;
	for(int i=0 ; i<9 ; i++) {
		meanFitness += parasites[i]*parasiteGenotypeFitness[i];
	}

	// Calculate parasite genotype freq after selection
	double parasitesAfterSelection[9];
	if(PARASITE_PLOIDY == 1) {
		for(int i=0; i<9 ; i++) {
			if(i<5 && i!=2) {
				parasitesAfterSelection[i] = parasites[i]*(parasiteGenotypeFitness[i]/meanFitness);
			}
			else {
				parasitesAfterSelection[i] = 0;
			}
		}
	}
	else if(PARASITE_PLOIDY == 2) {
		for(int i=0; i<9 ; i++) {
			parasitesAfterSelection[i] = parasites[i]*(parasiteGenotypeFitness[i]/meanFitness);
		}
	}

	// Calculate allele freq after selection
	double freqOfA1=0, freqOfB1=0;
	double freqOfA2, freqOfB2;
	if(PARASITE_PLOIDY == 1) {
		freqOfA1 = parasitesAfterSelection[0] + parasitesAfterSelection[1];
		freqOfB1 = parasitesAfterSelection[0] + parasitesAfterSelection[3];
	}
	else if (PARASITE_PLOIDY == 2) {
		freqOfA1 = parasitesAfterSelection[0]+parasitesAfterSelection[1]+parasitesAfterSelection[2] + 0.5*(parasitesAfterSelection[3]+parasitesAfterSelection[4]+parasitesAfterSelection[5]);
		freqOfB1 = parasitesAfterSelection[0]+parasitesAfterSelection[3]+parasitesAfterSelection[6] + 0.5*(parasitesAfterSelection[1]+parasitesAfterSelection[4]+parasitesAfterSelection[7]);
	}
	freqOfA2 = 1 - freqOfA1;
	freqOfB2 = 1 - freqOfB1;

	// Calculate the linkage disequilibrium after selection ASSUMING parasites are haploid
	double LD = parasitesAfterSelection[0]*parasitesAfterSelection[4] - parasitesAfterSelection[1]*parasitesAfterSelection[3];

	// Calculate genotype frequencies for parasites after recombination
	double freq_AB_afterR = freqOfA1*freqOfB1 - EFFECTIVE_PARASITE_R*LD;
	double freq_Ab_afterR = freqOfA1*freqOfB2 + EFFECTIVE_PARASITE_R*LD;
	double freq_aB_afterR = freqOfA2*freqOfB1 + EFFECTIVE_PARASITE_R*LD;
	double freq_ab_afterR = freqOfA2*freqOfB2 - EFFECTIVE_PARASITE_R*LD;

	// Calculate genotype frequencies for parasites after mutation (final genotype frequencies)
	parasites[0] = pow((1-U_PARASITE),2.0)*freq_AB_afterR + (1-U_PARASITE)*U_PARASITE*(freq_Ab_afterR + freq_aB_afterR) + pow(U_PARASITE, 2.0)*freq_ab_afterR;
	parasites[1] = pow((1-U_PARASITE),2.0)*freq_Ab_afterR + (1-U_PARASITE)*U_PARASITE*(freq_AB_afterR + freq_ab_afterR) + pow(U_PARASITE, 2.0)*freq_aB_afterR;
	parasites[3] = pow((1-U_PARASITE),2.0)*freq_aB_afterR + (1-U_PARASITE)*U_PARASITE*(freq_AB_afterR + freq_ab_afterR) + pow(U_PARASITE, 2.0)*freq_Ab_afterR;
	parasites[4] = pow((1-U_PARASITE),2.0)*freq_ab_afterR + (1-U_PARASITE)*U_PARASITE*(freq_Ab_afterR + freq_aB_afterR) + pow(U_PARASITE, 2.0)*freq_AB_afterR;
}

void studyPop::addMutationToImmunityGenes(int whichLocus) {
	int avgNumberOfHostMutations = 2*POPULATION*U_HOST_IMMUNITY;
	int numOfHostMutations;

	// GSL objects to generate random numbers
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_default_seed += time(NULL);
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r, gsl_rng_default_seed);

	numOfHostMutations = gsl_ran_poisson(r,avgNumberOfHostMutations);

	// Find the alleles in the population to be mutated
	vector<int> indexOfMutations;
	for(int i=0; i<numOfHostMutations ;) {
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set(r, gsl_rng_default_seed);
		int tempIndex = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
		int q=0;
		for(; q<indexOfMutations.size() ; q++) {
			if(indexOfMutations.at(q) == tempIndex)
				break;
		}
		if(q==indexOfMutations.size()) {
			indexOfMutations.push_back(tempIndex);
			i++;
		}
	}

	// Now we apply mutations to indiv with indices found above
	for(int i=0; i<indexOfMutations.size() ; i++) {
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set(r, gsl_rng_default_seed);
		int whichAllele = floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5);
		if(whichLocus == 1) {
			thePop.at(indexOfMutations.at(i)).immunityGeneA[(int) whichAllele] += 1;
			thePop.at(indexOfMutations.at(i)).immunityGeneA[(int) whichAllele] %= 2;
		}
		else if(whichLocus == 2) {
			thePop.at(indexOfMutations.at(i)).immunityGeneB[(int) whichAllele] += 1;
			thePop.at(indexOfMutations.at(i)).immunityGeneB[(int) whichAllele] %= 2;
		}
	}
	gsl_rng_default_seed += time(NULL);

    // free up the space taken by the random number generator
	gsl_rng_free(r);
}





// few functions used to get data from studyPop class
int studyPop::get_population() {
	return population;
}



vector<indiv>* studyPop::get_thePop() {
	return &thePop;
}

vector<indiv>* studyPop::get_newGen() {
	return &newGen;
}
