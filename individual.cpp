#include "individual.h"

// a helper function to compare two integers by their absolute values
bool isSmaller(int x, int y) {
	return abs(x)<abs(y);
}


// This constructor is called only for the individuals in the first generation
indiv::indiv(int currentRound) {

	// Initialize selfing rate values
	selfingGene[0] = INITIALSELFING;
	selfingGene[1] = INITIALSELFING;

	// Initialize the random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_taus;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed);

    // if the constant "INITIAL_SELFING_VARIATION" is set to 1, apply one mutation to each copy of the selfing allele in the population
	if(INITIAL_SELFING_VARIATION == 1) {
		double selfingGeneMutation;

        // find a mutation (a real value). If the mutation makes the selfing <0 then put 0. If it makes a value >1 then put 1
		selfingGeneMutation = gsl_ran_gaussian(r, SIGMAS);
		selfingGene[0] = ((selfingGene[0] + selfingGeneMutation) > 0 ? ((selfingGene[0] + selfingGeneMutation) < 1 ? (selfingGene[0] + selfingGeneMutation) : 1) : 0);
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
		selfingGeneMutation = gsl_ran_gaussian(r, SIGMAS);
		selfingGene[1] = ((selfingGene[1] + selfingGeneMutation) > 0 ? ((selfingGene[1] + selfingGeneMutation) < 1 ? (selfingGene[1] + selfingGeneMutation) : 1) : 0);
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	// Initialize recombination rate values
	recombinationGene[0] = INITIALRECOMBINATION;
	recombinationGene[1] = INITIALRECOMBINATION;

	// Initialize neutral gene values
	neutralGene[0] = INITIALNEUTRALGENE;
	neutralGene[1] = INITIALNEUTRALGENE;

	// Initialize immunity gene values
	for(int i=0; i<2 ; i++) {
		double randomMutationChecker = gsl_ran_flat(r, 0, 1);
		// randomely make 50% of the alleles 0 and 50% of them 1
		if(randomMutationChecker<0.5) {
			immunityGeneA[i] = 1;
		}
		else {
			immunityGeneA[i] = 0;
		}
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	for(int i=0; i<2 ; i++) {
		double randomMutationChecker = gsl_ran_flat(r, 0, 1);
		// randomely make 50% of the alleles 0 and 50% of them 1
		if(randomMutationChecker<0.5) {
			immunityGeneB[i] = 1;
		}
		else {
			immunityGeneB[i] = 0;
		}
		gsl_rng_default_seed += time(NULL);
		gsl_rng_set (r, gsl_rng_default_seed);
	}

	// initialize chromosomes
	// we get initial mutation number from a poisson distribution
	double meanInitialDelMutation;
	double meanInitialBenMutation;
	double P_ben;
    P_ben = PROB_BEN;

	meanInitialDelMutation = (MUTATIONRATE*(1-P_ben))/(H*S);
	meanInitialBenMutation = (MUTATIONRATE*P_ben)/(H*S);


	int numOfDelMutations = gsl_ran_poisson(r, meanInitialDelMutation);
	// change seed for more randomization
	gsl_rng_default_seed += time(NULL)*numOfDelMutations;
	gsl_rng_set (r, gsl_rng_default_seed);
	int numOfBenMutations = gsl_ran_poisson(r, meanInitialBenMutation);
	// change seed for more randomization
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);


	for(int i=0 ; i<numOfDelMutations ; i++) {
	    // decide on which locus the mutation will be
		int newMutation = floor(gsl_ran_flat(r, 0.5, L+0.4999) + 0.5);
		gsl_rng_default_seed += time(NULL)*newMutation;
		gsl_rng_set (r, gsl_rng_default_seed);
		// Decide on which chromosome the new mutation will be
		int whichChromosome = floor(gsl_ran_flat(r, 0.5, 2.4999) + 0.5);

		gsl_rng_default_seed += time(NULL)*whichChromosome;
		gsl_rng_set (r, gsl_rng_default_seed);

		// check that this mutation does not already exist on this chromosome
		// if not add it, if it does, get a new locus for the mutations
		if(whichChromosome == 1) {

			int pushMutation=1;

			//check if the mutation exists
			for(int k=0; k<theChromosome1.size() ; k++) {
				if(theChromosome1.at(k) == newMutation) {
					i--;
					pushMutation = 0;
					break;
				}
			}

			// if it does not exist already, add it
			if(pushMutation == 1)
				theChromosome1.push_back(newMutation);
		}


		// Check same things for chromosome2
		if(whichChromosome == 2) {
			int pushMutation=1;

			//check if the mutation exists
			for(int k=0; k<theChromosome2.size() ; k++) {
				if(theChromosome2.at(k) == newMutation) {
					i--;
					pushMutation = 0;
					break;
				}
			}

			// if it does not exist already, add it
			if(pushMutation == 1)
				theChromosome2.push_back(newMutation);
		}
	}


	// Add beneficial mutations
	for(int i=0 ; i<numOfBenMutations ; i++) {
	    // decide on which locus the mutation will be
		int newMutation = (-1)*floor(gsl_ran_flat(r, 0.5, L+0.4999) + 0.5);
		gsl_rng_default_seed += time(NULL)*(-1)*newMutation;
		gsl_rng_set (r, gsl_rng_default_seed);
		// Decide on which chromosome the new mutation will be
		int whichChromosome = floor(gsl_ran_flat(r, 0.5, 2.4999) + 0.5);

		gsl_rng_default_seed += time(NULL)*whichChromosome;
		gsl_rng_set (r, gsl_rng_default_seed);

		// check that this mutation does not already exist on this chromosome
		// if it does not, then add it
		// if the deleterious allele exists, then replace it
		if(whichChromosome == 1) {

			int pushMutation=1;

			//check if the mutation exists
			for(int k=0; k<theChromosome1.size() ; k++) {
				if(theChromosome1.at(k) == newMutation) {
					i--;
					pushMutation = 0;
					break;
				}
			}

			// if it does not exist already, add it
			// if the deleterious allele exists, replace it
			if(pushMutation == 1) {
				int replaced=0;
				for(int k=0;k<theChromosome1.size();k++) {
					if(theChromosome1.at(k) == abs(newMutation)) {
						theChromosome1.at(k) *= (-1);
						replaced = 1;
						break;
					}
				}
				if(replaced == 0)
					theChromosome1.push_back(newMutation);
			}
		}

		// Check same things for chromosome2
		if(whichChromosome == 2) {
			int pushMutation=1;

			//check if the mutation exists
			for(int k=0; k<theChromosome2.size() ; k++) {
				if(theChromosome2.at(k) == newMutation) {
					i--;
					pushMutation = 0;
					break;
				}
			}

			// if it does not exist already, add it
			// if the deleterious allele exists, replace it
			if(pushMutation == 1) {
				int replaced=0;
				for(int k=0;k<theChromosome2.size();k++) {
					if(theChromosome2.at(k) == abs(newMutation)) {
						theChromosome2.at(k) *= (-1);
						replaced = 1;
						break;
					}
				}
				if(replaced == 0)
					theChromosome2.push_back(newMutation);
			}
		}
	}

    // sort the mutations on the chromosomes by their absolute values
	sort(theChromosome1.begin(), theChromosome1.end(), isSmaller);
	sort(theChromosome2.begin(), theChromosome2.end(), isSmaller);

    // free up the space space taken by the random number generator
	gsl_rng_free(r);

    // use the "calcFitness" function to find the number of homozygous and heterozygous ben and del alleles (used to calc fitness)
	this->calcFitness();

}

// This constructor is called for indivs in second generation and so on. It is used to make new indiv from parent(s)
indiv::indiv(vector<int> chromosome1, vector<int> chromosome2, double selfing1, double selfing2, double recombination1, double recombination2, double neutral1, double neutral2, int A1, int A2, int B1, int B2) {
	// save the chromosomes that are inherited from parent(s)
	theChromosome1 = chromosome1;
	theChromosome2 = chromosome2;
	// save selfing rates inherited from parent(s)
	selfingGene[0] = selfing1;
	selfingGene[1] = selfing2;
	// save recombination rates inherited from parent(s)
	recombinationGene[0] = recombination1;
	recombinationGene[1] = recombination2;
	// save neutral gene values inherited from parents
	neutralGene[0] = neutral1;
	neutralGene[1] = neutral2;
	// save immunity genes values inherited from parents
	immunityGeneA[0] = A1;
	immunityGeneA[1] = A2;
	immunityGeneB[0] = B1;
	immunityGeneB[1] = B2;
	// use the "calcFitness" function to find the number of homozygous and heterozygous ben and del alleles (used to calc fitness)
	calcFitness();
}



// This function returns the selfing rate of this individual
float indiv::calcSelfingRate() {
	return ((selfingGene[0] +selfingGene[1])/2);
}



// This function returns the recombination rate of this individual
float indiv::calcRecombinationRate() {
	return ((recombinationGene[0] + recombinationGene[1])/2);
}




void indiv::calcFitness() {

	// number of deleterious heterozygous mutations
	int numOfDelHetero = 0;
	// number of deleterious homozygous mutations
	int numOfDelHomo = 0;
	// number of beneficial heterozygous mutations
	int numOfBenHetero = 0;
	// number of beneficial homozygous mutations
	int numOfBenHomo = 0;

	// this is used in the following loop
	int k=0;

	// in this part we only find the number of Homozygous mutations
	for(int i=0; i < theChromosome1.size(); i++) {
		for(int j=k; j < theChromosome2.size(); j=k) {
			if (abs(theChromosome1.at(i)) < abs(theChromosome2.at(j))) {
				// jump to check the next mutation on theChromosome1
				break;
			}
			else if (abs(theChromosome1.at(i)) > abs(theChromosome2.at(j))) {
				// continue moving on theChromosome2 and search if there is a homozygous mutation
				k++;
				continue;
			}
			else if (abs(theChromosome1.at(i)) == abs(theChromosome2.at(j))) {
				// found a new homozygote mutation
				if(theChromosome1.at(i)>0 && theChromosome2.at(j)>0)
					numOfDelHomo++;
				else if(theChromosome1.at(i)<0 && theChromosome2.at(j)<0)
					numOfBenHomo++;

				// jump to the next mutation loci on both chromosomes
				k++;
				break;
			}
		}

		if (k >= theChromosome2.size())
			break;
	}

	// we know the number of homozygous mutations
	// so we can calculate the number of heterozygous mutations
	int numOfHetero = (theChromosome1.size() + theChromosome2.size()) - (2*numOfDelHomo) - (2*numOfBenHomo);

	// Find the number of Del heterozygotes
	int numOfDel=0;
	for(int p=0; p<theChromosome1.size() ; p++) {
		if(theChromosome1.at(p) > 0)
			numOfDel++;
	}
	for(int p=0; p<theChromosome2.size() ; p++) {
		if(theChromosome2.at(p) > 0)
			numOfDel++;
	}

	numOfDelHetero = numOfDel - (2*numOfDelHomo);
	numOfBenHetero = numOfHetero - numOfDelHetero;

	this->heteroDeleterious = numOfDelHetero;
	this->heteroBeneficial = numOfBenHetero;
	this->homoBeneficial = numOfBenHomo;
	this->homoDeleterious = numOfDelHomo;
}



/* This functions returns the number of mutation on both chromosomes between loci
 determined by the arguments */
 // NOTE: the role of this function changed from what it was supposed to do. So, a difference between the name and function
int indiv::copyMutationsToGamete(vector<int>& newChromosome, vector<int> chromosome, int startLocus, int finishLocus, int whereToStart) {

	// this value will be returned so that the next times we don't search through the whole chromosome again
	int i;
	// loop through the chromosome and find the mutations in the range asked
	for(i=whereToStart; i < chromosome.size(); i++) {
		// if we have passed the range, don't search anymore
		if(abs(chromosome.at(i)) > finishLocus)
			break;
		// if the mutation is in the range add 1 to the returning value
		if(abs(chromosome.at(i)) < finishLocus && abs(chromosome.at(i)) >= startLocus)
			newChromosome.push_back(chromosome.at(i));
	}
	return i;
}



// This function performs a recombination between the chromosomes and returns a gamete to be inherited by the child
vector<int>* indiv::makeGamete(double meanPopulationRecombination, int isasex, int asexc) {

	// This is the gamete that will be returned
	vector<int>* gamete;
	gamete = new vector<int>;

	int numOfRecombination;
	vector<int> recombinationLoci;

    // these integers determine which alleles are inherited (for selfing, recombination, neutral gene, immunity A and B)
	int indexOfInheritedS, indexOfInheritedR, indexOfInheritedN, indexOfInheritedImmuneA, indexOfInheritedImmuneB;

	// These are GSL classes variables used to make
	// random number from a poisson distribution
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_default_seed += time(NULL);

	// Initialize the random number generator
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, gsl_rng_default_seed);

	// get the number of recombinations from a poisson ditrib. with mean=avg of recombination alleles (indirect selection on R exists)
	if(ONLY_DIRECT_SELECTION_ON_RECOMBINATION == 0 && isasex != 1) {
		numOfRecombination = gsl_ran_poisson(r,this->calcRecombinationRate()>0 ? this->calcRecombinationRate() : 0);
	}
	else if(ONLY_DIRECT_SELECTION_ON_RECOMBINATION != 0 && isasex != 1) {
		// if there is only direct selection on recombination, then the mean value for poisson distribution is independent
        // of the individual's recombination alleles and it equals mean recombination of population
		numOfRecombination = gsl_ran_poisson(r, meanPopulationRecombination);
	}else if(isasex == 1){
		// MH: With asex, no recombination occurs
		numOfRecombination = 0;
	}
	/*printf("With asex flag %d, no. of recs is %d\n",isasex,numOfRecombination);*/


	// get the random recombination loci
	double tempLocus;

	for(int i=0; i<numOfRecombination; i++) {
		gsl_rng_default_seed += time(NULL)^(i);
		gsl_rng_set (r, gsl_rng_default_seed);
		tempLocus = gsl_ran_flat(r, 1.5, L+0.4999);
		// check that the new recombination locus is not already in the vector
		int k;
		for(k=0 ; k<recombinationLoci.size() ; k++) {
			if(recombinationLoci.at(k) == floor(tempLocus + 0.5)) {
				i--;
				break;
			}
		}

		// If the recombination site is new, add it to vector
		if(k == recombinationLoci.size())
			recombinationLoci.push_back(floor(tempLocus + 0.5));
	}

	// sort the recombination loci
	sort(recombinationLoci.begin(), recombinationLoci.end());


	// add two new elements to beginning and at the end of the recombination loci vector (this is needed only for our algorithm)
	recombinationLoci.insert(recombinationLoci.begin(), 0);
	recombinationLoci.push_back(L+1);

	// renew the seed for more randomization
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

	// determine which chromosome to start with in recombination process
	// 0 means chromosome1; 1 means chromosome2
	int chromosomeToStart = 0;
	if(isasex != 1){
		chromosomeToStart = floor(gsl_ran_flat(r, -0.5, 1.4999) + 0.5);
	} else {
		chromosomeToStart = asexc;
	}
	/*printf("With asex flag %d, starting chromosome is %d\n",isasex,chromosomeToStart);*/
	
	// declare two numbers to keep track how much we have moved on each chromosome
	int placeOne=0, placeTwo=0;

	// start with chromosome1 if the value of "chromosomeToStart" is 0
	if(chromosomeToStart == 0) {

		// start making gamete by copying different parts from each chromosome
		for(int i=0; i<recombinationLoci.size()-1;) {
		    // copy mutations from chromosome1 to the gamete
			placeOne = this->copyMutationsToGamete(*gamete, this->theChromosome1, recombinationLoci.at(i), recombinationLoci.at(i+1), placeOne);
			// Check if R, S, N, immuneA & immuneB loci are going to be inherited from chromosome1 after recombination
			if(recombinationLoci.at(i)<=((double)L/3.0) && recombinationLoci.at(i+1)>((double)L/3.0))
				indexOfInheritedR = 0;
			if(recombinationLoci.at(i)<=((double)(2*L)/3.0) && recombinationLoci.at(i+1)>((double)(2*L)/3.0))
				indexOfInheritedS = 0;
			if(recombinationLoci.at(i)<=((double)L/2.0) && recombinationLoci.at(i+1)>((double)L/2.0))
				indexOfInheritedN = 0;
			if(recombinationLoci.at(i)<=((double)L/6.0) && recombinationLoci.at(i+1)>((double)L/6.0))
				indexOfInheritedImmuneA = 0;
			if(recombinationLoci.at(i)<=((double)(5*L)/6.0) && recombinationLoci.at(i+1)>((double)(5*L)/6.0))
				indexOfInheritedImmuneB = 0;

			i++;
			if(i<recombinationLoci.size()-1) {
			    // copy mutations from chromosome2 to the gamete
				placeTwo = this->copyMutationsToGamete(*gamete, this->theChromosome2, recombinationLoci.at(i), recombinationLoci.at(i+1), placeTwo);
				// Check if R, S, N, immuneA & immuneB loci are going to be inherited from chromosome2 after recombination
				if(recombinationLoci.at(i)<=((double)L/3.0) && recombinationLoci.at(i+1)>((double)L/3.0))
					indexOfInheritedR = 1;
				if(recombinationLoci.at(i)<=((double)(2*L)/3.0) && recombinationLoci.at(i+1)>((double)(2*L)/3.0))
					indexOfInheritedS = 1;
				if(recombinationLoci.at(i)<=((double)L/2.0) && recombinationLoci.at(i+1)>((double)L/2.0))
					indexOfInheritedN = 1;
				if(recombinationLoci.at(i)<=((double)L/6.0) && recombinationLoci.at(i+1)>((double)L/6.0))
					indexOfInheritedImmuneA = 1;
				if(recombinationLoci.at(i)<=((double)(5*L)/6.0) && recombinationLoci.at(i+1)>((double)(5*L)/6.0))
					indexOfInheritedImmuneB = 1;

				i++;
			}
		}
	}

	// here we start with chromosome2 because the value is 1
	else if(chromosomeToStart == 1) {

		// start making gamete by copying different parts from each chromosome
		for(int i=0; i<recombinationLoci.size()-1;) {
		    // copy mutations from chromosome2 to the gamete
			placeTwo = this->copyMutationsToGamete(*gamete, this->theChromosome2, recombinationLoci.at(i), recombinationLoci.at(i+1), placeTwo);
			// Check if R, S, N, immuneA & immuneB loci are going to be inherited from chromosome2 after recombination
			if(recombinationLoci.at(i)<=((double)L/3.0) && recombinationLoci.at(i+1)>((double)L/3.0))
				indexOfInheritedR = 1;
			if(recombinationLoci.at(i)<=((double)(2*L)/3.0) && recombinationLoci.at(i+1)>((double)(2*L)/3.0))
				indexOfInheritedS = 1;
			if(recombinationLoci.at(i)<=((double)L/2.0) && recombinationLoci.at(i+1)>((double)L/2.0))
				indexOfInheritedN = 1;
			if(recombinationLoci.at(i)<=((double)L/6.0) && recombinationLoci.at(i+1)>((double)L/6.0))
				indexOfInheritedImmuneA = 1;
			if(recombinationLoci.at(i)<=((double)(5*L)/6.0) && recombinationLoci.at(i+1)>((double)(5*L)/6.0))
				indexOfInheritedImmuneB = 1;

			i++;
			if(i<recombinationLoci.size()-1) {
			    // copy mutations from chromosome1 to the gamete
				placeOne = this->copyMutationsToGamete(*gamete, this->theChromosome1, recombinationLoci.at(i), recombinationLoci.at(i+1), placeOne);
				// Check if R, S, N, immuneA & immuneB loci are going to be inherited from chromosome1 after recombination
				if(recombinationLoci.at(i)<=((double)L/3.0) && recombinationLoci.at(i+1)>((double)L/3.0))
					indexOfInheritedR = 0;
				if(recombinationLoci.at(i)<=((double)(2*L)/3.0) && recombinationLoci.at(i+1)>((double)(2*L)/3.0))
					indexOfInheritedS = 0;
				if(recombinationLoci.at(i)<=((double)L/2.0) && recombinationLoci.at(i+1)>((double)L/2.0))
					indexOfInheritedN = 0;
				if(recombinationLoci.at(i)<=((double)L/6.0) && recombinationLoci.at(i+1)>((double)L/6.0))
					indexOfInheritedImmuneA = 0;
				if(recombinationLoci.at(i)<=((double)(5*L)/6.0) && recombinationLoci.at(i+1)>((double)(5*L)/6.0))
					indexOfInheritedImmuneB = 0;

				i++;
			}
		}
	}

    // free up the memory taken by the random number generator
	gsl_rng_free(r);

	gamete->push_back(indexOfInheritedN);
	gamete->push_back(indexOfInheritedS);
	gamete->push_back(indexOfInheritedR);
	gamete->push_back(indexOfInheritedImmuneA);
	gamete->push_back(indexOfInheritedImmuneB);

	return gamete;
}


// return the fitness of the individual
long double indiv::get_fitness() {
	return fitness;
}

// return chromosome 1 or 2 of the individual
vector<int>* indiv::get_theChromosome(int chromosomeNumber) {
	if (chromosomeNumber == 1)
		return &(theChromosome1);
	else if (chromosomeNumber == 2)
		return &(theChromosome2);
}
