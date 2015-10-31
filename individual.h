// This class represents individuals in the model
#include "included-files.h"

// This is just a small declaration to avoid errors from compiler
class studyPop;
class theWorld;



// Declaration of the class "indiv"
class indiv {

private:
	// vectors to save mutations on each chromosomes
	vector<int> theChromosome1;
	vector<int> theChromosome2;
	// floats to store some values for each individual
	float selfingGene[2];
	float recombinationGene[2];
	float neutralGene[2];
	int immunityGeneA[2];
	int immunityGeneB[2];
	long double fitness;

	int heteroDeleterious;
	int homoDeleterious;
	int heteroBeneficial;
	int homoBeneficial;

public:
	// Claculate the selfing rate of the individual
	float calcSelfingRate();
	// Calculate the number of homozygous and heterozygous del and ben mutations in the individual
	void calcFitness();
	// Calculate the recombination rate for individual's genome
	float calcRecombinationRate();
	// Search for the mutations between startLocus & finishLocus and rwturn the number of them
	int copyMutationsToGamete(vector<int>& newChromosome, vector<int> chromosome, int startLocus, int finishLocus, int whereToStart);
	// This function performs recombination between the 2 chromosomes and returns the chromosome that is inherited by child
	vector<int>* makeGamete(double meanPopulationRecombination, int isasex, int asexc);
	// The following functions are used to get info
	long double get_fitness();
	vector<int>* get_theChromosome(int chromosomeNumber);

	// Constructors
	indiv(int currentRound);
	indiv(vector<int> chromosome1, vector<int> chromosome2, double selfing1, double selfing2, double recombination1, double recombination2, double neutral1, double neutral2, int A1, int A2, int B1, int B2);

	// The studyPop class can access all the information of an individual
	friend class studyPop;
	friend class theWorld;
};
