// This class represents a single population with elements of type "indiv"

#include "individual.h"

class theWorld;



// Declaration of the studyPop class
class studyPop {
private:
	// Number of individuals in this studyPop
	int population;
	// The main vector to save population info in each generation
	vector<indiv> thePop;
	// This array contains a few genotypes of a parasite and stores the frequency of each genotype
	double parasites[9];
	// This array stores the genotypes of indivs at the immunnity loci
	double immunityGenotypeFreq[9];
	double tempHostFreq[9]; // This is for temp use
	// This table (2-D array) holds the reduction in fitness of the host according to its genotype, in the presence of parasites
	double hostFitnessTable[9][9];
	// this table (2-D array) determines the fitness of each gentype of parasite based on its genotype
	double parasiteFitnessTable[9][9];
	// A temp vector to save info about new indivs when making a new generation
	vector<indiv> newGen;
	// This value shows the highest fitness among the individuals in our population
	float highestW;
	// This value shows the highest ability among the individuals to outcross (this is used for fitness check for dad)
	float highestDadW;
	// This value shows the highest fitness among the individuals in the new temp generation
	float highestNewW;
	// This value shows the highest ability among the individuals to outcross in our new temp generation (used for fitness check for dad)
	float highestNewDadW;
	// This array saves the number of individuals made by selfing and also the total fitness of these individuals
	// First member is the # of indivs and second is the total fitness
	double selfedChildren[2];
	// This array saves the number of individuals made by outcrossing and also the total fitness of these individuals
	// First member is the # of indivs and second is the total fitness
	double outcrossedChildren[2];

	double meanRecombination;
	double newMeanRecombination;


public:
	// This function is called to make a new generation with the same population from the existing one
	void makeNewGeneration(int currentGen, int currentRound);
	// This functions picks a parent as mom and performs "fitness check" to determine if it can be a parent
	// to the child we want to add to population. It returns the index of a picked parent in the "thePop" vector.
	int pickAMom();
	// This function is called only when the mom outcrosses to reproduce. It picks an indiv from
	// studyPop and checks if it can be a father. It then returns the index of dad in "thePop" vectror.
	// It also checks to avoid picking the chosen mom as the dad again
	int pickADad(int chosenMom);
	// This functions make a new "indiv" using other functions
	// 1- picks a mom        		2-selfing or outcrossing
	// 3- picks a dad if outcrossing and gets a gamete from each parent
	// 4- add mutations to picked chromosomes		5- make a new indiv and add it to "newGen"
	void makeIndiv(int currentGen, int currentRound);
	// This function makes an indiv knowing mom and dad of the child
	indiv* makeForcedIndiv(int mom, int dad, int currentRound);
	// This function adds the effects of parasites on population in case they are included in the model
	void addParasitesEffect(indiv& newIndiv);
	// This function makes a new generation of parasites based on the method on Agrawal 2009
	void makeNewParasiteGeneration();
	// This function adds mutations to immunity genes in the host population
	void addMutationToImmunityGenes(int whichLocus);
	// This function add mutations to chromosomes that are picked from parents
	void addMutation(vector<int>& parent1, vector<int>& parent2, int currentRound);
	// The following three functions return info about the studyPop contents
	int get_population();
	vector<indiv>* get_thePop();
	vector<indiv>* get_newGen();

	// Constructor: makes the first generation
	studyPop(int currentRound);
	studyPop(vector<indiv> initPop);

	friend class theWorld;
};
