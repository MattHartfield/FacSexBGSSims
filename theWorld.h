// This class contains the collection of all populations in the model
// and controls for how many generations the model should run

#include "studyPop.h"



// Declaration of the class "theWorld"
class theWorld {
private:
    // collection of all populations (we had only one in our simulations)
	vector<studyPop> popCollection;
	// number of current generation for each population
	vector<int> generationNo;
	// vector to save a random sample of indivs (it is only used in our calculation for inbreeding dep)
	vector<int> sampledIndivs;
	// vectors to save children that are produced by forced inbreeding or outcrossing
	// These are used every 20 generations to calculate ID
	vector<indiv> forcedInbred;
	vector<indiv> forcedOutbred;

	int numOfBenFixed;
	int numOfDelFixed;


public:

	double S_ben;
	double P_ben;

	// When this function is called it runs the model several times according to the parameters
	// "generations" determine how many generations we want this model to run
	void runModel(long int generations, int whichRounds);
	// This function helps in finding some of the averages in the population. whichPop gives the index of the population
	// under syudy and whichAvg determines which average value we are looking for.
	// whichAvg: 1 means avg del alleles ; 2 means avg ben alleles ; 3 means avg homo del alleles ; 4 means avg homo ben alleles
	// 5 means avg fitness of population ; 6 means avg selfing rate of population
	double findAvg(int whichAvg, int whichPop);
	// This function gets a random sample of individuals in a population and stores the sample in sampledIndivs
	void takeSample(int whichPop, int sampleSize);
	// This function clears the fixed mutations in the population (it is called every 20 generations)
	void clearFixedMutations(int whichPopulation);
	// This function stores information from all the individuals in the generation when it is called
	void saveThisGeneration(int whichRound);
	// This function performs some regression analysis and finds covariences btwn SelfingRate and: w, X1, X2, e
	vector<double>* regressionAnalyzer(int whichPop, int whichRound);
	vector<double>* regressionAnalyzer_noBeneficial(int whichPop, int whichRound);
	// This function finds the linkage disequilibrium between the loci
	void findLD(long double& avgFreqOfDelAlleles, long double& avgPolymorphOfDelAlleles, double& avgDij, long double& avgD_prime);
	// This function finds linkage disequilibrium between two loci
	double findLD_betweenTwo(int firstLocus, int secondLocus, vector<double> delFrequencies);
	// This class find the frequencies of mutaions over the alleles
	// If which type is 1, then del mutations, if which type is 2 ben mutations are considered
	void findFrequencies(vector<double>& freqArray, int whichType);
	// This function return 1 if the population selfing rate is at equilibrium and 0 otherwise
	int isAtEquil(vector<double> selfingRates);
	// The following two functions return info about the contents of "theWorld"
	vector<studyPop>* get_popCollection();
	vector<int>* get_generationNo();

	// Constructor: This constructs the "studyPop"s
	theWorld(int pops, int currentRound);

};
