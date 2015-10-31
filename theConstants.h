/* this files defines the values for the constants.
by changing the values here you change them throughout the entire
code 

NOTE THOSE SURROUNDED BY ASTERISKS ARE NEW CONSTANTS INTRODUCED BY MH
TO MODEL ASEXUAL REPRODUCTION

*/

/***************** List Of Constants Used ************************
*                                                                *
*  (1) L -->  means chromosome length (a positive integer)       *
*  (2) STUDYPOPNO --> # of separate populations in the model     *
*  (3) POPULATION --> # of individuals in each separate studyPop *
*  (4) INITIALSELFING --> the selfing rate of an individual with *
*      no mutations in S loci                                    *
*  (5) INITIALRECOMBINATION --> recombination rate for an        *
*      inidividual with no mutations in M loci                   *
*  (6) INITIALNEUTRALGENE --> initial value of neutral gene for  *
*      individuals in the first generation                       *
*  (7) NEUTRALMUTATIONRATE --> mutation rate at neutral locus    *
*  (8) RECOMBMUTATIONRATE --> mutation rate at recomb. locus     *
*  (9) SELFINGMUTATIONRATE --> mutation rate at slefing locus    *
* (10) MUTATIONRATE --> mutations for each diploid genome        *
*      during each division                                      *
* (11) PROB_BEN --> the ratio of the mutations that are          *
*      beneficial                                                *
* (12) SIGMAN --> sd of the normal distribution from which       *
*      neutral gene mutation values are randomely picked         *
* (13) SIGMAS --> sd of the normal distribution from which       *
*      selfing gene mutation values are randomely picked         *
* (14) SIGMAR --> sd of the normal distribution from which       *
*      recombination gene mutation values are randomely picked   *
* (15) S --> this constant shows the fitness cost of having      *
*      a homozygous mutation at a locus                          *
* (16) BEN_S --> selection for beneficial alleles                *
* (17) H --> this is a constant to determine the relative cost   *
*      of heterozygous mutation to homozygous mutation           *
* (18) H_BEN --> this constant determines the relative benefit   *
*      of homozygous to heterozygous ben mutation                *
* (19) K --> pollen discounting: determines by how much, selfing *
*      reduces the chance for outcrossing                        *
* (20) RECOMBCOST --> the cost of higher recombination. This     *
*      constant puts direct selection on receombination if the   *
*      value is larger than 0                                    *
* (21) *ASEXP* --> probability than an individual clonally 		 *
* 	   reproduces (as opposed to outcross)						 *
*****************************************************************/


#define L 5000
#define STUDYPOPNO 1
#define POPULATION 10000
#define INITIALSELFING 0.0
#define INITIALRECOMBINATION 1.0
#define	INITIALNEUTRALGENE 0.0
#define NEUTRALMUTATIONRATE 1.0
#define RECOMBMUTATIONRATE 0.0
#define SELFINGMUTATIONRATE 0.0
#define MUTATIONRATE 0.02
#define PROB_BEN 0.0
#define SIGMAN 1
#define SIGMAS 0
#define SIGMAR 0
#define S 0.02
#define BEN_S 0.0
#define H 0.25
#define H_BEN 0.0
#define K 0.0
#define RECOMBCOST 0.0
#define ASEXP 0.0

// These two constants determine the desired number of replicates and
// the desired number of generations in each replicate
// NOTE: the number og generations may increase if the simulation does not arrive at an equilibrium
#define NUMOFROUNDS 1
#define NUMOFGENERATIONS 40000

// These constants determine whether the simulation should save the last
// generation, perform regression analysis and/or find LD at the last generation
#define SAVE_LAST_GENERATION 1
#define REGRESSION_ANALYSIS 0
#define FIND_LD 0


// Constants related to parasites (names are taken from Agrawal 2009)
#define PARASITES_EXIST 0
#define U_HOST_IMMUNITY 0.000001
#define U_PARASITE 0.000001
#define MODEL_OF_PARASITE_INFECTION "GFG"
#define PARASITE_PLOIDY 1
#define PLEIOTROPIC_COST_FOR_HOST 0.1
#define PLEIOTROPIC_COST_FOR_PARASITE 0.3
#define D1 0.5
#define D2 1.0
#define V 0.3
#define T_COST 1.0
#define EFFECTIVE_PARASITE_R 0.1


/**** Misc constants ****/
#define ONLY_DIRECT_SELECTION_ON_RECOMBINATION 0
#define INITIAL_SELFING_VARIATION 0
#define EXPERIMENTAL_INBREEDING_DEP 0
// these constant determine the information in the name of the text file that saves the data
#define R_IN_NAME 1
#define K_IN_NAME 0
#define P_BEN_IN_NAME 0
#define S_BEN_IN_NAME 0
#define ASEX_IN_NAME 1
