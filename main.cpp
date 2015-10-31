#include "theWorld.h"




int main() {

    // perform simulations until the desired number of replicates is reached
	for(int currentRound = 1; currentRound<=NUMOFROUNDS ; currentRound++) {
        // create an object that initializes and controls all the populations and individuals
        theWorld modelWorld(1, currentRound);

        // run the simulation until the desired number of generations is reached
        modelWorld.runModel(NUMOFGENERATIONS, currentRound);
	}

	return 0;
}
