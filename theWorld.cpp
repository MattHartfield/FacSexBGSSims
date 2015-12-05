#include "theWorld.h"



// The constructer of this class makes all the populations, each of which stores its first generation
// NOTE: in our simulations, we only have one population (i.e. the "pops" argument is set to 1)
theWorld::theWorld(int pops, int currentRound){

	for(int i=0; i<pops; i++) {
		studyPop newPop(currentRound);
		popCollection.push_back(newPop);
		generationNo.push_back(1);
	}

    // set values related to beneficial mutations, based on the constants
    
	S_ben = BEN_S;
    P_ben = PROB_BEN;

    
	// calculate the fitness for every individual in the population
    
	for(int i=0;i<POPULATION; i++) {
        popCollection.at(0).thePop.at(i).fitness = (pow((1 + S_ben),2*popCollection.at(0).thePop.at(i).homoBeneficial+H_BEN*popCollection.at(0).thePop.at(i).heteroBeneficial))*(pow((1 - S),popCollection.at(0).thePop.at(i).homoDeleterious)) * (pow((1 - H*S),popCollection.at(0).thePop.at(i).heteroDeleterious));
    }

    // find he highest fitness in the population
    popCollection.at(0).highestW = popCollection.at(0).thePop.at(0).get_fitness();
    for(int i=1; i<POPULATION ; i++) {
        if(popCollection.at(0).thePop.at(i).get_fitness() > popCollection.at(0).highestW)
            popCollection.at(0).highestW = popCollection.at(0).thePop.at(i).get_fitness();
    }

    // find the highest ability to involve in an outcrossing reproduction
    popCollection.at(0).highestDadW = (1 - K*popCollection.at(0).thePop.at(0).calcSelfingRate())*popCollection.at(0).thePop.at(0).get_fitness();
    for(int i=1; i<POPULATION ; i++) {
        if((1 - K*popCollection.at(0).thePop.at(i).calcSelfingRate())*popCollection.at(0).thePop.at(i).get_fitness() > popCollection.at(0).highestDadW)
            popCollection.at(0).highestDadW = (1 - K*popCollection.at(0).thePop.at(i).calcSelfingRate())*popCollection.at(0).thePop.at(i).get_fitness();
    }

	numOfBenFixed = 0;
	numOfDelFixed = 0;
}



// This is a helper function. It is used in other functions below to find some average values
// Different values of "whichAvg" argument correspond to different means
// whichAvg==1 --> find mean del mutations per individual, whichAvg==2 --> find mean ben mutations per individual
// whichAvg==3 --> find mean homozygous del alleles per individual, whichAvg==4 --> find mean homozygous ben alleles per individual
// whichAvg==5 --> find mean fitness, whichAvg==6 --> find mean selfing rate in the population
double theWorld::findAvg(int whichAvg, int whichPop) {
	double returningVal;

	if(whichAvg == 1) {
		int totalDel=0;
		// find average of deleterious alleles per indiv
		for(int i=0; i<POPULATION ; i++) {
			totalDel += this->popCollection.at(whichPop).thePop.at(i).heteroDeleterious + 2*this->popCollection.at(whichPop).thePop.at(i).homoDeleterious;
		}
		returningVal = ((double)totalDel)/((double)POPULATION);
	}

	else if(whichAvg == 2) {
		int totalBen=0;
		// find average of beneficial alleles per indiv
		for(int i=0; i<POPULATION ; i++) {
			totalBen += this->popCollection.at(whichPop).thePop.at(i).heteroBeneficial + 2*this->popCollection.at(whichPop).thePop.at(i).homoBeneficial;
		}
		returningVal = ((double)totalBen)/((double)POPULATION);
	}

	else if(whichAvg == 3) {
		int totalHomoDel=0;
		// find average of homozygote deleterious alleles
		for(int i=0; i<POPULATION ; i++) {
			totalHomoDel += this->popCollection.at(whichPop).thePop.at(i).homoDeleterious;
		}
		returningVal = ((double)totalHomoDel)/((double)POPULATION);
	}

	else if(whichAvg == 4) {
		int totalHomoBen=0;
		// find average of homozygote beneficial alleles
		for(int i=0; i<POPULATION ; i++) {
			totalHomoBen += this->popCollection.at(whichPop).thePop.at(i).homoBeneficial;
		}
		returningVal = ((double)totalHomoBen)/((double)POPULATION);
	}

	else if(whichAvg == 5) {
		double totalFitness=0;
		// find average fitness in the population
		for(int i=0; i<POPULATION ; i++) {
			totalFitness += this->popCollection.at(whichPop).thePop.at(i).get_fitness();
		}
		returningVal = totalFitness/((double)POPULATION);
	}

	else if(whichAvg == 6) {
		double totalSelfing=0;
		// find average selfing rate in the population
		for(int i=0; i<POPULATION ; i++) {
			totalSelfing += this->popCollection.at(whichPop).thePop.at(i).calcSelfingRate();
		}
		returningVal = totalSelfing/((double)POPULATION);
	}

	return returningVal;
}





// This function performs the simulations and saves the data
void theWorld::runModel(long int generations, int whichRounds) {

    // create a text file to save data into it (every 20 generations)
	ofstream modelResults;

    // create a string for the name of the text file
	stringstream sstr;
	sstr << whichRounds;

	stringstream sstr2;
	sstr2 << popCollection.at(0).thePop.at(0).calcRecombinationRate();
	stringstream sstr1;
	sstr1 << K;
	stringstream sstr3;
	sstr3 << P_ben;
	stringstream sstr4;
	sstr4 << S_ben;
	stringstream sstr5;
	sstr5 << ASEXP;

	string fileName("model results");
	fileName.append(sstr.str());

	fileName.append(" (");
	if(R_IN_NAME) {
		fileName.append("R=");
		fileName.append(sstr2.str());
	}
	if(K_IN_NAME) {
		fileName.append(",K=");
		fileName.append(sstr1.str());
	}
	if(P_BEN_IN_NAME) {
		fileName.append(",P_ben=");
		fileName.append(sstr3.str());
	}
	if(S_BEN_IN_NAME) {
		fileName.append(",S_ben=");
		fileName.append(sstr4.str());
	}
	if(ASEX_IN_NAME) {
		fileName.append(",Asp=");
		fileName.append(sstr5.str());
	}

	fileName.append(")");

	fileName.append(".txt");

	
	
	// open the text file in which the data will be saved
	modelResults.open(fileName.c_str(), ios::out);
	
	/*

    // create a text file in which regression analysis results will be saved (every 20 generations)
	ofstream regressionAnalyzer;

    // create a name for the file in which regression results will be saved
	string fileName2("Regression_Analysis");
	fileName2.append(sstr.str());

	fileName2.append(" (");
	if(R_IN_NAME) {
		fileName2.append("R=");
		fileName2.append(sstr2.str());
	}
	if(K_IN_NAME) {
		fileName2.append(",K=");
		fileName2.append(sstr1.str());
	}
	if(P_BEN_IN_NAME) {
		fileName2.append(",P_ben=");
		fileName2.append(sstr3.str());
	}
	if(S_BEN_IN_NAME) {
		fileName2.append(",S_ben=");
		fileName2.append(sstr4.str());
	}
	if(ASEX_IN_NAME) {
		fileName.append(",Asp=");
		fileName.append(sstr5.str());
	}


	fileName2.append(")");
	fileName2.append(".txt");

	regressionAnalyzer.open(fileName2.c_str(), ios::out);
	regressionAnalyzer << "Cov[SelfingRate & n_hom_del],Cov[SelfingRate & n_hom_ben],Cov[SelfingRate & n_het_del],Cov[SelfingRate & n_het_ben],Cov[SelfingRate & w],Cov[SelfingRate & X1],Cov[SelfingRate & X2],Cov[SelfingRate & e],B1,B2,B3,B4" << '\n';
	*/

	// Save the constants of the simulation into the text file created above
	modelResults << "Chromosome Length," << L << "\n";
	modelResults << "Population," << POPULATION << "\n";
	modelResults << "Initial Selfing," << INITIALSELFING << '\n';
	modelResults << "Initial Recombination," << popCollection.at(0).thePop.at(0).calcRecombinationRate() << "\n";
	modelResults << "Initial Neutral Value," << INITIALNEUTRALGENE << "\n";
	modelResults << "U neutral," << NEUTRALMUTATIONRATE << "\n";
	modelResults << "U recombination," << RECOMBMUTATIONRATE << "\n";
	modelResults << "U selfing," << SELFINGMUTATIONRATE << "\n";
	modelResults << "U deleterious," << MUTATIONRATE << "\n";
	modelResults << "Beneficial Mutation Probability," << P_ben << "\n";
	modelResults << "SIGMA (for neutral mutations)," << SIGMAN << "\n";
	modelResults << "SIGMA (for selfing mutations)," << SIGMAS << "\n";
	modelResults << "SIGMA (for recombination mutations)," << SIGMAR << "\n";
	modelResults << "Selection Coefficient," << S << "\n";
	modelResults << "S_beneficial," << S_ben << "\n";
	modelResults << "H," << H << "\n";
	modelResults << "H_ben," << H_BEN << "\n";
	modelResults << "K," << K << "\n";
	modelResults << "Cost Of Recombination," << RECOMBCOST << "\n";
	modelResults << "Probability of asexual reproduction," << ASEXP << "\n";	
	modelResults << "Paraqsites exist?," << (PARASITES_EXIST==1 ? "YES" : "NO") << "\n";
	modelResults << "U host immunity," << U_HOST_IMMUNITY << "\n";
	modelResults << "U parasite," << U_PARASITE << "\n";
	modelResults << "Model for host-parasite interaction," << MODEL_OF_PARASITE_INFECTION << "\n";
	modelResults << "Parasite ploidy," << PARASITE_PLOIDY << "\n";
	modelResults << "Pleiotropic cost for host," << PLEIOTROPIC_COST_FOR_HOST << "\n";
	modelResults << "Pleiotropic cost for parasite," << PLEIOTROPIC_COST_FOR_PARASITE << "\n";
	modelResults << "Probability of heterozygote host being universally matched," << D1 << "\n";
	modelResults << "Probability that A/a behaves like A/A," << D2 << "\n";
	modelResults << "Parasite virulence," << V << "\n";
	modelResults << "Selection against unsuccessful parasites," << T_COST << "\n";
	modelResults << "Effective parasite recombination rate," << EFFECTIVE_PARASITE_R << "\n";
	modelResults << "\n\n";
	modelResults << "Generation #,Mean Del mutations per indiv,Mean Ben mutations per indiv,Variance of Del mutations/indiv,Variance of Ben mutations/indiv,Mean fitness,Mean Neutral Gene Values,Variance of neutral gene values,Mean Selfing rate,Variance of selfing rates,Mean Recombination Rate,Variance for recombination rate,Inbreeding Depression,Freq Of 'a' Allele (Host),Freq Of 'b' Allele (Host),Freq Of 'a' Allele (Parasite),Freq Of 'b' Allele (Parasite),# Fixed Deleterious Mutations,# Fixed Beneficial Mutations,avg Freq of Del Alleles,avg Polymorphism of Del Alleles,avg disequilibria,avg D_prime" << "\n";

	// A variable to keep track of time
	time_t startTime = time(0);

	vector<double> selfingRates;

    // This variable will show how many generatins the simulation should run. It will increase by 5000 if mean selfing rate does not reach equilibrium
	int generationsToRun = generations;

	vector<double>* regressionResults;

	// make new generations and save their data
	for(long int i=1; i<generationsToRun; i++) {
		for(int j=0; j<popCollection.size(); j++) {

			// save data every 20 generations
			// MH: Every 100 in my case, and ONLY after 5N gens
			if((((i-1)%(5*POPULATION/200) == 0) && ((i-1) >= (5*POPULATION)))) {

				// calculate mean values & variances
				vector<int> Delmutations;
				vector<int> Benmutations;
				vector<double> neutralGeneValues;
				vector<double> selfingRateValues;
				vector<double> recombinationValues;
				double meanW=0;
				double meanRelativeW=0;
				// save values of: # of ben & del mutatoins, neutral gene values, selfing rates, fitness, relative fitness recombination rate
				for(int k=0; k<POPULATION ; k++) {
					Delmutations.push_back(popCollection.at(j).thePop.at(k).heteroDeleterious + 2*popCollection.at(j).thePop.at(k).homoDeleterious);
					Benmutations.push_back(popCollection.at(j).thePop.at(k).heteroBeneficial + 2*popCollection.at(j).thePop.at(k).homoBeneficial);
					neutralGeneValues.push_back((popCollection.at(j).thePop.at(k).neutralGene[0]+popCollection.at(j).thePop.at(k).neutralGene[1])/2);
					selfingRateValues.push_back(popCollection.at(j).thePop.at(k).calcSelfingRate());
					recombinationValues.push_back(popCollection.at(j).thePop.at(k).calcRecombinationRate());
					meanW += popCollection.at(j).thePop.at(k).fitness;
					meanRelativeW += popCollection.at(j).thePop.at(k).fitness/popCollection.at(j).highestW;
				}

				// calculate mean W and mean relative fitness
				meanW /= POPULATION;
				meanRelativeW /= POPULATION;

				// calculate total of all neutral gene values in the population & save squares of all values
				long double totalNeutralValues=0;
				long double totalNeutralValuesSquares=0;
				for(int k=0; k<neutralGeneValues.size() ; k++) {
					totalNeutralValues += neutralGeneValues.at(k);
					totalNeutralValuesSquares += pow(neutralGeneValues.at(k), 2.0);
				}
				// calculate mean neutral values & the variance for neutral genes
				double meanNeutralValues = totalNeutralValues/((double)neutralGeneValues.size());
				double varianceForNeutralValues = (totalNeutralValuesSquares)/((double)neutralGeneValues.size()) - pow(meanNeutralValues, 2.0);




				// calculate total of all selfing rates in the population & save squares of all rates
				long double totalSelfingRates=0;
				long double totalSelfingRatesSquares=0;
				for(int k=0; k<POPULATION; k++) {
					totalSelfingRates += selfingRateValues.at(k);
					totalSelfingRatesSquares += pow(selfingRateValues.at(k), 2.0);
				}
				// calculate mean selfing rate and variance for selfing rates
				double meanSelfingRates = totalSelfingRates/POPULATION;
				double varianceForSelfingRates = (totalSelfingRatesSquares)/((double)POPULATION) - pow(meanSelfingRates, 2.0);



				// calculate total of all recombination rates in the population & save squares of all rates
				long double totalRecombinationRates=0;
				long double totalRecombinationRatesSquares=0;
				for(int k=0; k<POPULATION ; k++) {
					totalRecombinationRates += recombinationValues.at(k);
					totalRecombinationRatesSquares += pow(recombinationValues.at(k), 2.0);
				}
				// calculate mean recombination rate and variance for recombination rates
				double meanRecombinationRates = totalRecombinationRates/((double)POPULATION);
				double varianceForRecombinationRates = totalRecombinationRatesSquares/((double)POPULATION) - pow(meanRecombinationRates, 2.0);




				// calculate total of all mutation counts in the population & save squares of all mutation counts
				long int totalDelMutations=0;
				long int totalDelSquares=0;
				long int totalBenMutations=0;
				long int totalBenSquares=0;
				for(int k=0; k<Delmutations.size() ; k++) {
					totalDelMutations += Delmutations.at(k);
					totalDelSquares += pow((double)Delmutations.at(k),2.0);
				}
				for(int k=0; k<Benmutations.size() ; k++) {
					totalBenMutations += Benmutations.at(k);
					totalBenSquares += pow((double)Benmutations.at(k),2.0);
				}
				// calculate mean mutation count and the variance for mutation counts
				double meanDelMutations = ((double)totalDelMutations)/((double)POPULATION);
				long double varianceForDelMutations = ((long double)totalDelSquares)/((long double)POPULATION) - (long double)pow(meanDelMutations, 2.0);
				double meanBenMutations = ((double)totalBenMutations)/((double)POPULATION);
				long double varianceForBenMutations = ((long double)totalBenSquares)/((long double)POPULATION) - (long double)pow(meanBenMutations, 2.0);


				// If selfing is zero or 1, or "EXPERIMENTAL_INBREEDING_DEP" is set to 1, do an experiment to make forced inbred and outbred children
				// This is used to calculate inbreeding depression
				// MH: I've altered these conditions since I don't want to perform this test by default in my simulations
				//if(popCollection.at(j).selfedChildren[0]==0 || popCollection.at(j).outcrossedChildren[0]==0 || EXPERIMENTAL_INBREEDING_DEP==1) {
				/*
				if(EXPERIMENTAL_INBREEDING_DEP==1) {

                    // take a random sample of 200 individuals from the population
                    // if population size<200, the sample size is half of the population size
					takeSample(j, 200 < POPULATION ? 200 : POPULATION/2);

					// Initialie random number generator
					const gsl_rng_type * T;
					gsl_rng * r;
					T = gsl_rng_default;
					r = gsl_rng_alloc (T);
					gsl_rng_default_seed += time(NULL);
					gsl_rng_set (r, gsl_rng_default_seed);


					// make forced inbred Children with 100 random indivs from our sample
					for(int forcedChild=0; forcedChild<100 ; forcedChild++) {
						int selferMom = floor(gsl_ran_flat(r, -0.5, (sampledIndivs.size()-1)+0.4999) + 0.5);
						indiv* newInbred = popCollection.at(j).makeForcedIndiv(sampledIndivs.at(selferMom), sampledIndivs.at(selferMom), whichRounds);
						newInbred->fitness = (pow((1 + S_ben),2*newInbred->homoBeneficial+newInbred->heteroBeneficial))*(pow((1 - S),newInbred->homoDeleterious)) * (pow((1 - H*S),newInbred->heteroDeleterious));
						this->forcedInbred.push_back(*newInbred);
						delete newInbred;
						gsl_rng_default_seed += time(NULL)*selferMom;
						gsl_rng_set (r, gsl_rng_default_seed);
					}


					// make forced outbred Children with randomely chosen parents from the sample
					for(int forcedChild=0; forcedChild<100 ; forcedChild++) {
						int mom = floor(gsl_ran_flat(r, -0.5, (sampledIndivs.size()-1)+0.4999) + 0.5);
						gsl_rng_default_seed += time(NULL)*mom;
						gsl_rng_set (r, gsl_rng_default_seed);
						int dad = floor(gsl_ran_flat(r, -0.5, (sampledIndivs.size()-1)+0.4999) + 0.5);
						indiv* newOutbred = popCollection.at(j).makeForcedIndiv(sampledIndivs.at(mom), sampledIndivs.at(dad), whichRounds);
						newOutbred->fitness = (pow((1 + S_ben),2*newOutbred->homoBeneficial+newOutbred->heteroBeneficial))*(pow((1 - S),newOutbred->homoDeleterious)) * (pow((1 - H*S),newOutbred->heteroDeleterious));
						this->forcedOutbred.push_back(*newOutbred);
						delete newOutbred;
						gsl_rng_default_seed += time(NULL)*dad;
						gsl_rng_set (r, gsl_rng_default_seed);
					}

                    // save the number of inbred and outbred children and the total fitness of each group (this will be used later to calculate inbreeding depression)
					popCollection.at(j).selfedChildren[0] = this->forcedInbred.size();
					popCollection.at(j).selfedChildren[1] = 0;
					for(int p=0; p<forcedInbred.size() ; p++) {
						popCollection.at(j).selfedChildren[1] += this->forcedInbred.at(p).get_fitness();
					}
					popCollection.at(j).outcrossedChildren[0] = this->forcedOutbred.size();
					popCollection.at(j).outcrossedChildren[1] = 0;
					for(int p=0; p<forcedOutbred.size() ; p++) {
						popCollection.at(j).outcrossedChildren[1] += this->forcedOutbred.at(p).get_fitness();
					}

					gsl_rng_free(r);

					sampledIndivs.clear();
				}
				*/

				// Calculate the frequency of immunity alleles both for hosts and parasites
				/*
				double freq_a_host=0, freq_b_host=0, freq_a_parasite=0, freq_b_parasite=0;
				for(int l=0; l<9 ; l++) {
					freq_a_host += ((double)(l/3)/2.0)*popCollection.at(j).immunityGenotypeFreq[l];
					freq_b_host += ( ((double)(l-(l/3)*3)) /2.0)*popCollection.at(j).immunityGenotypeFreq[l];
					if(PARASITE_PLOIDY == 1) {
						if(l<5 && l!=2) {
							freq_a_parasite += ((l/3))*popCollection.at(j).parasites[l];
							freq_b_parasite += (l-(l/3)*3)*popCollection.at(j).parasites[l];
						}
					}
				}


				// Save the mean selfing rates for the last 3000 generations (to check for equiliberium)
				if((generationsToRun-i+1) <= 3000) {
					selfingRates.push_back(meanSelfingRates);
				}
		
				*/		
				// if this is the last generation, check if we have reached an equilibrium (based on the selfing rates in the last 3000 generations)
				// if equiliberium is not reached, run for another 5000 generations and check for equilibrium (repeat this until equiliberium is reached).
				// MH: Wiping this since I do NOT want extensions
				/*
				if((generationsToRun-i == 1) && !(isAtEquil(selfingRates))) {
					generationsToRun += 5000;
					selfingRates.clear();
				}
				*/


				// Finding LD (Linkage Disequilibrium)
				
				long double avgFreqOfDelAlleles=0, avgPolymorphOfDelAlleles=0, avgD_prime=0;
				double avgDij=0;

				if(generationsToRun-i == 1 && FIND_LD == 1) {
					this->findLD(avgFreqOfDelAlleles, avgPolymorphOfDelAlleles, avgDij, avgD_prime);
				}
				/*
                // perform regression analysis
				if(P_ben == 0)
					regressionResults = this->regressionAnalyzer_noBeneficial(j, whichRounds);
				else
					regressionResults = this->regressionAnalyzer(j, whichRounds);

                // write the regression analysis results into a text file
				regressionAnalyzer << regressionResults->at(0) << ",";
				regressionAnalyzer << regressionResults->at(1) << ",";
				regressionAnalyzer << regressionResults->at(2) << ",";
				regressionAnalyzer << regressionResults->at(3) << ",";
				regressionAnalyzer << regressionResults->at(4) << ",";
				regressionAnalyzer << regressionResults->at(5) << ",";
				regressionAnalyzer << regressionResults->at(6) << ",";
				regressionAnalyzer << regressionResults->at(7) << ",";
				regressionAnalyzer << regressionResults->at(8) << ",";
				regressionAnalyzer << regressionResults->at(9) << ",";
				regressionAnalyzer << regressionResults->at(10) << ",";
				regressionAnalyzer << regressionResults->at(11) << "\n";

				regressionResults->clear();
				delete regressionResults;
				*/

				// print out info during simulation
				cout << generationNo.at(j) << ", " << meanDelMutations << ", " << varianceForDelMutations << ", " << meanBenMutations << ", " << meanW << ", " << meanNeutralValues << ", " << varianceForNeutralValues << ", " << meanSelfingRates << ", ";
				cout << varianceForSelfingRates << ", " << meanRecombinationRates << ", " << varianceForRecombinationRates << ", ";
				cout << 1- (popCollection.at(j).selfedChildren[1]/popCollection.at(j).selfedChildren[0])/(popCollection.at(j).outcrossedChildren[1]/popCollection.at(j).outcrossedChildren[0]) << ", ";
				/* cout << freq_a_host << ", " << freq_b_host << ", " << numOfDelFixed << ", " << numOfBenFixed << ", ";	*/
				cout << numOfDelFixed << ", " << numOfBenFixed << ", ";
				cout << avgFreqOfDelAlleles << ", " << avgPolymorphOfDelAlleles << ", " << avgDij << ", " << avgD_prime << "\n";

				/* cout << difftime(time(0), startTime) << endl;*/
				startTime = time(0);

				// Write info in a text file to save them
				modelResults << generationNo.at(j) << "," << meanDelMutations << "," << meanBenMutations << "," << varianceForDelMutations << "," << varianceForBenMutations << "," << meanW << "," << meanNeutralValues << "," << varianceForNeutralValues << ",";
				modelResults << meanSelfingRates << "," <<varianceForSelfingRates << "," << meanRecombinationRates << "," << varianceForRecombinationRates << ",";
				/* modelResults << 1- (popCollection.at(j).selfedChildren[1]/popCollection.at(j).selfedChildren[0])/(popCollection.at(j).outcrossedChildren[1]/popCollection.at(j).outcrossedChildren[0]) << "," << freq_a_host << "," << freq_b_host << ","; */
				modelResults << 1- (popCollection.at(j).selfedChildren[1]/popCollection.at(j).selfedChildren[0])/(popCollection.at(j).outcrossedChildren[1]/popCollection.at(j).outcrossedChildren[0]) << ","; 
				modelResults << numOfDelFixed << "," << numOfBenFixed << "," << avgFreqOfDelAlleles << "," << avgPolymorphOfDelAlleles << "," << avgDij << "," << avgD_prime << "\n";

                // clear some saved data and prepare for the next round of saving data
				forcedInbred.clear();
				forcedOutbred.clear();

				this->clearFixedMutations(j);

			}
			popCollection.at(j).selfedChildren[0] = 0.0;
			popCollection.at(j).selfedChildren[1] = 0.0;
			popCollection.at(j).outcrossedChildren[0] = 0.0;
			popCollection.at(j).outcrossedChildren[1] = 0.0;

			// save all the info in the last generation if the constant to save last generation is set to 1
			if(i == generationsToRun-1 && SAVE_LAST_GENERATION == 1) {
				this->saveThisGeneration(whichRounds);
			}

			// create a new generation for the study pop
			popCollection.at(j).makeNewGeneration(generationNo.at(j), whichRounds);

            // calculate the fitness of every individual in the new generation
			for(int n=0;n<POPULATION; n++) {
				popCollection.at(0).thePop.at(n).fitness = (pow((1 + S_ben),2*popCollection.at(0).thePop.at(n).homoBeneficial+H_BEN*popCollection.at(0).thePop.at(n).heteroBeneficial))*(pow((1 - S),popCollection.at(0).thePop.at(n).homoDeleterious)) * (pow((1 - H*S),popCollection.at(0).thePop.at(n).heteroDeleterious));
			}

            // find the highest fitness in the population (new generation)
			popCollection.at(0).highestW = popCollection.at(0).thePop.at(0).get_fitness();
			for(int n=1; n<POPULATION ; n++) {
				if(popCollection.at(0).thePop.at(n).get_fitness() > popCollection.at(0).highestW)
					popCollection.at(0).highestW = popCollection.at(0).thePop.at(n).get_fitness();
			}

            // find the highest ability to engage in an outcrossing (for the new generation)
			popCollection.at(0).highestDadW = (1 - K*popCollection.at(0).thePop.at(0).calcSelfingRate())*popCollection.at(0).thePop.at(0).get_fitness();
			for(int n=1; n<POPULATION ; n++) {
				if((1 - K*popCollection.at(0).thePop.at(n).calcSelfingRate())*popCollection.at(0).thePop.at(n).get_fitness() > popCollection.at(0).highestDadW)
					popCollection.at(0).highestDadW = (1 - K*popCollection.at(0).thePop.at(n).calcSelfingRate())*popCollection.at(0).thePop.at(n).get_fitness();
			}

			generationNo.at(j) = i+1;

		}
	}
	// close the text files which have the data and the regression analysis
	/*regressionAnalyzer.close();*/
	modelResults.close();
}



// This function takes a random sample from the population
void theWorld::takeSample(int whichPop, int sampleSize) {

	// Initialie random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

    // this loop keeps on picking individuals until we have the desired sample size
	for(int i=0; i<sampleSize; i++) {
		int tempIndiv;
		int pass;
		// pick a new random individual
		do {
			pass = 1;
			tempIndiv = floor(gsl_ran_flat(r, -0.5, (POPULATION-1)+0.4999) + 0.5);
			// check that the individual is new (not in the sample already
			for(int l=0; l<this->sampledIndivs.size() ; l++) {
				if(sampledIndivs.at(l) == tempIndiv)
					pass--;
			}

			gsl_rng_default_seed += time(NULL);
			gsl_rng_set (r, gsl_rng_default_seed);
		} while(pass<1);

		this->sampledIndivs.push_back(tempIndiv);
	}
    // free the space taken up by the random number generator
	gsl_rng_free(r);

}



void theWorld::clearFixedMutations(int whichPop) {
	vector<int> fixedMutations;
	// we take a copy of the first chromosome of the first individual in the population
	// then check it against every chromosome. If an allele is not present in the chromosome being checked, we remove it from the copy we started with
	// the remaining alleles at the end, are the fixed alleles
	fixedMutations = *(popCollection.at(whichPop).thePop.at(0).get_theChromosome(1));
	for(int i=0; i<POPULATION; i++) {
		vector<int> newChromosome1;
		vector<int> newChromosome2;
		newChromosome1 = *(popCollection.at(whichPop).thePop.at(i).get_theChromosome(1));
		newChromosome2 = *(popCollection.at(whichPop).thePop.at(i).get_theChromosome(2));
		int mutationExists;
		// check which mutations in fixedMutations vector are in the newChromosome1
		for(int j=0; j<fixedMutations.size() ;) {
			mutationExists = 0;
			// go through one of the chromosomes and see if the mutations in the first indiv exist in this chromosome
			for(int k=0; k<newChromosome1.size() ; k++) {
				if(fixedMutations.at(j)==newChromosome1.at(k)) {
					mutationExists = 1;
					break;
				}
			}
			// if mutation is not fixed, delete it from the list
			if(mutationExists != 1) {
				fixedMutations.erase(fixedMutations.begin() + j);
			}
			// if mutation exists, keep it and check the next mutation
			else {
				j++;
			}
		}
		// Check the second chromosome
		for(int j=0; j<fixedMutations.size() ;) {
			mutationExists = 0;
			// go through the chromosome and see if the mutations in the first indiv exist in this chromosome
			for(int k=0; k<newChromosome2.size() ; k++) {
				if(fixedMutations.at(j)==newChromosome2.at(k)) {
					mutationExists = 1;
					break;
				}
			}
			// if mutation is not fixed, delete it from the list
			if(mutationExists != 1) {
				fixedMutations.erase(fixedMutations.begin() + j);
			}
			// if mutation exists, keep it and check the next mutation
			else {
				j++;
			}
		}
	}

	// now fixedMutation vector contains all the fixed mutations in the population

	// keep record of number of the fixed mutations
	for(int i=0; i<fixedMutations.size() ; i++) {
		if(fixedMutations.at(i) > 0) {
			numOfDelFixed++;
		}
		else {
			numOfBenFixed++;
		}
	}

	// Now clear these fixed mutations from the population
	for(int i=0; i<POPULATION ; i++) {
		vector<int>* firstChromosome = popCollection.at(whichPop).thePop.at(i).get_theChromosome(1);
		vector<int>* secondChromosome = popCollection.at(whichPop).thePop.at(i).get_theChromosome(2);
		int placeOnChromosome = 0;
		for(int j=0; j<fixedMutations.size() ; j++) {
			for(int k=placeOnChromosome; k<firstChromosome->size() ; k=placeOnChromosome) {
				if(firstChromosome->at(k) == fixedMutations.at(j)) {
					firstChromosome->erase(firstChromosome->begin() + k);
					break;
				}
				else {
					placeOnChromosome++;
				}
			}
		}
		placeOnChromosome = 0;
		for(int j=0; j<fixedMutations.size() ; j++) {
			for(int k=placeOnChromosome; k<secondChromosome->size() ; k=placeOnChromosome) {
				if(secondChromosome->at(k) == fixedMutations.at(j)) {
					secondChromosome->erase(secondChromosome->begin() + k);
					break;
				}
				else {
					placeOnChromosome++;
				}
			}
		}
	}
}



void theWorld::saveThisGeneration(int whichRound) {

    // create a text file to save some info about every individual in the last population
	ofstream generationSaver;

    // create a name for the text file
	stringstream sstr;
	sstr << whichRound;
	string fileName("lastGenerationInfo");
	fileName.append(sstr.str());

	fileName.append(" (");
	if(R_IN_NAME) {
		fileName.append("R=");
		stringstream sstr2;
		sstr2 << popCollection.at(0).thePop.at(0).calcRecombinationRate();
		fileName.append(sstr2.str());
	}
	if(K_IN_NAME) {
		fileName.append(",K=");
		stringstream sstr1;
		sstr1 << K;
		fileName.append(sstr1.str());
	}
	if(P_BEN_IN_NAME) {
		fileName.append(",P_ben=");
		stringstream sstr3;
		sstr3 << P_ben;
		fileName.append(sstr3.str());
	}
	if(S_BEN_IN_NAME) {
		fileName.append(",S_ben=");
		stringstream sstr4;
		sstr4 << S_ben;
		fileName.append(sstr4.str());
	}
	if(ASEX_IN_NAME) {
		fileName.append(",Asp=");
		stringstream sstr5;
		sstr5 << ASEXP;
		fileName.append(sstr5.str());
	}
	fileName.append(")");
	fileName.append(".txt");

	generationSaver.open(fileName.c_str(), ios::out);

    // title of the data being saved for each individual
	generationSaver << "Selfing Rate,Fitness,Del_Hetero,Del_Homo,Ben_Hetero,Ben_Homo" << "\n\n";

    // save the data
	for(int i=0; i<POPULATION ; i++) {
		generationSaver << this->popCollection.at(0).thePop.at(i).calcSelfingRate() << ",";
		generationSaver << (double)this->popCollection.at(0).thePop.at(i).get_fitness() << ",";
		generationSaver << this->popCollection.at(0).thePop.at(i).heteroDeleterious << ",";
		generationSaver << this->popCollection.at(0).thePop.at(i).homoDeleterious << ",";
		generationSaver << this->popCollection.at(0).thePop.at(i).heteroBeneficial << ",";
		generationSaver << this->popCollection.at(0).thePop.at(i).homoBeneficial << "\n";
	}

	generationSaver.close();

}



vector<double>* theWorld::regressionAnalyzer_noBeneficial(int whichPop, int whichRound) {

	vector<double>* resultOfRegression = new vector<double>;

    // This is a vector to save mutations and their related values for every individual in the population
    // x1 = del mutations - mean del mutations
    // x2 = homo del - mean homo del
    // x3 = x1^2
    // x4 =x5 = x1*x2
    // x6 = x2^2
	gsl_matrix* X = gsl_matrix_alloc(POPULATION, 7);
	// a vector to save fitness of each individual relative to mean fitness
	gsl_vector* w = gsl_vector_alloc(POPULATION);
	gsl_vector* B = gsl_vector_calloc(7);
	gsl_vector* e = gsl_vector_alloc(POPULATION);
	gsl_vector* n_homo_del = gsl_vector_alloc(POPULATION);
	gsl_vector* n_het_del = gsl_vector_alloc(POPULATION);
	gsl_vector* n_homo_ben = gsl_vector_alloc(POPULATION);
	gsl_vector* n_het_ben = gsl_vector_alloc(POPULATION);
	gsl_vector* selfingRate = gsl_vector_alloc(POPULATION);
	gsl_vector* X1vals = gsl_vector_alloc(POPULATION);
	gsl_vector* X2vals = gsl_vector_alloc(POPULATION);
	gsl_matrix* covars = gsl_matrix_alloc(7, 7);
	gsl_multifit_linear_workspace* fittingWorkspace = gsl_multifit_linear_alloc(POPULATION, 7);


	double chiSq;

	int meanDel = this->findAvg(1, whichPop);
	int meanHomoDel = this->findAvg(3, whichPop);
	double meanFitness = this->findAvg(5, whichPop);


	for(int i=0; i<POPULATION ; i++) {
		int placeOnMatrix=0;
		// we put a 1 at the first place in order to calculate the intercept
		gsl_matrix_set(X, i, placeOnMatrix++, 1.0);
		double X1 = (this->popCollection.at(whichPop).thePop.at(i).heteroDeleterious + 2*this->popCollection.at(whichPop).thePop.at(i).homoDeleterious) - meanDel;
		gsl_matrix_set(X, i, placeOnMatrix++, X1);
		double X2 = this->popCollection.at(whichPop).thePop.at(i).homoDeleterious - meanHomoDel;
		gsl_matrix_set(X, i, placeOnMatrix++, X2);

		// finding X11, X12, X21, X22
		for(int j=1; j<=2 ; j++) {
			for(int k=j; k<=2 ; k++) {
				double newX1;
				if(j==1)
					newX1 = X1;
				else
					newX1 = X2;

				double newX2;
				if(k==1)
					newX2 = X1;
				else
					newX2 = X2;

				gsl_matrix_set(X, i, placeOnMatrix++, newX1*newX2);
			}
		}

		double relativeW = this->popCollection.at(whichPop).thePop.at(i).get_fitness() / meanFitness;
		gsl_vector_set(w, i, relativeW);

	}

    // find linear fit between relative fitness and the
	gsl_multifit_linear(X, w, B, covars, &chiSq, fittingWorkspace);

	// find the errors in the linear fit found above, and set the values on selfingRate & X1vals & X2vals vectors
	for(int i=0; i<POPULATION ; i++) {
		double error = gsl_vector_get(w, i);
		for(int j=0; j<7 ; j++) {
			error -= gsl_matrix_get(X, i, j)*gsl_vector_get(B,j);
		}
		gsl_vector_set(e, i, error);

		gsl_vector_set(selfingRate, i, this->popCollection.at(whichPop).thePop.at(i).calcSelfingRate());
		gsl_vector_set(X1vals, i , gsl_matrix_get(X, i, 1));
		gsl_vector_set(X2vals, i , gsl_matrix_get(X, i, 2));
		gsl_vector_set(n_homo_del, i , this->popCollection.at(whichPop).thePop.at(i).homoDeleterious);
		gsl_vector_set(n_het_del, i , this->popCollection.at(whichPop).thePop.at(i).heteroDeleterious);
		gsl_vector_set(n_homo_ben, i , this->popCollection.at(whichPop).thePop.at(i).homoBeneficial);
		gsl_vector_set(n_het_ben, i , this->popCollection.at(whichPop).thePop.at(i).heteroBeneficial);
	}

    // find covariances
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, n_homo_del->data, n_homo_del->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, n_homo_ben->data, n_homo_ben->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, n_het_del->data, n_het_del->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, n_het_ben->data, n_het_ben->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, w->data, w->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, X1vals->data, X1vals->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, X2vals->data, X2vals->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, e->data, e->stride, POPULATION));

	resultOfRegression->push_back(gsl_vector_get(B, 1));
	resultOfRegression->push_back(gsl_vector_get(B, 2));
	resultOfRegression->push_back(gsl_vector_get(B, 3));
	resultOfRegression->push_back(gsl_vector_get(B, 4));



    // free up the space taken by the vectors used in the calculations
	gsl_multifit_linear_free(fittingWorkspace);
	gsl_matrix_free(X);
	gsl_matrix_free(covars);
	gsl_vector_free(w);
	gsl_vector_free(B);
	gsl_vector_free(e);
	gsl_vector_free(selfingRate);
	gsl_vector_free(X1vals);
	gsl_vector_free(X2vals);

	return resultOfRegression;
}



vector<double>* theWorld::regressionAnalyzer(int whichPop, int whichRound) {

	vector<double>* resultOfRegression = new vector<double>;

	gsl_matrix* X = gsl_matrix_alloc(POPULATION, 70);
	gsl_vector* w = gsl_vector_alloc(POPULATION);
	gsl_vector* B = gsl_vector_calloc(70);
	gsl_vector* e = gsl_vector_alloc(POPULATION);
	gsl_vector* n_homo_del = gsl_vector_alloc(POPULATION);
	gsl_vector* n_het_del = gsl_vector_alloc(POPULATION);
	gsl_vector* n_homo_ben = gsl_vector_alloc(POPULATION);
	gsl_vector* n_het_ben = gsl_vector_alloc(POPULATION);
	gsl_vector* selfingRate = gsl_vector_alloc(POPULATION);
	gsl_vector* X1vals = gsl_vector_alloc(POPULATION);
	gsl_vector* X2vals = gsl_vector_alloc(POPULATION);
	gsl_matrix* covars = gsl_matrix_alloc(70, 70);
	gsl_multifit_linear_workspace* fittingWorkspace = gsl_multifit_linear_alloc(POPULATION, 70);

	double chiSq;

	int meanDel = this->findAvg(1, whichPop);
	int meanBen = this->findAvg(2, whichPop);
	int meanHomoDel = this->findAvg(3, whichPop);
	int meanHomoBen = this->findAvg(4, whichPop);
	double meanFitness = this->findAvg(5, whichPop);


	for(int i=0; i<POPULATION ; i++) {
		int placeOnMatrix=0;
		// we put a 1 at the first place in order to calculate the intercept
		gsl_matrix_set(X, i, placeOnMatrix++, 1.0);
		double X1 = (this->popCollection.at(whichPop).thePop.at(i).heteroDeleterious + 2*this->popCollection.at(whichPop).thePop.at(i).homoDeleterious) - meanDel;
		gsl_matrix_set(X, i, placeOnMatrix++, X1);
		double X2 = (this->popCollection.at(whichPop).thePop.at(i).heteroBeneficial + 2*this->popCollection.at(whichPop).thePop.at(i).homoBeneficial) - meanBen;
		gsl_matrix_set(X, i, placeOnMatrix++, X2);
		double X3 = this->popCollection.at(whichPop).thePop.at(i).homoDeleterious - meanHomoDel;
		gsl_matrix_set(X, i, placeOnMatrix++, X3);
		double X4 = this->popCollection.at(whichPop).thePop.at(i).homoBeneficial - meanHomoBen;
		gsl_matrix_set(X, i, placeOnMatrix++, X4);

		// finding X11, X12, X13,...
		for(int j=1; j<=4 ; j++) {
			for(int k=j; k<=4 ; k++) {
				double newX1;
				if(j==1)
					newX1 = X1;
				else if(j==2)
					newX1 = X2;
				else if(j==3)
					newX1 = X3;
				else
					newX1 = X4;

				double newX2;
				if(k==1)
					newX2 = X1;
				else if(k==2)
					newX2 = X2;
				else if(k==3)
					newX2 = X3;
				else
					newX2 = X4;

				gsl_matrix_set(X, i, placeOnMatrix++, newX1*newX2);

			}
		}

		// finding X111, X112, X234, X223, ...
		for(int j=1; j<=4 ; j++) {
			for(int k=j; k<=4 ; k++) {
				for(int l=k; l<=4 ; l++) {
					double newX1;
					double newX2;
					double newX3;
					if(j==1)
						newX1 = X1;
					else if(j==2)
						newX1 = X2;
					else if(j==3)
						newX1 = X3;
					else
						newX1 = X4;

					if(k==1)
						newX2 = X1;
					else if(k==2)
						newX2 = X2;
					else if(k==3)
						newX2 = X3;
					else
						newX2 = X4;

					if(l==1)
						newX3 = X1;
					else if(l==2)
						newX3 = X2;
					else if(l==3)
						newX3 = X3;
					else
						newX3 = X4;

					gsl_matrix_set(X, i, placeOnMatrix++, newX1*newX2*newX3);
				}
			}
		}


		// finding X1111, X1233, X1123, X1112, X2234, ...
		for(int j=1; j<=4 ; j++) {
			for(int k=j; k<=4 ; k++) {
				for(int l=k; l<=4 ; l++) {
					for(int m=l; m<=4 ; m++) {
						double newX1;
						double newX2;
						double newX3;
						double newX4;

						if(j==1)
							newX1 = X1;
						else if(j==2)
							newX1 = X2;
						else if(j==3)
							newX1 = X3;
						else
							newX1 = X4;

						if(k==1)
							newX2 = X1;
						else if(k==2)
							newX2 = X2;
						else if(k==3)
							newX2 = X3;
						else
							newX2 = X4;

						if(l==1)
							newX3 = X1;
						else if(l==2)
							newX3 = X2;
						else if(l==3)
							newX3 = X3;
						else
							newX3 = X4;

						if(m==1)
							newX4 = X1;
						else if(m==2)
							newX4 = X2;
						else if(m==3)
							newX4 = X3;
						else
							newX4 = X4;

						gsl_matrix_set(X, i, placeOnMatrix++, newX1*newX2*newX3*newX4);
					}
				}
			}
		}

		double relativeW = this->popCollection.at(whichPop).thePop.at(i).get_fitness() / meanFitness;
		gsl_vector_set(w, i, relativeW);

	}

	gsl_multifit_linear(X, w, B, covars, &chiSq, fittingWorkspace);

	// find the errors in the linear fit found above, and set the values on selfingRate & X1vals & X2vals vectors
	for(int i=0; i<POPULATION ; i++) {
		double error = gsl_vector_get(w, i);
		for(int j=0; j<70 ; j++) {
			error -= gsl_matrix_get(X, i, j)*gsl_vector_get(B,j);
		}
		gsl_vector_set(e, i, error);

		gsl_vector_set(selfingRate, i, this->popCollection.at(whichPop).thePop.at(i).calcSelfingRate());
		gsl_vector_set(X1vals, i , gsl_matrix_get(X, i, 1));
		gsl_vector_set(X2vals, i , gsl_matrix_get(X, i, 2));
		gsl_vector_set(n_homo_del, i , this->popCollection.at(whichPop).thePop.at(i).homoDeleterious);
		gsl_vector_set(n_het_del, i , this->popCollection.at(whichPop).thePop.at(i).heteroDeleterious);
		gsl_vector_set(n_homo_ben, i , this->popCollection.at(whichPop).thePop.at(i).homoBeneficial);
		gsl_vector_set(n_het_ben, i , this->popCollection.at(whichPop).thePop.at(i).heteroBeneficial);
	}

    // find covariances
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, n_homo_del->data, n_homo_del->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, n_homo_ben->data, n_homo_ben->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, n_het_del->data, n_het_del->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, n_het_ben->data, n_het_ben->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, w->data, w->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, X1vals->data, X1vals->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, X2vals->data, X2vals->stride, POPULATION));
	resultOfRegression->push_back(gsl_stats_covariance(selfingRate->data, selfingRate->stride, e->data, e->stride, POPULATION));

	resultOfRegression->push_back(gsl_vector_get(B, 1));
	resultOfRegression->push_back(gsl_vector_get(B, 2));
	resultOfRegression->push_back(gsl_vector_get(B, 3));
	resultOfRegression->push_back(gsl_vector_get(B, 4));

    // free up the memory taken by the vectors used in calculations
	gsl_multifit_linear_free(fittingWorkspace);
	gsl_matrix_free(X);
	gsl_matrix_free(covars);
	gsl_vector_free(w);
	gsl_vector_free(B);
	gsl_vector_free(e);
	gsl_vector_free(selfingRate);
	gsl_vector_free(X1vals);
	gsl_vector_free(X2vals);

	return resultOfRegression;
}


// This is a helper function that calculates frequency of mutations
// If the argument "whichType" is 1, then del mutations, if "whichType" is 2 ben mutations are calculated
void theWorld::findFrequencies(vector<double>& freqArray, int whichType) {

	for(int i=0; i<L ; i++) {
		freqArray[i] = 0;
	}
	// find deleterios frequencies (del mutations have positive values)
	if(whichType==1) {
		for(int i=0; i<POPULATION; i++) {
			for(int j=0; j<this->popCollection.at(0).thePop.at(i).theChromosome1.size() ; j++) {
				if(this->popCollection.at(0).thePop.at(i).theChromosome1.at(j) > 0) {
					freqArray[popCollection.at(0).thePop.at(i).theChromosome1.at(j)-1]++;
				}
			}
			for(int j=0; j<this->popCollection.at(0).thePop.at(i).theChromosome2.size() ; j++) {
				if(this->popCollection.at(0).thePop.at(i).theChromosome2.at(j) > 0) {
					freqArray[popCollection.at(0).thePop.at(i).theChromosome2.at(j)-1]++;
				}
			}
		}
	}

	// find beneficial mutations (ben mutations have negative values)
	if(whichType==2) {
		for(int i=0; i<POPULATION; i++) {
			for(int j=0; j<this->popCollection.at(0).thePop.at(i).theChromosome1.size() ; j++) {
				if(this->popCollection.at(0).thePop.at(i).theChromosome1.at(j) < 0) {
					freqArray[abs(popCollection.at(0).thePop.at(i).theChromosome1.at(j))-1]++;
				}
			}
			for(int j=0; j<this->popCollection.at(0).thePop.at(i).theChromosome2.size() ; j++) {
				if(this->popCollection.at(0).thePop.at(i).theChromosome2.at(j) < 0) {
					freqArray[abs(popCollection.at(0).thePop.at(i).theChromosome2.at(j))-1]++;
				}
			}
		}
	}

    // devide the total number of mutations in each locus by the total number of chromosomes, to find the frequency
	for(int i=0; i<L; i++) {
		freqArray[i] /= (double)(2.0*POPULATION);
	}

}


// find linkage disequilibrium between two loci
double theWorld::findLD_betweenTwo(int firstLocus, int secondLocus, vector<double> delFrequencies) {

	double Dij=0;

	for(int i=0; i<POPULATION ; i++) {

		// First work with Chromosome1
		// Xi(Xj) shows whether the first(second) mutation exists on chromosome1
		int Xi=0;
		int Xj=0;
		// find Xi
		for(int j=0; j<this->popCollection.at(0).thePop.at(i).theChromosome1.size() ; j++) {
			if(this->popCollection.at(0).thePop.at(i).theChromosome1.at(j) == firstLocus) {
				Xi=1;
				break;
			}
			else if(abs(this->popCollection.at(0).thePop.at(i).theChromosome1.at(j)) > firstLocus) {
				break;
			}
		}
		// find Xj
		for(int j=0; j<this->popCollection.at(0).thePop.at(i).theChromosome1.size() ; j++) {
			if(this->popCollection.at(0).thePop.at(i).theChromosome1.at(j) == secondLocus) {
				Xj=1;
				break;
			}
			else if(abs(this->popCollection.at(0).thePop.at(i).theChromosome1.at(j)) > secondLocus) {
				break;
			}
		}

		Dij += (double)((Xi - delFrequencies[firstLocus-1])*(Xj - delFrequencies[secondLocus-1]));

		// Now work with Chromosome2
		// Xi(Xj) shows whether the first(second) mutation exists on chromosome2
		Xi=0;
		Xj=0;
		// find Xi
		for(int j=0; j<this->popCollection.at(0).thePop.at(i).theChromosome2.size() ; j++) {
			if(this->popCollection.at(0).thePop.at(i).theChromosome2.at(j) == firstLocus) {
				Xi=1;
				break;
			}
			else if(abs(this->popCollection.at(0).thePop.at(i).theChromosome2.at(j)) > firstLocus) {
				break;
			}
		}
		// find Xj
		for(int j=0; j<this->popCollection.at(0).thePop.at(i).theChromosome2.size() ; j++) {
			if(this->popCollection.at(0).thePop.at(i).theChromosome2.at(j) == secondLocus) {
				Xj=1;
				break;
			}
			else if(abs(this->popCollection.at(0).thePop.at(i).theChromosome2.at(j)) > secondLocus) {
				break;
			}
		}

		Dij += (double)((Xi - delFrequencies[firstLocus-1])*(Xj - delFrequencies[secondLocus-1]));

	}

	return Dij;

}



// find Linkage Disequilibrium in the population by choosing 5000 random pair of loci and finding LD between them
void theWorld::findLD(long double& avgFreqOfDelAlleles, long double& avgPolymorphOfDelAlleles, double& avgDij, long double& avgD_prime) {
	vector<double> delFrequencies(L, 0.0);
	vector<double> benFrequencies(L, 0.0);

	this->findFrequencies(delFrequencies, 1);
	this->findFrequencies(benFrequencies, 2);

	avgD_prime=0;
	double helperToFindD_Prime=0;
	long double total_Dij=0;

	// Initialie random number generator
	const gsl_rng_type * T;
	gsl_rng * r;
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_default_seed += time(NULL);
	gsl_rng_set (r, gsl_rng_default_seed);

	for(int n=0; n<5000 ; n++) {
	    // find 5000 random pairs of loci
		int i = floor(gsl_ran_flat(r, 0.5, L + 0.499999) + 0.5);
		gsl_rng_default_seed += time(NULL)*i;
		int j = floor(gsl_ran_flat(r, 0.5, L + 0.499999) + 0.5);
		gsl_rng_default_seed += time(NULL)*j;
        // check that the two loci are not the same
		if(i == j) {
			n--;
			continue;
		}

		double Dij = this->findLD_betweenTwo(i,j, delFrequencies);
		double maxDij;
		if(Dij<0) {
			maxDij = delFrequencies[i-1]*delFrequencies[j-1] < (1-delFrequencies[i-1])*(1-delFrequencies[j-1]) ? delFrequencies[i-1]*delFrequencies[j-1] : (1-delFrequencies[i-1])*(1-delFrequencies[j-1]);
		}
		else {
			maxDij = delFrequencies[i-1]*(1-delFrequencies[j-1]) < delFrequencies[j-1]*(1-delFrequencies[i-1]) ? delFrequencies[i-1]*(1-delFrequencies[j-1]) : delFrequencies[j-1]*(1-delFrequencies[i-1]);
		}
		if(maxDij != 0) {
			avgD_prime += abs(Dij/maxDij);
			helperToFindD_Prime += 1.0;
		}
		total_Dij += Dij;
	}

    // free up the space taken by the random number generator
	gsl_rng_free(r);

	long double avgFreqOfDel=0, avgPolyOfDel=0;
	for(int i=0; i<L ; i++) {
		avgFreqOfDel += delFrequencies[i];
		avgPolyOfDel += delFrequencies[i] * (1-delFrequencies[i]);
	}
	avgFreqOfDel /= (long double)(L);
	avgPolyOfDel /= (long double)(L);

	avgFreqOfDelAlleles = avgFreqOfDel;
	avgPolymorphOfDelAlleles = avgPolyOfDel;
	avgDij = total_Dij/5000.0;
	avgD_prime /= helperToFindD_Prime;
}


// This function takes the mean slefing rates in the last 3000 generations (saved every 20 generations) and checks if the population is at equilibrium
// Checking for equilibrium: 1. divide the last 3000 generations into 6 sections of 500 generations
// 2. find the mean value in each of the 6 sections
// 3. If the difference between the max and min values is less than 0.02, equilibrium is reached
int theWorld::isAtEquil(vector<double> selfingRates) {

	double maxAvg=0;
	double minAvg=1;

	for(int i=0; i<6 ; i++) {
		double avgSelfing=0;
		if((int)(selfingRates.size()) - 25*i >= 25) {
			for(int j=0; j<25; j++) {
				avgSelfing += selfingRates.at((int) i*25 + j);
			}
			avgSelfing /= 25.0;
			if(maxAvg < avgSelfing) {
				maxAvg = avgSelfing;
			}
			if(minAvg > avgSelfing) {
				minAvg = avgSelfing;
			}
		}
	}
	if(maxAvg - minAvg < 0.02)
		return 1;

	return 0;
}



// two functions give us some of the info stored in a "theWorld" object
vector<studyPop>* theWorld::get_popCollection() {
	return &popCollection;
}

vector<int>* theWorld::get_generationNo() {
	return &generationNo;
}
