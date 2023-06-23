#include "zHardCoding.h"

int zIsFivePrimeOverhang(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "RGenomic1") == 0 || 
		strcmp(zStrIdx2Char(state), "RCDna1")    == 0 );
}

int zIsThreePrimeOverhang(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "RGenomic2") == 0 || 
		strcmp(zStrIdx2Char(state), "RCDna2")    == 0 );
}

int zIsGenomicOnly(zStrIdx state) {
	return (zIsIntron(state) || 
	        strcmp(zStrIdx2Char(state), "RGenomic1") == 0 || 
	        strcmp(zStrIdx2Char(state), "RGenomic2") == 0 || 
	        strcmp(zStrIdx2Char(state), "Genomic")   == 0 );
}

int zIsCDnaOnly(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "RCDna1")    == 0 || 
	        strcmp(zStrIdx2Char(state), "RCDna2")    == 0 || 
	        strcmp(zStrIdx2Char(state), "CDna")      == 0 );
}

int zIsGenomicInsertion(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "Genomic")   == 0 ); 
}

int zIsCDnaInsertion(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "CDna")      == 0 ); 
}

int zIsU2Entry(zStrIdx state, strand_t strand) {
	return ((strcmp(zStrIdx2Char(state), "DonorU2")  == 0 && strand == '+') ||
	        (strcmp(zStrIdx2Char(state), "AccU2")    == 0 && strand == '-'));
}

int zIsU2Exit(zStrIdx state, strand_t strand) {
	return ((strcmp(zStrIdx2Char(state), "DonorU2")  == 0 && strand == '-') ||
	        (strcmp(zStrIdx2Char(state), "AccU2")    == 0 && strand == '+'));
}

int zIsU12Entry(zStrIdx state, strand_t strand) {
	return ((strcmp(zStrIdx2Char(state), "DonorU12")  == 0 && strand == '+') ||
	        (strcmp(zStrIdx2Char(state), "AccU12")    == 0 && strand == '-'));
}

int zIsU12Exit(zStrIdx state, strand_t strand) {
	return ((strcmp(zStrIdx2Char(state), "DonorU12")  == 0 && strand == '-') ||
	        (strcmp(zStrIdx2Char(state), "AccU12")    == 0 && strand == '+'));
}

int zIsIntronEntry(zStrIdx state, strand_t strand) {
	return (zIsU2Entry(state, strand) ||
	        strcmp(zStrIdx2Char(state), "Entry")  == 0 || 
	        strcmp(zStrIdx2Char(state), "EntryU12")  == 0 || 
	        strcmp(zStrIdx2Char(state), "EntryU2")  == 0 || 
	        zIsU12Entry(state, strand));
}

int zIsIntronExit(zStrIdx state, strand_t strand) {
	return (zIsU2Exit(state, strand) ||
	        strcmp(zStrIdx2Char(state), "Exit")  == 0 || 
	        strcmp(zStrIdx2Char(state), "ExitU2")  == 0 || 
	        strcmp(zStrIdx2Char(state), "ExitU12")  == 0 || 
	        zIsU12Exit(state, strand));
}

int zIsBranchPoint(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "BranchU2")  == 0 ||
	        strcmp(zStrIdx2Char(state), "BranchU12") == 0 );
}

int zIsInsideU2Intron(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "IntronU2")  == 0 ||
		strcmp(zStrIdx2Char(state), "IntronC")  == 0 || 
	        strcmp(zStrIdx2Char(state), "BranchU2")  == 0 || 
	        strcmp(zStrIdx2Char(state), "BrAccU2")   == 0 ); 
}

int zIsInsideU12Intron(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "IntronU12") == 0 ||
	        strcmp(zStrIdx2Char(state), "BranchU12") == 0 || 
	        strcmp(zStrIdx2Char(state), "BrAccU12")  == 0 ); 
}

int zIsU2(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "IntronU2") == 0 ||
	        strcmp(zStrIdx2Char(state), "DonorU2" ) == 0 || 
	        strcmp(zStrIdx2Char(state), "AccU2"  ) == 0 || 
	        strcmp(zStrIdx2Char(state), "BranchU2") == 0 || 
	        strcmp(zStrIdx2Char(state), "BrAccU2")  == 0 ); 
}

int zIsU12(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "IntronU12") == 0 ||
	        strcmp(zStrIdx2Char(state), "DonorU12")  == 0 || 
	        strcmp(zStrIdx2Char(state), "AccU12")   == 0 || 
	        strcmp(zStrIdx2Char(state), "BranchU12") == 0 || 
	        strcmp(zStrIdx2Char(state), "BrAccU12")  == 0 ); 
}

int zIsInsideIntron(zStrIdx state) {
	return (zIsInsideU2Intron(state) ||
		zIsInsideU12Intron(state));
}

int zIsMatch(zStrIdx state) {
	return (strcmp(zStrIdx2Char(state), "Match") == 0); 
}

int zIsIntron(zStrIdx state) {
	return (zIsIntronEntry(state, '+')  || 
		zIsInsideIntron(state) || 
		zIsIntronExit(state, '+'));
}

int zIsExon(zStrIdx state) {
	return (zIsMatch(state) || zIsGenomicInsertion(state) || zIsCDnaInsertion(state));
}

int zIsOverhang(zStrIdx state) {
	return (zIsFivePrimeOverhang(state) || 
		zIsThreePrimeOverhang(state)); 
}

int zIsGenomicOverhang(zStrIdx state) {
	return zIsGenomicOnly(state) && zIsOverhang(state);
}

int zIsCDnaOverhang(zStrIdx state) {
	return zIsCDnaOnly(state) && zIsOverhang(state);
}

int zGetStateByName(zHMM* hmm, const char* name) {
	int i;
	for (i = 0; i < hmm->states; i++) {
		if (strcmp(zStrIdx2Char(hmm->state[i].name), name) == 0) {
			return i;
		}
	}
	return -1;
}

int zGetU2Donor(zHMM* hmm) {
	return zGetStateByName(hmm, "DonorU2");
}

int zGetU2Acceptor(zHMM* hmm) {
	return zGetStateByName(hmm, "AccU2");
}

int zGetU12Donor(zHMM* hmm) {
	return zGetStateByName(hmm, "DonorU12");
}

int zGetU12Acceptor(zHMM* hmm) {
	return zGetStateByName(hmm, "AccU12");
}

int zGetU2Entry(zHMM* hmm) {
	return zGetStateByName(hmm, "EntryU2");
}

int zGetU2Exit(zHMM* hmm) {
	return zGetStateByName(hmm, "ExitU2");
}

int zGetU12Entry(zHMM* hmm) {
	return zGetStateByName(hmm, "EntryU12");
}

int zGetU12Exit(zHMM* hmm) {
	return zGetStateByName(hmm, "ExitU12");
}

int zGetFivePrimeGenomic(zHMM* hmm) {
	return zGetStateByName(hmm, "RGenomic1");
}

int zGetFivePrimeCDna(zHMM* hmm) {
	return zGetStateByName(hmm, "RCDna1");
}

int zGetThreePrimeGenomic(zHMM* hmm) {
	return zGetStateByName(hmm, "RGenomic2");
}

int zGetThreePrimeCDna(zHMM* hmm) {
	return zGetStateByName(hmm, "RCDna2");
}

int zGetMatch(zHMM* hmm) {
	int i;
	for (i = 0; i < hmm->states; i++) {
		if (zIsMatch(hmm->state[i].name)) {
			return i;
		}
	}
	return -1;
}

int zGetU12BranchPoint(zHMM* hmm) {
	return zGetStateByName(hmm, "BranchU12");
}

/* Set the models according to the strand. Should be used only in pairHMM */
int zSetHMMStrand(zHMM* hmm, strand_t strand) {
	int     entry_state,  exit_state, i;
	zModel *entry_model, *exit_model;

	if (hmm->mode != GPAIRHMM) {
		zDie("Cannot set HMM strand for non-GPAIRHMM mode");
	}

	if (strand == '.') {
		zDie("Cannot set HMM strand to '.'");
	}

	/* strand already in place */
	entry_state = zGetMatch(hmm);
	if (hmm->state[entry_state].strand == strand) {
		return 1;
	}

	/* Anti all models */

	/* DO NOT SHARE MODELS BETWEEN STATES */
	/* Remember that by design states can share models. If the same model
         * is used twice, you will anti it twice resulting in no change.
         * This is ok with models that are the same as NULL model when using
         * NULL model. However, in my opinion, the best solution is to use
         * different models for different states, duplicating identical models */

	for (i = 0; i < hmm->states; i++) {
		zModel* model = hmm->mmap[hmm->state[i].model];
		if (zAntiModel(model) == -1) return -1;
		hmm->state[i].strand = strand;
	}

	if (zGetU12BranchPoint(hmm) == -1) {

		/* No Branch point models, so it is easy to swap them */

		int buffer;

		/* Swap Entry/Exit for U2 introns */

		/* Swap the models */

		entry_state = zGetU2Entry(hmm);
		exit_state  = zGetU2Exit(hmm);
		entry_model = hmm->mmap[hmm->state[entry_state].name];
		exit_model  = hmm->mmap[hmm->state[exit_state].name];
		hmm->mmap[hmm->state[entry_state].name] = exit_model;
		hmm->mmap[hmm->state[exit_state].name]  = entry_model;

		/* Swap the increments */

		buffer      = hmm->increments[entry_state][0];
		hmm->increments[entry_state][0] = hmm->increments[exit_state][0];
		hmm->increments[exit_state][0]  = buffer;

		buffer      = hmm->increments[entry_state][1];
		hmm->increments[entry_state][1] = hmm->increments[exit_state][1];
		hmm->increments[exit_state][1]  = buffer;

		/* Swap Entry/Exit for U12 introns */

		entry_state = zGetU12Entry(hmm);
		exit_state  = zGetU12Exit(hmm);
		entry_model = hmm->mmap[hmm->state[entry_state].name];
		exit_model  = hmm->mmap[hmm->state[exit_state].name];
		hmm->mmap[hmm->state[entry_state].name] = exit_model;
		hmm->mmap[hmm->state[exit_state].name]  = entry_model;

		/* Swap the increments */

		buffer      = hmm->increments[entry_state][0];
		hmm->increments[entry_state][0] = hmm->increments[exit_state][0];
		hmm->increments[exit_state][0]  = buffer;

		buffer      = hmm->increments[entry_state][1];
		hmm->increments[entry_state][1] = hmm->increments[exit_state][1];
		hmm->increments[exit_state][1]  = buffer;

	} else {
		if (strand == '+') {
			/*
			All tr to entry will be from exit
			All tr from exit will be to entry
			Intron will be reversed
			*/
			int j;
			int state;
			state = zGetU2Acceptor(hmm);
			for (i = 0; i < hmm->states; i++) {
				if (hmm->tmap[i][state][0] != MIN_SCORE) {
					if (hmm->tmap[state][i][0] != MIN_SCORE) {
						zDie("Something is wrong with transitions");
					} 
					hmm->tmap[state][i][0] = hmm->tmap[i][state][0];
					hmm->tmap[i][state][0] = MIN_SCORE;
				}
			}
			state = zGetU12Acceptor(hmm);
			for (i = 0; i < hmm->states; i++) {
				if (hmm->tmap[i][state][0] != MIN_SCORE) {
					if (hmm->tmap[state][i][0] != MIN_SCORE) {
						zDie("Something is wrong with transitions");
					} 
					hmm->tmap[state][i][0] = hmm->tmap[i][state][0];
					hmm->tmap[i][state][0] = MIN_SCORE;
				}
			}
			state = zGetU2Donor(hmm);
			for (i = 0; i < hmm->states; i++) {
				if (hmm->tmap[state][i][0] != MIN_SCORE) {
					if (hmm->tmap[i][state][0] != MIN_SCORE) {
						zDie("Something is wrong with transitions");
					} 
					hmm->tmap[i][state][0] = hmm->tmap[state][i][0];
					hmm->tmap[state][i][0] = MIN_SCORE;
				}
			}
			state = zGetU12Donor(hmm);
			for (i = 0; i < hmm->states; i++) {
				if (hmm->tmap[state][i][0] != MIN_SCORE) {
					if (hmm->tmap[i][state][0] != MIN_SCORE) {
						zDie("Something is wrong with transitions");
					} 
					hmm->tmap[i][state][0] = hmm->tmap[state][i][0];
					hmm->tmap[state][i][0] = MIN_SCORE;
				}
			}

			/* Now we are done with the x -> entry and exit -> x transitions. We need to fix the
			 * entry -> * -> exit by reversing the direction. Since the tmaps are half-way fixed,
			 * we shouldnt use that to reverse introns. Use jmaps instead. */

			exit_state = zGetU2Donor(hmm);
			if (hmm->jmap[exit_state]->size != 1) zDie("Something is wrong with transitions");
			while (exit_state != zGetU2Acceptor(hmm)) {
				int pre;
				if (hmm->jmap[exit_state]->size > 2) zDie("Something is wrong with transitions for %d", exit_state);
				pre = hmm->jmap[exit_state]->elem[0];
				if (pre == exit_state) {
					pre = hmm->jmap[exit_state]->elem[1];
				}
				hmm->tmap[exit_state][pre][0] = hmm->tmap[pre][exit_state][0];
				hmm->tmap[pre][exit_state][0] = MIN_SCORE;
				exit_state = pre;
			}

			exit_state = zGetU12Donor(hmm);
			if (hmm->jmap[exit_state]->size != 1) zDie("Something is wrong with transitions");
			while (exit_state != zGetU12Acceptor(hmm)) {
				int pre;
				if (hmm->jmap[exit_state]->size > 2) zDie("Something is wrong with transitions for %d", exit_state);
				pre = hmm->jmap[exit_state]->elem[0];
				if (pre == exit_state) {
					pre = hmm->jmap[exit_state]->elem[1];
				}
				hmm->tmap[exit_state][pre][0] = hmm->tmap[pre][exit_state][0];
				hmm->tmap[pre][exit_state][0] = MIN_SCORE;
				exit_state = pre;
			}

			for (i = 0; i < hmm->states; i++) {
				zFreeIVec(hmm->jmap[i]);
				zInitIVec(hmm->jmap[i], 1);
				for (j = 0; j < hmm->states; j++) {
					if (hmm->tmap[j][i][0] != MIN_SCORE) zPushIVec(hmm->jmap[i], j); /* hard coding tmap[][][0] is OK since we
													    are only checking for presence of 
													    transitions */
													 /* check this assumption : 05/21/03 */
				}
			}
		} else if (strand == '-') {
			/*
			All tr to entry will be from exit
			All tr from exit will be to entry
			Intron will be reversed
			*/
			int j;
			int state;
			state = zGetU2Donor(hmm);
			for (i = 0; i < hmm->states; i++) {
				if (hmm->tmap[i][state][0] != MIN_SCORE) {
					if (hmm->tmap[state][i][0] != MIN_SCORE) {
						zDie("Something is wrong with transitions");
					} 
					hmm->tmap[state][i][0] = hmm->tmap[i][state][0];
					hmm->tmap[i][state][0] = MIN_SCORE;
				}
			}
			state = zGetU12Donor(hmm);
			for (i = 0; i < hmm->states; i++) {
				if (hmm->tmap[i][state][0] != MIN_SCORE) {
					if (hmm->tmap[state][i][0] != MIN_SCORE) {
						zDie("Something is wrong with transitions");
					} 
					hmm->tmap[state][i][0] = hmm->tmap[i][state][0];
					hmm->tmap[i][state][0] = MIN_SCORE;
				}
			}
			state = zGetU2Acceptor(hmm);
			for (i = 0; i < hmm->states; i++) {
				if (hmm->tmap[state][i][0] != MIN_SCORE) {
					if (hmm->tmap[i][state][0] != MIN_SCORE) {
						zDie("Something is wrong with transitions");
					} 
					hmm->tmap[i][state][0] = hmm->tmap[state][i][0];
					hmm->tmap[state][i][0] = MIN_SCORE;
				}
			}
			state = zGetU12Acceptor(hmm);
			for (i = 0; i < hmm->states; i++) {
				if (hmm->tmap[state][i][0] != MIN_SCORE) {
					if (hmm->tmap[i][state][0] != MIN_SCORE) {
						zDie("Something is wrong with transitions");
					} 
					hmm->tmap[i][state][0] = hmm->tmap[state][i][0];
					hmm->tmap[state][i][0] = MIN_SCORE;
				}
			}

			/* Now we are done with the x -> entry and exit -> x transitions. We need to fix the
			 * entry -> * -> exit by reversing the direction. Since the tmaps are half-way fixed,
			 * we shouldnt use that to reverse introns. Use jmaps instead. */

			exit_state = zGetU2Acceptor(hmm);
			if (hmm->jmap[exit_state]->size != 1) zDie("Something is wrong with transitions");
			while (exit_state != zGetU2Donor(hmm)) {
				int pre;
				if (hmm->jmap[exit_state]->size > 2) zDie("Something is wrong with transitions for %d", exit_state);
				pre = hmm->jmap[exit_state]->elem[0];
				if (pre == exit_state) {
					pre = hmm->jmap[exit_state]->elem[1];
				}
				hmm->tmap[exit_state][pre][0] = hmm->tmap[pre][exit_state][0];
				hmm->tmap[pre][exit_state][0] = MIN_SCORE;
				exit_state = pre;
			}

			exit_state = zGetU12Acceptor(hmm);
			if (hmm->jmap[exit_state]->size != 1) zDie("Something is wrong with transitions");
			while (exit_state != zGetU12Donor(hmm)) {
				int pre;
				if (hmm->jmap[exit_state]->size > 2) zDie("Something is wrong with transitions for %d", exit_state);
				pre = hmm->jmap[exit_state]->elem[0];
				if (pre == exit_state) {
					pre = hmm->jmap[exit_state]->elem[1];
				}
				hmm->tmap[exit_state][pre][0] = hmm->tmap[pre][exit_state][0];
				hmm->tmap[pre][exit_state][0] = MIN_SCORE;
				exit_state = pre;
			}

	/*
			for (i = 0; i < hmm->states; i++) {
				for (j = 0; j < hmm->states; j++) {
					if (hmm->tmap[i][j][0] == MIN_SCORE) printf("\t.");
					else printf("\t%5.0f", hmm->tmap[i][j][0]);
				}
				printf("\n");
			}
	*/
			for (i = 0; i < hmm->states; i++) {
				zFreeIVec(hmm->jmap[i]);
				zInitIVec(hmm->jmap[i], 1);
				for (j = 0; j < hmm->states; j++) {
					if (hmm->tmap[j][i][0] != MIN_SCORE) zPushIVec(hmm->jmap[i], j); /* hard coding tmap[][][0] is OK since we
													    are only checking for presence of 
													    transitions */
													 /* check this assumption : 05/21/03 */
				}
			}
		}
	}

	/* Swap 5'/3' unaligned cDNA (sorry about the variable names) */

	entry_state = zGetFivePrimeCDna(hmm);
	exit_state  = zGetThreePrimeCDna(hmm);
	entry_model = hmm->mmap[hmm->state[entry_state].name];
	exit_model  = hmm->mmap[hmm->state[exit_state].name];
	hmm->mmap[hmm->state[entry_state].name] = exit_model;
	hmm->mmap[hmm->state[exit_state].name]  = entry_model;

	return 1;
}
