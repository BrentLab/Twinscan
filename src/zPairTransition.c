/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zTransition.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_TRANSITION_C
#define ZOE_TRANSITION_C

#include "zTransition.h"
#include "zPairTrellis.h"
#include "zFeatureFactory.h" /* for PSC definition */

/******************************************************************************\
 Transition Functions
 
\******************************************************************************/

score_t zGetScannerScore (zPairTrellis *trellis, zScanner *scanner, int state, coor_t genomic, coor_t cdna) {
	score_t total_score = 0;
/* The order here is GENOMIC, DNA, PAIR since the number of corresponding states are descending in that order */
/* Calling zGetUScore since it has been precomputed */
 	if (scanner->model->seq_type == GENOMIC) {
	 	total_score = zGetUScore(scanner, genomic);
	} else if (scanner->model->seq_type == DNA) {
	 	total_score = zGetUScore(scanner, cdna);
	} else if (scanner->model->seq_type == PAIR) {
		if (scanner->pairscore != NULL) {
			total_score = scanner->pairscore(scanner, trellis->genomic, trellis->cdna, genomic, cdna);
		} else {
			zDie("PairHMM scoring function absent in scanner: %s", zStrIdx2Char(trellis->hmm->state[state].model));
		}
	}
	return total_score;
}

/*****************************************************************************************
 * The GHMM version of this function (zScoreInternalState) takes in a bool first_base    *
 * argument to add the exit probability for this state if it is the first base of this   *
 * state. There are two cases: this is the first base in the sequence itself as well     *
 * as the state, this is the first base in the state after transitioning from a previous *
 * state. The exit probabilities for GPAIRHMM are added to the initial probabilities of  *
 * GEOMETRIC/Internal/GInternal states and transition to GEOMETRIC/Internal/GInternal    *
 * states by calling zGetFixedInitProb instead of zGetInitProb, and                      *
 * zFixInternalTransitions(hmm) once before Viterbi, resp. This saves an additional      *
 * comparison and scoring of the exit probability inside Viterbi which turns out to be   *
 * very expensive in GPAIRHMM. zFixInternalTransitions is called by zReadHMM if it is    *
 * a GPAIRHMM.                                                                           *
 *                                                                                       *
 * See also: zHMM.c:zGetFixedInitProb(), zHMM.c:zFixInternalTransitions() and comments   *
 * therein.                                                                              *
 *****************************************************************************************/
/* this function scores a single base of an internal state */

void zInternalPairTransHelper (zPairTrellis *trellis, int from_state, int to_state, coor_t genomic, coor_t cdna, score_t trans) {
	score_t    total_score, cell_score;
	zPairTrellisCell *cell, *prev_cell;

	
	/* check phase */
	if (EXTERNAL == trellis->hmm->state[from_state].type) {
		zPhase_t int_phase, ext_phase;
		int_phase = trellis->hmm->state[to_state].phase;
		ext_phase = trellis->hmm->state[from_state].phase;
		if (ext_phase != int_phase) return;
	}
	
	total_score = zGetScannerScore(trellis,trellis->scanner[trellis->hmm->state[to_state].model],to_state,genomic,cdna);
	total_score += trans;

	if (total_score == MIN_SCORE) return;

#ifdef FORWARD
	if (trellis->forward != NULL) {
		trellis->forward[genomic][cdna][to_state] = 
			zFloatwiseScoreAdd(total_score + trellis->forward[genomic-zGetGenomicIncrement(trellis->hmm, to_state)][cdna-zGetCDnaIncrement(trellis->hmm, to_state)][from_state], 
								trellis->forward[genomic][cdna][to_state]);
	}
#endif
	
	cell      = zGetCurrentCell(trellis, genomic, cdna, to_state);
	prev_cell = zGetPreviousCell(trellis, genomic, cdna, to_state, from_state);

	cell_score = total_score + prev_cell->score;
	if (cell_score > cell->score) {
		cell->score = cell_score;
		if (from_state == to_state) {
			cell->length = prev_cell->length + 1;
			cell->trace = prev_cell->trace;
		} else {
			cell->trace = from_state;
			cell->length = 1;
		}
	}
	return; 
}

static void zInternalPairTrans (zPairTrellis *trellis, int from_state, int to_state, coor_t genomic, coor_t cdna) {
	zInternalPairTransHelper(trellis, from_state, to_state, genomic, cdna, zGetTransitionScore(trellis->hmm,from_state,to_state,trellis->tiso_group));
}

static void zInternalPairTransBack (zPairTrellis *trellis, int from_state, int to_state, coor_t genomic, coor_t cdna) {
	score_t    total_score;
	int        gincrement = zGetGenomicIncrement(trellis->hmm, to_state);
	int        cincrement = zGetCDnaIncrement(trellis->hmm, to_state);

	if (EXTERNAL == trellis->hmm->state[from_state].type) {
		zPhase_t int_phase, ext_phase;
		int_phase = trellis->hmm->state[to_state].phase;
		ext_phase = trellis->hmm->state[from_state].phase;
		if (ext_phase != int_phase) return;
	}

	total_score = zGetScannerScore(trellis,trellis->scanner[trellis->hmm->state[to_state].model],to_state,genomic,cdna);
	total_score += zGetTransitionScore(trellis->hmm, from_state, to_state, trellis->tiso_group);
	if (total_score == MIN_SCORE) return;
	
	trellis->backward[genomic][cdna][from_state] = 
		  zFloatwiseScoreAdd(total_score + trellis->backward[genomic+gincrement][cdna+cincrement][to_state],
							trellis->backward[genomic][cdna][from_state]);
	return;
}

/* Assumption: genomic sequence moves in this state: increments for genomic > 0 */
void zExplicitPairTrans (zPairTrellis *trellis, int from_state, int state, coor_t genomic, coor_t cdna) {
	zStrIdx    name;
	score_t    tscore, scan_score, total_score, best_score;
	coor_t     best_length = 0, length;
	zPairTrellisCell *cell; 
	zScanner *scanner;
	int gincrement, cincrement;

	zPhase_t int_phase, ext_phase;
	zDurationGroup *group;
	zDistribution *d;
	coor_t gmin, gmax, cmin, cmax;
	coor_t steps;

	if (EXTERNAL == trellis->hmm->state[from_state].type) 
	{ /*  check phase */
			int_phase = trellis->hmm->state[state].phase;
			ext_phase = trellis->hmm->state[from_state].phase;
			if (ext_phase != int_phase) return;
	} 

	gincrement = zGetGenomicIncrement(trellis->hmm, state);
	cincrement = zGetCDnaIncrement(trellis->hmm, state);
	group     = trellis->hmm->dmap[trellis->hmm->state[state].duration];
	d = &group->duration[0].distribution[0];

	/* Set min and max start coordinates for the range */
	steps = d->end - d->start + 1;
	if ((int)(genomic - steps*gincrement) < (int)trellis->blocks->gb_start - 1) {
		steps = (genomic - trellis->blocks->gb_start + 1)/gincrement;
	}
	if ((int)(cdna - steps*cincrement) < PADDING - 1) {
		steps = (cdna - PADDING + 1)/cincrement;
	}
	
	name      = trellis->hmm->state[state].model;
	scanner   = trellis->scanner[name];
                        
	scan_score = 0.;
	best_score  = MIN_SCORE;
	tscore  = zGetTransitionScore(trellis->hmm, from_state, state, trellis->tiso_group);

	for (length = 1, gmax = genomic, cmax = cdna, gmin = gmax - gincrement, cmin = cmax - cincrement;
	     length <= steps;
	     length++, gmax = gmin, cmax = cmin, gmin = gmax - gincrement, cmin = cmax - cincrement) {

		coor_t i, j;

		total_score = tscore
				+ zScoreDistribution(d, length)
				+ zGetCurrentCell(trellis, gmin, cmin, from_state)->score;
			
		for (i = gmin+1, j = cmin+1; i <= gmax; i+=gincrement, j+=cincrement) {
			scan_score += zGetScannerScore(trellis, scanner, state, i, j);
		}
		total_score += scan_score;

		if (total_score >= best_score) {
			best_score = total_score;
			best_length = length;
		}
	  
	}
               
	cell = zGetCurrentCell(trellis, genomic, cdna, state);
	if (best_score > cell->score) {
		cell->score = best_score;
		cell->length = best_length;
		cell->trace = from_state;
	}
}

#define STATE_TYPES 4
static zPairTransFunc zPairTransLookup[STATE_TYPES] = {
	zInternalPairTrans, /* INTERNAL */
	zInternalPairTrans, /* GINTERNAL */
	zInternalPairTrans, /* EXTERNAL */
	zExplicitPairTrans  /* EXPLICIT */
};

zPairTransFunc zGetPairTransFunc(int type) {
 	if (type >= STATE_TYPES || type < 0) {
		zWarn("Attempt to get invalid state type");
		return NULL;
	}

	return zPairTransLookup[type]; 
}

static zPairTransFunc zPairBackTransLookup[STATE_TYPES] = {
	zInternalPairTransBack, /* INTERNAL */
	zInternalPairTransBack, /* INTERNAL */
	zInternalPairTransBack, /* INTERNAL */
	zExplicitPairTrans
};

zPairTransFunc zGetPairBackTransFunc(int type) {
	if (type >= STATE_TYPES || type < 0) {
		zWarn("Attempt to get invalid state type");
		return NULL;
	}
	return zPairBackTransLookup[type];
}

#endif
