/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zTransition.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_TRANSITION_H
#define ZOE_TRANSITION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zMath.h"
#include "zSfeature.h"
#include "zTools.h"
#include "zDNA.h"
#include "zEstseq.h"

/******************************************************************************\
 zTransition

zTransition is used for hmm state transitions.

\******************************************************************************/

struct zTransition {
        zStrIdx  from;  /* state id */
        zStrIdx  to;    /* state id */

        /* the following arrays are indexed by the isochore groups. */
        /* prob[0] is isochore group 0 and so on. */

        float    *prob;  /* probability of transition */
        score_t  *score; /* score of transition */
};
typedef struct zTransition zTransition;

int  zReadTransition (FILE*, zTransition*, int);
void zWriteTransition (FILE*, const zTransition*, int);

/******************************************************************************\
 Transition Functions
 
Zoe was designed to use internal transitions and external
transitions.  That functionality can be redefined here, using
zRegisterTransFunc.

A transition function has a signature looks like:

	void zExternalTrans(zTrellis* trellis, int from_state, int to_state
					coor_t pos) {

For a transition from state i to state j, zTrellis calls the
transition function of state i.  The responsibility of the
transition function is to update the trellis cell specified by
its parameters, if and only if the transition is better than the
one currently in the cell.  If the traceback is to an HMM state,
then the traceback entry should just be the index of the state. If
the traceback is to a feature, then the traceback entry sholud be
-1 * (index of the feature in trellis->keep).  Note:  any time a
feature is added to the keep, a corresponding entry should be added
to the trellis->jump vector for the state prior to the feature.

For use in creating other scoring methods, two helper functions
are available:

	zCompatibleJumpFromExon() - checks for phase compatibility
	zCompatibleJumpToExon() - nature doesn't seem to like stop
		codons in the middle of genes.

Note: With the addition of the forward-backward algorithm, this method
is straining.  The logic of transitions is scattered and repeated in 
several different functions.  Any ideas to consolidate this would be 
greatly appreciated.

\******************************************************************************/

struct zTrellis;
typedef void (*zTransFunc)(struct zTrellis* trel, int from, int to, coor_t j);

zTransFunc zGetTransFunc(int type); 
zTransFunc zGetBackTransFunc(int type);

int zCompatibleJumpFromExon (zDNA *dna, zPhase_t phase, zSfeature *exon);
int zCompatibleJumpToExon (zDNA *dna, zPhase_t phase, zSfeature *exon);

score_t zScoreExternalState(struct zTrellis *trellis, int ext_state, 
							coor_t pos, coor_t length, score_t cscore);
score_t zScoreInternalState(struct zTrellis *trellis, int state, 
							coor_t pos, bool first_base);

#endif





