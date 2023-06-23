/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zPairTrellis.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_PAIR_TRELLIS_H
#define ZOE_PAIR_TRELLIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zAlnFeature.h"
#include "zDNA.h"
#include "zFastaFile.h"
#include "zFeatureFactory.h"
#include "zHMM.h"
#include "zModel.h"
#include "zPairTransition.h"
#include "zScanner.h"
#include "zSeedUtils.h"
#include "zSfeature.h"
#include "zStopSeq.h"
#include "zTools.h"

struct zPairTrellisCell {
	/* Regular Viterbi */
	score_t       score; 
	int           length:24;
	int           trace:7;  /* 0 state index */
	unsigned int  keep:1;
};
typedef struct zPairTrellisCell zPairTrellisCell;


/******************************************************************************\
 zPairTrellis

zPairTrellis collects many of the zoe components into a single entity used for the
gene prediciton algorithms. A trellis contains a specific zDNA and zHMM.

\******************************************************************************/

struct zPairTrellis {
	zDNA             *genomic;
	zDNA             *cdna; /* Placeholder for cdna sequence: no memory allocation involved */
	zDNA             *fcdna; /* cdna */
	zDNA             *rcdna; /* reverse complemented cdna */

	zHMM             *hmm;
	int               tiso_group; /* isocore group for transitions */
	int               iiso_group; /* isocore group for init. prob. */
	score_t           seq_prob;  /* result of forward alg */
	score_t           forward_score;
	score_t           backward_score;
	coor_t            padding;
	
	/* members used for optimized memory usage with blocks from seed alignments */

	zSeedAlignment   *seed;
	zSeedAlignment   *blocks;
	zSeedAlignment   *mem_blocks;
	int               allocated_mem_blocks;
	int               allocated_fwd_blocks;
	int               allocated_bak_blocks;

	/* indexed by pos, pos and hmm state */
	zPairTrellisCell  ***cell;      /* cell[state_idx][seq_pos] */

	score_t         ***forward;   /* forward algorithm values */
	score_t         ***backward;  /* backward algorithm values */
	zSFList        **external;   /* array of lists of external feature on 
									+ strand ending at current pos */

	/* indexed by zStrIdx */
	zScanner        **scanner;   /* map hmm models to scanners here */

	/* These are not used, but I left them there so someone might use them later */

	zFeatureFactory **factory;   /* external feature factory */
	zIVec           **fmap5;
	zIVec           **fmap3;
	zAFVec          **fexternal; /* external forward links */
	zIVec            *extpos;     /* backward links to external       */
                                     /* features scoring above MIN_SCORE */
                                     /* This is used to compute explicit */
                                     /* state features                   */

	zStopSeq         *stopseq;
};
typedef struct zPairTrellis zPairTrellis;

/*********************************************\
 zPairTrellis Utilities
\*********************************************/

void    zInitPairTrellis (zPairTrellis*, zSeedAlignment*, zDNA*, zDNA*, zHMM*);
void    zFreePairTrellis (zPairTrellis*);

/*********************************************\
 Regular Viterbi Decoding
\*********************************************/

zAFVec* zRunPairViterbiAndForward (zPairTrellis*,score_t*);
void    zTracePairTrellis (zPairTrellis*, int, zAFVec*);
void    zRunPairBackward (zPairTrellis *trellis);
void    zComputePosteriorProbability(zPairTrellis* trellis, zAFVec* path);

/*********************************************\
 Functions for Evan's Memory Optimized Version
\*********************************************/

int     zFindPairTrellisPin(zPairTrellis*,coor_t);
zAFVec* zRunPairViterbi(zPairTrellis*, score_t*);
zPtrList* zRunSNPPairViterbi(zPairTrellis*, score_t*);
zSFVec* zRunPinPairViterbi(zPairTrellis*, score_t*, coor_t, coor_t);

/*********************************************\
 General Utilities
\*********************************************/

void    ShowPairTrellis(zPairTrellis*);
void    zShowPairTrellisCell(zPairTrellis *trellis, coor_t i, coor_t j, int k); 
score_t zGetPairPathScore(zPairTrellis*, const zAFVec*, zAFVec*);

/*********************************************\
 zPairTrellis[Cell] Utilities
\*********************************************/
#ifdef DEBUG

zPairTrellisCell* zGetCurrentCell(zPairTrellis *trellis, coor_t genomic, coor_t cdna, int state); 
zPairTrellisCell* zGetPreviousCell(zPairTrellis *trellis, coor_t genomic, coor_t cdna, int current_state, int previous_state); 
zPairTrellisCell* zGetNextCell(zPairTrellis *trellis, coor_t genomic, coor_t cdna, int next_state); 
zPairTrellisCell* zGetCellArray(zPairTrellis *trellis, coor_t genomic, coor_t cdna); 

#else

#define      zGetCurrentCell(trellis, genomic, cdna, current_state)                 (&trellis->cell[genomic][cdna][current_state]) 
#define        zGetCellArray(trellis, genomic, cdna)                                (trellis->cell[genomic][cdna]) 
#define         zGetNextCell(trellis, genomic, cdna, next_state)                    (&trellis->cell[genomic+zGetGenomicIncrement(trellis->hmm, next_state)][cdna+zGetCDnaIncrement(trellis->hmm, next_state)][next_state])
#define     zGetPreviousCell(trellis, genomic, cdna, current_state, previous_state) (&trellis->cell[genomic-zGetGenomicIncrement(trellis->hmm, current_state)][cdna-zGetCDnaIncrement(trellis->hmm, current_state)][previous_state])

#endif
#endif /* ZOE_PAIR_TRELLIS_H */

