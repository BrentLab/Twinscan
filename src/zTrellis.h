/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zTrellis.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_TRELLIS_H
#define ZOE_TRELLIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zDNA.h"
#include "zFeatureFactory.h"
#include "zHMM.h"
#include "zModel.h"
#include "zScanner.h"
#include "zSfeature.h"
#include "zTools.h"
#include "zConseq.h"
#include "zEstseq.h"
#include "zEstseqSearch.h"
#include "zAlignmentScanner.h"
#include "zStopSeq.h"

struct zTBTreeNode;

struct zTrellisCell {
	score_t       score; 
	int           length;
	int           trace;  /* 0 state index */
  struct zTBTreeNode*  tb; 
	zSfeature*    jump;	  

	frame_t       frag_frame;   /* keep lfrag, rfrag and frame   */
	/* information in bitwise manner */
	
	zSFVec        *submax; /* for states of EXPLICIT and INTERNAL   */
	/* types. Keeps highest-scoring features */ 
	/* for each range of length distribution */
	/* see zExplicitTrans function in        */
	/* zTransition.c for details             */

};
typedef struct zTrellisCell zTrellisCell;

/******************************************************************************\
 zTrellis

zTrellis collects many of the zoe components into a single entity used for the
gene prediciton algorithms. A trellis contains a specific zDNA and zHMM.

\******************************************************************************/

struct zTrellis {
	zDNA             *dna;
	zDNA             *rdna; /* reverse complemented dna */
	zDNA             *unmasked_dna;
	zDNA             *unmasked_rdna; /* reverse complemented dna */
	zConseq          *conseq;
	zConseq          *rconseq;
	zEstseq          *estseq;
	zEstseq          *restseq;
	zAlignment       *alignment;
	zAlignment       *ralignment;
	bool              est_para_mode;
	bool              phylo_enabled;

#ifdef DEBUG
	bool              show_trellis; /* indicates whether --show_trellis */
	                                /* option has been used             */
#endif

	zHMM             *hmm;
	int               tiso_group; /* isocore group for transitions */
	int               iiso_group; /* isocore group for init. prob. */
	score_t           seq_prob;  /* result of forward alg */
	coor_t            padding;
	
	/* indexed by an hmm state */
	zSFVec          **fexternal; /* external forward links */
	zTrellisCell    **cell;      /* cell[state_idx][seq_pos] */
	score_t         **forward;   /* forward algorithm values */
	score_t         **backward;  /* backward algorithm values */
	zIVec           *extpos;     /* backward links to external       */
                                     /* features scoring above MIN_SCORE */
                                     /* This is used to compute explicit */
                                     /* state features                   */
	zSFList        **external;   /* array of lists of external feature on 
									+ strand ending at current pos */
	zSFList        **rexternal;  /* array of lists of external feature on 
									- strand ending at current pos */

	/* indexed by zStrIdx */
	zScanner        **scanner;   /* map hmm models to scanners here */
	zScanner        **rscanner;  /* scanners for reverse dna */
	zScanner        **conseqscanner;
	zScanner        **rconseqscanner;
	zScanner        **estseqscanner;
	zScanner        **restseqscanner;

	zAlignmentScanner **alignmentscanner;
	zAlignmentScanner **ralignmentscanner;


	zFeatureFactory **factory;   /* external feature factory */
	zFeatureFactory **rfactory;  /* external feature factory */
	zIVec           **fmap5;
	zIVec           **rfmap5;
	zIVec           **fmap3;
	zIVec           **rfmap3;

	zStopSeq         *stopseq;
	zStopSeq         *rstopseq;
	
	zEstseqSearch    *estsearch;
	zEstseqSearch    *restsearch;
};
typedef struct zTrellis zTrellis;

void    zFreeViterbiVars (zTrellis*);
void    zFreeForwardVars (zTrellis*);
void    zFreeBackwardVars (zTrellis*);
void    zFreeBackLinks (zTrellis*);
void    zFreeForwardLinks (zTrellis*);
void    zFreeTrellis (zTrellis*);
void    zInitTrellis (zTrellis*, char*, zHMM*, char*, char*, char*, bool, char*, char*);
zSFVec* zRunViterbiAndForward (zTrellis*,score_t*);
void    zRunPartialViterbiAndForward (zTrellis*, coor_t, coor_t);
void    zTraceTrellis (zTrellis*, int, zSFVec*);
void    zTracePartialTrellis (zTrellis*, coor_t, coor_t, int, zSFVec*);
void    zRunBackward (zTrellis*);
float   zForwardBackwardProb(zTrellis*, zSfeature*);
void    ShowTrellis(zTrellis*);

score_t zGetIntergenicContinueScore(zTrellis*);
score_t zGetPathScore(zTrellis*, const zSFVec*, zSFVec*);

frame_t zFragFrame2Char(frame_t, frame_t, frame_t);
void zChar2FragFrame(frame_t, zSfeature *);


int zFindTrellisPin(zTrellis*,coor_t);

zSFVec* zRunViterbi(zTrellis*, score_t*);
zPtrList* zRunSNPViterbi(zTrellis*, score_t*);
zSFVec* zRunPinViterbi(zTrellis*, score_t*, coor_t, coor_t);

#endif

