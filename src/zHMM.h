/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zHMM.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_HMM_H
#define ZOE_HMM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zAlnFeature.h"
#include "zBNTreeModel.h"
#include "zDuration.h"
#include "zHMM_State.h"
#include "zMath.h"
#include "zModel.h"
#include "zPhasePref.h"
#include "zTransition.h"
#include "zTools.h"
#include "zGTF.h"

/******************************************************************************\
 zHMM

zHMM integrates zHMM_States, zTransitions between those states,
zDurations of those states, and the sequence models.

zHMM Files
==========

A zHMM is read in from a file that has 6 sections: a header, states,
transitions, phase preferences, state durations, and sequence models.
An example can be found at zoe/T/iscan.zhmm.

The header looks like:
  zHMM <name> <# states> <# transitions> <# durations> <# models>
An example would be:
  zHMM my_hmm 4 5 4 2

Each of the remaining sections is introduced by the appropriate tag:
<STATES>, <STATE_TRANSITIONS>, <PHASE_PREFERENCES>,
<STATE_DURATIONS>, and <SEQUENCE_MODELS>.  After each tag, zHMM
expects exactly the number of objects specified in the header.  The
file format of each type of object is documented in each objects
header file.

\******************************************************************************/

typedef enum {GHMM, GPAIRHMM} hmm_mode_t;

struct zHMM {

	/* general attributes */

	hmm_mode_t  mode;          /* = GHMM or GPAIRHMM */
	char       *name;          /* completely arbitrary */
	int         orig_states;   /* number of states in specification */
	int         states;        /* number of states in the HMM */
	int         transitions;   /* number of transtions in the HMM */
	int         durations;     /* number of duration models (zDurationGroup structs) */
	int         models;        /* number of sequence models */
	int         cons_models;   /* number of conservation sequence models */
	int         est_models;    /* number of conservation sequence models */
	int         phylo_models;  /* number of phylogenetic models for alignment*/

	int         feature_count;  /* number of zStrIdx for HMM */
	int         iso_states;     /* number of isochore groups for states*/
	int         iso_transitions;/* number of isochore groups for transitions*/
	float      *iso_state;      /* isochore group boundaries for states*/
	float      *iso_transition; /* isochore group boundaries for transitions*/
	
	/* object storage ------------------------------------------------ */

	zHMM_State     *orig_state;
	zHMM_State     *state;      /* array of zState */
	zTransition    *transition; /* array of zTransition */
	zDurationGroup *duration;   /* array of zDuration */
	zModel         *model;      /* array of zModel */
	zModel         *cons_model; /* array of zModel */
	zModel         *est_model;  /* array of zModel */
	zBNTreeModel   *phylo_model;/* array of zBNTreeModel */

/* 	zPhasePref      phasepref;  /\* exon-intron phase preferences *\/ */
	zGTFConvInfo   *gtf_conv;   /* information on how to convert to GTF, only for GHMM */
	int           **increments; /* array of genomic/cdna coordinate increments per state, only for GPAIRHMM */
	
	/* object mapping ------------------------------------------------ */

	zIVec       *simap; /* origstate index --> state index */

	/*  -- indices here are zStrIdx and < feature_count */
	int              *somap; /* original states */
	zStrIdx          *reverse_somap;
	zModel          **mmap;  /* models */
	zModel          **cmmap; /* conservation models */
	zModel          **emmap; /* est  models */
	zBNTreeModel    **ammap; /* alignment models */
	zDurationGroup  **dmap;  /* durations */

	/*  -- indices here correspond to the indices into state[] */
	zIVec      **jmap;  /* jump list (reverse arrows) */
	zIVec      **fmap;  /* forward list (forward arrows) */
	score_t   ***tmap;  /* transition score tmap[from][to][iso_group] */
	score_t     *inter_continue; /* intergenic continue scores */

};
typedef struct zHMM zHMM;

void zFreeHMM (zHMM*);
int  zReadHMM (FILE*, zHMM*, hmm_mode_t);
void zWriteHMM (FILE*, const zHMM*);
int zMapHMM (zHMM *hmm);

int zGetStateIdxFromStrIdx(const zHMM*, zStrIdx);
int zFillInStateIdxFromFeature(const zHMM* hmm, zDNA* dna, zDNA* rdna, zSfeature *f);

void zGTFVec2SFVec(const zHMM*, zDNA*, zDNA* rdna, zGTFVec*, zSFVec*);
void zSFVec2GTFVec(const zHMM*, const zSFVec*, zGTFVec*, char*);

bool zUsedInExplicitState(const zHMM *, zStrIdx);

score_t zGetInitProb(zHMM*, int, int);
score_t zGetFixedInitProb(zHMM*, int, int, float);

/* Specific to GPAIRHMM */

int zNullifyHMM(zHMM *hmm); 
int zFixInternalTransitions(zHMM* hmm);

#ifdef DEBUG

int     zGetGenomicIncrement(const zHMM*, int);
int     zGetCDnaIncrement(const zHMM*, int);
score_t zGetTransitionScore(zHMM*, int, int, int);

#else

#define zGetGenomicIncrement(hmm, state) (hmm->increments[state][0])
#define    zGetCDnaIncrement(hmm, state) (hmm->increments[state][1])
#define  zGetTransitionScore(hmm, from_state, to_state, iso_group) (hmm->tmap[from_state][to_state][iso_group])

#endif


#endif
