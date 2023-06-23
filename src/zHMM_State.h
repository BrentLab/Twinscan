/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zHMM_State.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_HMM_STATE_H
#define ZOE_HMM_STATE_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zSfeature.h"
#include "zTools.h"
#include "zFeatureFactory.h"
#include "zTransition.h"

/******************************************************************************\
 zHMM_State

Zoe recognizes two basic types of HMM states: External (Burge's
C-states) and Internal (Burge's D-states).  In addition, Zoe can
expand a single specified state into several internal states that
track phase, with the Tracking type of state.

Stream Format
=============

There are seven parameters to a state:

Name  Type  Strand Phase Duration_Name Model_name init_probability

Where: Name      - is a string
       Type      - is {Internal, External, Tracking}
       Strand    - is {+, -}
	   Phase     - is {0, 1, 1T, 2, 2TA, 2TG}
	   Duration  - is a string that matches a duration's name
	   Model     - is a string that matches a sequenc model's name
	   init      - is a probabilit (between 0.0 and 1.0, inclusive)

\******************************************************************************/

typedef enum {
	INTERNAL,       /* Traditional type HMM state type */
	GINTERNAL,      /* Genscan style  (uses intergenic continue score) */
	EXTERNAL,		/* Explicit duration state type */
	EXPLICIT
} zHMM_StateType; /* Hidden dependancy: transLookup in zTransition.c */

struct zHMM_State {
	zHMM_StateType   type;      /* enumerated above */
	zStrIdx          name;
	strand_t         strand;
	zPhase_t         phase;

	zStrIdx          duration; /* which duration to use */
	zStrIdx          model;    /* which model to use */

	zFFInit          ffactory;
	
	score_t*          init;      /* initial score of state indexed by the isochore group*/
                                     /* init[0] is the initial score of the state in isochore group 0 */
};
typedef struct zHMM_State zHMM_State;

int  zReadHMM_State (FILE*, zHMM_State*, int);
void zWriteHMM_State (FILE*, const zHMM_State*, int);

#endif
