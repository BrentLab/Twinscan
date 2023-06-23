/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zPhasePref.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_PHASEPREF_H
#define ZOE_PHASEPREF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zMath.h"
#include "zSfeature.h"
#include "zTools.h"

/******************************************************************************\
 zPhasePref

Phase preferences are exon-intron and intron-exon probabilities (scores
actually).

\******************************************************************************/

#define PHASECOUNT 18

typedef enum {
	Ei_I0, /* transition from Einit to Int0 */
	Ei_I1, /* transition from Einit to Int1 */
	Ei_I2, /* transition from Einit to Int2 */
	E0_I0, /* transition from Exon0 to Int0 */
	E0_I1, /* transition from Exon0 to Int1 */
	E0_I2, /* transition from Exon0 to Int2 */
	E1_I0, /* transition from Exon1 to Int0 */
	E1_I1, /* transition from Exon1 to Int1 */
	E1_I2, /* transition from Exon1 to Int2 */
	E2_I0, /* transition from Exon2 to Int0 */
	E2_I1, /* transition from Exon2 to Int1 */
	E2_I2, /* transition from Exon2 to Int2 */
	
	I0_E0, /* transition from Int0 to Exon0 */
	I0_Et, /* transition from Int0 to Eterm */
	I1_E1, /* transition from Int1 to Exon1 */
	I1_Et, /* transition from Int1 to Eterm */
	I2_E2, /* transition from Int2 to Exon2 */
	I2_Et  /* transition from Int2 to Eterm */
} zPhasePrefName;

struct zPhasePref {
	score_t score[PHASECOUNT];
	float   prob[PHASECOUNT];

	zStrIdx Esngl, Einit, Exon, Eterm;
};
typedef struct zPhasePref zPhasePref;

int     zReadPhasePref (FILE*, zPhasePref*);
void    zWritePhasePref (FILE*, const zPhasePref*);
score_t zScorePhaseFromE (zPhasePref*, zStrIdx, zPhase_t, int);
score_t zScorePhaseToE (zPhasePref*, zPhase_t, zStrIdx);

#endif
