/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zTransition.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_PAIR_TRANSITION_H
#define ZOE_PAIR_TRANSITION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zMath.h"
#include "zSfeature.h"
#include "zTools.h"

struct zPairTrellis;
struct zScanner;
typedef void (*zPairTransFunc)(struct zPairTrellis* trellis, int from, int to, coor_t genomic, coor_t cdna);

zPairTransFunc zGetPairTransFunc(int type); 
zPairTransFunc zGetPairBackTransFunc(int type);

int zCompatibleJumpFromExon (zDNA *dna, zPhase_t phase, zSfeature *exon);
int zCompatibleJumpToExon (zDNA *dna, zPhase_t phase, zSfeature *exon);

score_t zGetScannerScore(struct zPairTrellis *trellis, struct zScanner *scanner, int state, coor_t genomic, coor_t cdna);

void zInternalPairTransHelper (struct zPairTrellis *trellis, int from_state, int to_state, coor_t genomic, coor_t cdna, score_t trans);
#endif





