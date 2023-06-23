/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zEstseq.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_Estseq_H
#define ZOE_Estseq_H

#include <stdio.h>
#include <string.h>

#include "zSequence.h"
#include "zFastaFile.h"
#include "zProtein.h"
#include "zSfeature.h"
#include "zTools.h"

/******************************************************************************\
 zEstseq

zEstseq is basically a striped down zDNA type.  It holds conservation sequence
the same way the zDNA object holds DNA sequence and performs some functions
common to both and in the future will provide some functions specific to 
conservation sequence.

\******************************************************************************/

struct zEstseq {
	char      *def;
	zSequence *seq;
	int        bits;
	int        digits;
	coor_t     length;
};
typedef struct zEstseq zEstseq;

void zFreeEstseq (zEstseq*);
void zInitEstseq (zEstseq*, int);
void zCopyEstseq (const zEstseq*, zEstseq*);
void zLoadEstseqFromFasta(zEstseq*, char*);
void zReverseEstseq (zEstseq*);
void zSetEstseqPadding (zEstseq*,coor_t);
char zGetEstseqSeq(zEstseq*,coor_t);
char zGetEstseqS10(zEstseq*,coor_t);
#endif
