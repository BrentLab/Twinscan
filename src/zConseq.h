/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zConseq.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_Conseq_H
#define ZOE_Conseq_H

#include <stdio.h>
#include <string.h>

#include "zSequence.h"
#include "zFastaFile.h"
#include "zProtein.h"
#include "zSfeature.h"
#include "zTools.h"

/******************************************************************************\
 zConseq

zConseq is basically a striped down zDNA type.  It holds conservation sequence
the same way the zDNA object holds DNA sequence and performs some functions
common to both and in the future will provide some functions specific to 
conservation sequence.

\******************************************************************************/

struct zConseq {
	char      *def;
	zSequence *seq;
	int        bits;
	int        digits;
	coor_t     length;
};
typedef struct zConseq zConseq;

void zFreeConseq (zConseq*);
void zInitConseq (zConseq*, int);
void zCopyConseq (const zConseq*, zConseq*);
void zLoadConseqFromFasta(zConseq*, char*);
void zReverseConseq (zConseq*);
void zSetConseqPadding (zConseq*,coor_t);
char zGetConseqSeq(zConseq*,coor_t);
char zGetConseqS10(zConseq*,coor_t);
#endif
