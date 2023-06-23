/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zModel.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_MODEL_H
#define ZOE_MODEL_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zMath.h"
#include "zTools.h"

/******************************************************************************\
 zModel

zModel is used for modeling sequence features like splice sites or
coding sequence. There are currently two primitive models (WMM and
LUT) and 5 compound models (SDT, SAM, CDS, MIX, and ISO).

IO Format:

Each model starts with a header line:

<name_str> <model_type> <DNA/Conseq> <length_int> <focus_int> <num_symbols> <num_submodels>

For example:

Start SDT DNA 3 0 4 2

In addition, the ISO and MIX models take additional parameters, as described below.

WMM: Weight Matrix Model.  Nucleotide positional odds (Not dinucleotide)

LUT: Lookup Table, or nth-order Markov Model.  

SDT: Sequence Decision Tree.  Decides which submodel to use based on
sequence. For example, an acceptor SDT might make a decision on
encountering the sequence "ag" to score with a WMM, but use LUT that
scores -Inf for "nn"

SAM: Sequence Array Model.  Scores the sequence under all submodels,
then multiplies the probability together.

CDS: Coding Sequence.  3 periodic with LUT submodels.

MIX: Mixture.  Weighted average of probabilities (not an average of
scores).  For each submodel, and additional parameter specifying the
weight mustbe provided in the model header. For example:

    Acceptor MIX DNA 0 0 4 2       0.7  0.3
        Acc1 LUT DNA 3 3 4 0
		   ...
		Acc2 WMM DNA 9 6 4 0
           ...

For a given sequence position, if x1 is the probability for Acc1 and x2 is the probablity for Acc2, then the Acceptor model will return the score corresponding to the probability 0.7 * x1 + 0.3 * x2.

ISO: Isochore.  Selects a submodel based on GC content.  If a per-base
GC content has been calculated for the sequence, then that is used.
Otherwise the overall GC content is used.  The header line for an
isochore model should include one extra parameter per submodel,
specifying the upper GC level for each sub model.  For example, under
the model

    Acceptor ISO DNA 0 0 4 4    40.0 55.0 68.0 100.0
	  Acc1 ...
	  Acc2 ...
	  Acc3 ...
	  Acc4 ...
	  
Acc1 is used for a GC level between 0% and 40%, Acc2 for GC levels 40% -
55%, Acc3 for GC levels 55% - 68%, and Acc4 for GC levels 68% - 100%.



\******************************************************************************/

typedef enum {
	WMM,  /* Weight Matrix Model */
	LUT,  /* Lookup Table (aka Markov Model) */
	WWAM, /* Windowed Weight Array Matrix */
	SAM,  /* Scoring Array Model */
	SDT,  /* Sequence Decision Tree */
	CDS,  /* 3 Periodic (LUT submodels) */
	MIX,  /* averages all submodels */
	ISO,  /* Isochore model */
    SIG  /* Codon level WMM (used for signal peptide) */
} zModelType;

typedef enum {
	DNA,
	CONSEQ,
	ESTSEQ,
	GENOMIC,
	PAIR
} zSeqType;

struct zModel {
	zModelType     type;        /* enumerated above */
	char*          name;        /* enumerated above, but also used by SDT */
	coor_t         length;      /* length of the model; how many bp covered */
	coor_t         focus;       /* which bp of model gets score */
	int            symbols;     /* number of symbols in the alphabet */
	int            order;       /* for WWAM */
	int            submodels;   /* number of sub-models */
	coor_t         window_size; /*for windowed GC with ISO models*/
	struct zModel *submodel;    /* sub-models, for complex models */
	score_t       *data;        /* scores */
	zSeqType       seq_type;
};
typedef struct zModel zModel;

void zFreeModel (zModel*);
int  zReadModel (FILE*, zModel*);
int  zReadModelHeader (FILE*, zModel*, int);
void zWriteModel (FILE*, const zModel*);
void zAmbiguateModel (zModel*, score_t);
void zDeambiguateModel (zModel*);
int zAntiModel (zModel* model); 
void zMarginalizeModel (zModel *model, int base); 

#endif
