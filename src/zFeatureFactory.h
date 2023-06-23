/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zFeatureFactory.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_FEATUREFACTORY_H
#define ZOE_FEATUREFACTORY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zScanner.h"
#include "zSfeature.h"
#include "zTools.h"
#include "zEstseq.h"
#include "zEstseqSearch.h"
#include "zAlignmentScanner.h"
#include "zStopSeq.h"

/* phylogenetic score coefficient */
/* multiplies all phylogenetic scores */
#define PSC 0.4f

/******************************************************************************\
 zFeatureFactory

Feature factories create sequence features like exons. There are several kinds
of factories. zEFactory is used only for exons (coding exons really). It
requires 5 specific scanners. Exon scores are caluclated by adding the score at
each end of the exon (start or acceptor + stop or donor). The zRFactory creates
repeats where there are long enough runs of N's. Repeat scores are all given a
constant score. The zSFactory is based on a scanner. The score comes directly
from the scanner, and the coordinates of the feature are based on the length
and offsets of the underlying model.

Examples:

\******************************************************************************/

struct zFeatureFactory;

typedef void (*zFFInit)(struct zFeatureFactory*, zDNA* dna, zEstseq* estseq,
						zScanner** scanners, zScanner** conseqscanners,
						zScanner** estseqscanners,
						zAlignmentScanner **alignmentscanners,
						int feature_count, strand_t strand, 
						bool conseq_enabled, bool est_para_mode,
						bool phylo_enabled, zStopSeq* stop_seq,
						zEstseqSearch* estsearch);
typedef void (*zFFCreate)(const struct zFeatureFactory*, coor_t, zSFList*);
typedef score_t (*zFFScoref)(const struct zFeatureFactory*, zSfeature*);

struct zFeatureFactory {
	/* used by all factories */
	zFFCreate create5;	/* create features starting (5' coord) at pos */
	zFFCreate create3;	/* create features ending (3' coord) at pos */
	zFFCreate createi;	/* create features crossing pos */
	zFFScoref scoref;
	zStrIdx   type;
	strand_t strand;	

	/* these are used by SFactories */
	int         length;  /* length of feature - 0 means variable */
	zScanner   *scanner; /* for those factories that use a scanner */
	
	/* these are used for RFactory */
	zDNA       *dna;	 /* for finding repeats */
	int         min_len;  /* minimum length */
	score_t     score;
	
	/* these are used for EFactory */
	coor_t      offset;
	int         end_offset; /* position offset for create functions */
	zStopSeq   *fstop;
	zEstseqSearch   *estsearch;
	
	/* New additions */
	zEstseq    *estseq;
	bool        est_para_mode;
	zScanner  **estscan;

	/* For phylogenetic scoring */
	zAlignmentScanner **alignscan;
	bool phylo_enabled;

	zScanner  **scan;
	zScanner **conscan;
	bool conseq_enabled;
};
typedef struct zFeatureFactory zFeatureFactory;

void zFreeFeatureFactory (zFeatureFactory*);

void     zRegisterFactory(const char* feature_name, zFFInit);
zFFInit  zGetFactory(const char* feature_name);
zScanner* zGetScanner( const zFeatureFactory *fac, int type); 
int not_in_intron_check(const zFeatureFactory *fac, coor_t start, coor_t end, zStrIdx name);
int overlap_exon_or_intron(const zFeatureFactory *fac, coor_t start, coor_t end);

#endif

