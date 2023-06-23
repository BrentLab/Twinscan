/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zScanner.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_SCANNER_H
#define ZOE_SCANNER_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zMath.h"
#include "zModel.h"
#include "zSfeature.h"
#include "zTools.h"
#include "zSequence.h"
#include "zEstseq.h"
#include "zConseq.h"
#include "zDNA.h"

/******************************************************************************\
 zScanner

 A scanner is the interface between a sequence (zSequence) and a model
 (zModel). The chief purpose of a scanner is that it allows one to ask the score
 of a feature at some position of the sequence. 

 The zScanner object is organized to hold a list of uscore blocks, each of which 
 hold the scores for a certain range of the scanner's sequence.  In general the
 scanner will have fewer blocks (zUScoreBlock) than it needs to hold the score 
 for the full sequence.  If a score is requested and that region is not currently
 loaded into a block then sequence range to load is chosen from the uscore_map 
 and the scores are loaded into the least recently used block.

 This scanner replaces the old zScanner, zConseqScanner, zESTScanner, 
 zAlignmentScanner, etc. because it is general to deal with scanners for any
 type of sequence.

\******************************************************************************/

/* the data for scanner scores corresponding to some block of sequence */
struct zUScoreBlock {
	score_t*  score;     /* array of scores for this block */
	score_t*  score_sum; /* array of sum of scores for this block */
	int*      var_map;   /* specify the index of the scannervariant 
							for this position */
    coor_t    pos;       /* position of the first base scored in score */
	int       map_idx;   /* index in uscore_map of this block */
	int       idx;       /* index in uscore of this block */
};
typedef struct zUScoreBlock zUScoreBlock;

/* The data for a  group of overlapping (in this scanners application)
   sequence variants */
struct zScannerVariant {
	zSeqVariant**     vars;       /* list of varying poisitions */
	int               var_count;  /* number of varying positions */
	score_t**         score;      /* scores for each posible value of the 
									 varying positions 
									 score[i][j] : i=var_val index
									 j=position index */
	score_t**         score_sum;  /* sum scores for each posible value of the 
									 varying positions 
									 score_sum[i][j] : i=var_val index
									 j=position index */
	int               var_vals;   /* total number of possible variant
									 value combinations (e.g. if each 3 vars
									 have 3 possible values, then this is
									 27 */
	coor_t            first_pos;  /* the first position covered by these 
									 scores */
	coor_t            last_pos;   /* the last position covered by theses
									 scores */
};
typedef struct zScannerVariant zScannerVariant;

/* Mapping between sequence coordinates and a score block */ 
struct zUScoreMap {
	coor_t            pos;        /* first position of this map */
	int               block_idx;  /* index in uscore of block holding scores
									 for this map.  -1 if not loaded */
	zScannerVariant** vars;       /* array of zScannerVariants for this map */
	short             var_count;  /* number of varying position blocks in this
									 map */
	short             tag;        /* tag value used in zScoreCoding */
};	
typedef struct zUScoreMap zUScoreMap;

struct zScanner {
	coor_t            min_pos;    /* minimum scoring position */
	coor_t            max_pos;    /* maximum scoring position */
	coor_t            min_gpos;    /* minimum scoring position in genomic for PAIR model scanners */
	zSequence        *seq;        /* scanners require a seq */
	zModel           *model;      /* scanners require a model */
	struct zScanner  *subscanner; /* for higher order models */
	struct zScanner  *parent;     /* parent scanner */
	char             *sig;        /* binary signature (SDT) */
	zUScoreMap       *uscore_map; /* the array storing the map from
									 the sequence to the scanner blocks */
	int               uscore_map_size; /* the number of zUscoreMap objects */
	zUScoreBlock     *uscore;     /* an array of all blocks */
	zPtrList          uscore_list; /* the list of all blocks */
                   	  /* this list is used to keep an order on 
						 the blocks */
	int               uscore_block_count;/* number of blocks */
	coor_t            uscore_block_size; /* number of bases covered
											by a single block */

	zPtrList         *vars;       

	zStrIdx          Einit, Exon, Eterm, Esngl; /* for CDS */
	
	score_t (*score)(struct zScanner*, coor_t); /* score a base */
	score_t (*scoref)(struct zScanner*, zSfeature*); /* score a feature */
	score_t (*scorer)(struct zScanner*, coor_t, coor_t); /* score a range of bases */

	score_t (*pairscore)(struct zScanner*, zDNA*, zDNA*, coor_t, coor_t); /* score a pair of base */

	long int        *temp_array;   

	int id;
};
typedef struct zScanner zScanner;

void zFreeScanner (zScanner*);
void zInitScanner (zScanner*, zSequence*, zModel*);
void zSetScannerScore (zScanner*, coor_t, score_t);
score_t zScoreFeature (zScanner*, zSfeature*);
void zPreComputeScanner (zScanner*);
void zDePreComputeScanner (zScanner*);
score_t zScoreSIG (zScanner *scanner, zSfeature *f);
int zGetScannerIso(zScanner *scanner, float gc);
zScanner* zGetScannerForIso(zScanner *scanner, float gc);
score_t zGetUScore(zScanner *scanner,coor_t i);
score_t zGetRangeScore(zScanner *scanner,coor_t start, coor_t end, int frame);
#endif



