/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zScanner.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_ALIGNMENTSCANNER_H
#define ZOE_ALIGNMENTSCANNER_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "zAlignment.h"
#include "zBNTreeModel.h"
#include "zMath.h"
#include "zSfeature.h"
#include "zTools.h"

/******************************************************************************\
 zAlignmentScanner

A scanner is the interface between an alignment (zAlignment) sequence and
a sequence model (zModel). The chief purpose of a scanner is that it allows 
one to ask the score of a feature at some position of the sequence. Scanners 
also provide a mechanism for defining specific scores at specific positions.

	int i;
	score_t s;
	zScanner scanner;
	zInitScanner(&scanner, &alignment, &model);
	
	for (i = 0; i < scanner.alignment.length) {
		score = scanner->score(&scanner, i);
	}

\******************************************************************************/

struct zAlignmentScanner {
	coor_t                min_pos;     /* minimum scoring position */
	coor_t                max_pos;     /* maximum scoring position */
	zDNA                  *dna;        /* target sequence */
	zAlignment            *alignment;  /* informant alignments */
	zBNTreeModel          *model;      /* scanners require a BNTree model */
	int                   subscanners;
	struct zAlignmentScanner *subscanner; /* for CDS */
	score_t               *uscore;      /* user-defined score */
	score_t               *uscore2;     /* user-defined sum score */
	zStrIdx               Einit, Exon, Eterm, Esngl;  /* for CDS */
	score_t (*score)(const struct zAlignmentScanner*, coor_t);
	score_t (*scoref)(const struct zAlignmentScanner*, zSfeature*);

	int id;
};
typedef struct zAlignmentScanner zAlignmentScanner;

void zFreeAlignmentScanner (zAlignmentScanner*);
void zInitAlignmentScanner (zAlignmentScanner*, zDNA*, zAlignment*, zBNTreeModel*);
void zSetAlignmentScannerScore (zAlignmentScanner*, coor_t, score_t);
void zPreComputeAlignmentScanner(zAlignmentScanner*);
score_t zGetUScoreA(const zAlignmentScanner *scanner,coor_t i);
score_t zGetUScoreARange(const zAlignmentScanner *scanner, coor_t start, coor_t end);
#endif
