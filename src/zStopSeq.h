
/******************************************************************************\
 zStopSeq.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002

\******************************************************************************/

#ifndef ZOE_StopSeq_H
#define ZOE_StopSeq_H

#include <stdio.h>
#include <string.h>

#include "zScanner.h"
#include "zDNA.h"
#include "zTools.h"

/******************************************************************************\
 zStopSeq

 A simple object which holds the distance from any position to the previous in 
 frame stop codon based on a DNA object.

\******************************************************************************/

struct zStopSeq {
  short id;
  zDNA*     dna;
  zList     value;
  zList     rev_value;
};
typedef struct zStopSeq zStopSeq;

void zInitStopSeq (zStopSeq*,zDNA*);
void zFreeStopSeq (zStopSeq*);

coor_t zGetStopSeqPrevPos(zStopSeq*, coor_t);
coor_t zGetStopSeqNextPos(zStopSeq*, coor_t);

coor_t zGetStopSeqMinPrevPos(zStopSeq*, coor_t);
coor_t zGetStopSeqMaxNextPos(zStopSeq*, coor_t);
#endif
