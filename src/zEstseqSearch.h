
/******************************************************************************\
 zEstseqSearch.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002

\******************************************************************************/

#ifndef ZOE_EstseqSearch_H
#define ZOE_EstseqSearch_H

#include <stdio.h>
#include <string.h>

#include "zScanner.h"
#include "zEstseq.h"
#include "zTools.h"

/******************************************************************************\
 zEstseqSearch

 A simple object which holds the distance from any position to the next/previous 
 occurance of a 0,1,2 in a zEstseq

\******************************************************************************/

struct zEstseqSearch {
  short id;
  zEstseq*  estseq;
  zList     value;
  short     symbols;
  coor_t*   temp;
};
typedef struct zEstseqSearch zEstseqSearch;

void zInitEstseqSearch (zEstseqSearch*,zEstseq*);
void zFreeEstseqSearch (zEstseqSearch*);

coor_t zGetEstseqPrevPos(zEstseqSearch*, char, coor_t);
coor_t zGetEstseqNextPos(zEstseqSearch*, char, coor_t);

coor_t zGetEstseqPrevPosBounded(zEstseqSearch*, char, coor_t, coor_t);
coor_t zGetEstseqNextPosBounded(zEstseqSearch*, char, coor_t, coor_t);
#endif
