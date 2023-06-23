/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zAlnFeature.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2004-2005 Manimozhiyan Arumugam

\******************************************************************************/

#ifndef ZOE_AFEATURE_H
#define ZOE_AFEATURE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zTools.h"
#include "zDNA.h"

/******************************************************************************\
 zAlnFeature

zAlnFeature is for alignment features like matches, gaps, etc. zAlnfeatureName is
a very important enumerated type.

The zAlnFVec copies value, so you can pass a pointer to a zAlnFeature in and it
will allocate storage for it. It is very important that you make sure that if
the group is unspecified, you give it a NULL value. You don't have to do this
for zReadAlnFeature, but if you're creating zAlnFeatures on the fly, be aware that
the zCopyAlnFeature command will segfault if the group points randomly.

\******************************************************************************/

struct zAlnFeature {
	zStrIdx   name;
	int       state;
	strand_t  strand;
	zDNA     *genomic;
	zStrIdx   genomic_def;
	coor_t    genomic_start;
	coor_t    genomic_end;
	zDNA     *cdna;
	zStrIdx   cdna_def;
	coor_t    cdna_start;
	coor_t    cdna_end;
	coor_t    length;
	score_t   score;
	score_t   percent_identity;
	coor_t    padding;
};

typedef struct zAlnFeature zAlnFeature;

struct zAFVec {
        zAlnFeature *elem; /* array of features */
        int        size;   /* number of features */
        int        limit;  /* number of elements currently allocated */
        zAlnFeature *last;
};
typedef struct zAFVec zAFVec;


/* zAlnFeature Stuff */

void zClearAlnFeature  (zAlnFeature*);
void zFreeAlnFeature   (zAlnFeature*);
int  zReadAlnFeature   (FILE*, zAlnFeature*, zDNA*, zDNA*);
void zWriteAlnFeature  (FILE*, const zAlnFeature*, int is_gtf, coor_t g_offset);
void zCopyAlnFeature   (const zAlnFeature*, zAlnFeature*);
int  zVerifyAlnFeature (const zAlnFeature*);
int  zAlnFeatureCmp    (const zAlnFeature*, const zAlnFeature*);
int  zAFPtrCmp         (const void*, const void*);
int  zAntiAlnFeature   (zAlnFeature*);

/**********************************************************************\
  zAFList - list of zAlnFeature objects 
\**********************************************************************/

struct zAFListNode;
typedef struct zAFListNode zAFListNode;
struct zAFListNode{
	zAlnFeature    data;
	zAFListNode* next;
	zAFListNode* prev;
};

struct zAFList {
	zAFListNode* head;
	zAFListNode* tail;
	zAFListNode* current;
	int          size; 
};
typedef struct zAFList zAFList;

void       zInitAFList(zAFList*);
void       zFreeAFList(zAFList*);
void       zResetAFList(zAFList*);

zAlnFeature* zAFListMoveNext(zAFList*);
zAlnFeature* zAFListMovePrev(zAFList*);
zAlnFeature* zAFListMoveFirst(zAFList*);
zAlnFeature* zAFListMoveLast(zAFList*);

zAlnFeature* zAFListGetCurrent(zAFList*);
zAlnFeature* zAFListInsert(zAFList*,zAlnFeature*);
zAlnFeature* zAFListInsertNext(zAFList*,zAlnFeature*);
zAlnFeature* zAFListInsertPrev(zAFList*,zAlnFeature*);
zAlnFeature* zAFListAppend(zAFList* l);
zAlnFeature* zAFListPrepend(zAFList* l);
void       zAFListRemoveFirst(zAFList*);
void       zAFListRemoveLast(zAFList*);
void       zAFListRemoveCurrent(zAFList*);
void       zAFList2AFVec(zAFList*,zAFVec*);
/* zAFVec Stuff */

void    zInitAFVec (zAFVec*, int);
void    zPushAFVec (zAFVec*, const zAlnFeature*);
void    zFreeAFVec (zAFVec*);
void    zFlipAFVec (zAFVec*);
void    zCopyAFVec (zAFVec* original, zAFVec* copy);
zAFVec* zReadAFVec (FILE*, zDNA*, zDNA*, int*, int*);
void    zWriteAFVec (FILE*, const zAFVec*, int is_gtf, coor_t g_offset);
void    zTranslateAFVec (zAFVec*, int);
void    zAFVec2GFF (const zAFVec* afv, zAFVec* gffv);
void    zInsertAFVec (zAFVec*, const zAlnFeature*, const int);
void    zAntiAFVec (zAFVec*);

/* Pretty Alignment */

int     zChopString(char *instring, int length);
void    zWriteAlignment (FILE*, const zAFVec*, coor_t g_offset);

#endif
