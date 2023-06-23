/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zGTF.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2003 Charles J. Vaske

\******************************************************************************/

#ifndef ZOE_GTF_H
#define ZOE_GTF_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zSfeature.h"

/******************************************************************************\
 zGTF

zGTF contains functions for working with the GTF format, which is specified at

http://genes.cs.wustl.edu/GTF2.html

Right now, the zGTFfeature structure exists primarily as a placeholder for
future functionality, and is immediately converted to a zSfeature.

\******************************************************************************/

typedef enum {
    zGTF_START, /* "start_codon" */
    zGTF_STOP,  /* "stop_codon"  */
    zGTF_CDS,   /* "CDS"         */
    zGTF_EXON,  /* "exon"        */
	zGTF_ESNGL, /* "zEsngl"      */
	zGTF_5UTR,  /* "5UTR"        */
    zGTF_UNDEF  /* anything else */
} gtf_t;

gtf_t  zText2GTF(const char*);
void   zGTF2Text(const gtf_t ftype, char* s);

struct zGTFfeature {
    char          *seqname;
    char          *source;
    gtf_t          type;
    coor_t         start;
    coor_t         end;
    score_t        score;
    strand_t       strand;
    frame_t        frame;
    char          *gene_id;
    char          *transcript_id;
};
typedef struct zGTFfeature zGTFfeature;

int  zReadGTFfeature (FILE *stream, zGTFfeature *f);
void zClearGTFfeature (zGTFfeature *f);
void zFreeGTFfeature (zGTFfeature *f);
void zWriteGTFfeature(FILE *stream, zGTFfeature *f);
void zCopyGTFfeature (const zGTFfeature*, zGTFfeature*);
int  zGTFfeatureCmp  (const zGTFfeature*, const zGTFfeature*);

struct zGTFVec {
	zGTFfeature *elem;   /* array of features */
	int         size;   /* number of features */
	int         limit;  /* number of elements currently allocated */
	zGTFfeature *last;
};
typedef struct zGTFVec zGTFVec;

void zInitGTFVec (zGTFVec *vec, size_t limit);
void zPushGTFVec (zGTFVec *vec, const zGTFfeature *feature);
void zFreeGTFVec (zGTFVec *vec);
zGTFVec* zReadGTFVec (FILE* stream);
void zWriteGTFVec (FILE* stream, const zGTFVec *vec);

void zSortGTFVec(zGTFVec *vec);

struct zGTFConvInfo {
	zStrIdx  gtft2StrIdx[6];
	gtf_t   *strIdx2gtft;
	int      maxStrIdx;
};
typedef struct zGTFConvInfo zGTFConvInfo;

void zInitGTFConvInfo(zGTFConvInfo *conv);
void zFreeGTFConvInfo(zGTFConvInfo *conv);
void zReadGTFConvInfo(FILE* stream, zGTFConvInfo *conv);
void zWriteGTFConvInfo(FILE* stream, const zGTFConvInfo *conv);
int zCDSMatchesStart(zGTFfeature *cds, zGTFfeature *start);
int zCDSMatchesStop(zGTFfeature *cds, zGTFfeature *stop);


#endif /* ZOE_GTF_H */
