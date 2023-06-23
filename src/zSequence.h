/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zSequence.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 

\******************************************************************************/

#ifndef ZOE_SEQUENCE_H
#define ZOE_SEQUENCE_H

#include <stdio.h>
#include <string.h>

#include "zTools.h"

/******************************************************************************\
 zSequence

\******************************************************************************/

struct zSeqVariant{
	coor_t pos;        /* position in actual fasta file */
	coor_t real_pos;   /* position in possibly reversed and paddeded sequence */
	char*  values;     /* array of possible values */
	float* freqs;      /* array of variant frequencies */
	short  variants;   /* number of possible values */
    short  cur;        /* index of current value */
    short  map_idx;    /* index of file_map to which this variant corresponds */
};
typedef struct zSeqVariant zSeqVariant;

/* EVAN haplotypes/snps are DNA specific and shouldn't really be called by 
   those names in the generic zSequence object */
struct zSeqHaplotype{
	int first_snp;     /* index of first var in this haplotype */
	int last_snp;      /* index of last var in this haplotype */
	int snp_count;     /* number of vars in this haplotype */
	short* snp_vals;   /* array of var value indicies */
};
typedef struct zSeqHaplotype zSeqHaplotype;

struct zSeqBlock{
	int      id;       /* unique among blocks in a single zSequence object */
	char    *data;     /* array of sequence data stored by this block */
	coor_t   pos;      /* position in sequence of first index of data */
	coor_t   size;     /* size of data array */
	short    map_idx;  /* index of file_map to which the block corrends */
};
typedef struct zSeqBlock zSeqBlock;

struct zSeqFileMap{
	coor_t      seq_pos;    /* sequence position of this map */
	long int    file_pos;	/* position in input file from which the sequence
							   at seq_pos will be loaded */
	zSeqVariant **vars;     /* list of variants in this section of sequence */ 
    int         var_count;  /* number of variants in this section of sequence */ 
    zSeqBlock*  block;      /* pointer to zSeqBlock containing the corresponding 
							   sequence.  NULL if the sequence is not loaded */
};
typedef struct zSeqFileMap zSeqFileMap;

struct zSequence;
typedef struct zSequence zSequence;

typedef void (*zGenericInitFunc)(void*);
typedef long (*zSequenceInitFunc)(zSequence*, void*);
typedef int (*zSequenceReadFunc)(zSequence*,zSeqBlock*);

struct zSequence {
	coor_t        real_length;
	coor_t        length;
	zList         seq;
	zSeqFileMap  *file_map;
	int           map_size;
	char*         filename;
	FILE*         fp;
	void*         parent;
	coor_t        block_size;
	int           block_count;
	bool          reverse;
	coor_t        padding;
	char          default_char;

	zSeqVariant  **variants;
	int            var_count;

	zSeqHaplotype *haps;
	int            hap_count;

	zSequenceInitFunc    init_func;
	zSequenceReadFunc    read_func;
	zListInitFunc        init_block_func;
};

int  zInitSequence (char*, void*, zSequenceInitFunc, zSequenceReadFunc);
long  zInitSequenceSpecialized (char*, void*, zSequenceInitFunc, 
							   zSequenceReadFunc, zListInitFunc, coor_t, coor_t, long, long);
int zInitMultiSequenceSpecialized (char*, zVec*, zGenericInitFunc ,
			size_t, zSequenceInitFunc, zSequenceReadFunc, zListInitFunc, coor_t, coor_t);
void zFreeSequence (zSequence*);
void zCopySequence (zSequence*, zSequence*);
char zGetSequencePos(zSequence*,coor_t);
void zSetSequencePadding(zSequence*,coor_t);
void zReverseSequence (zSequence*);
int zGetSeqBlockID(zSequence*,coor_t);
void* zCreateSeqBlockWithSize(coor_t);

void zLoadSeqVariants(zSequence*,char*);
void zSetSequenceBlockVariants(zSequence*,zSeqBlock*);
void zSetSequenceVariants(zSequence*);

bool zSeqCheckChar(zSequence*,coor_t,char);
#endif

