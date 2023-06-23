/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */  
/******************************************************************************\
zAlignment.h - part of the ZOE library for genomic analysis

\******************************************************************************/

#ifndef ZOE_Alignment_H
#define ZOE_Alignment_H

#include <stdio.h>
#include <string.h>

#include "zFastaFile.h"
#include "zProtein.h"
#include "zSfeature.h"
#include "zTools.h"

/******************************************************************************\
 zAlignment

zAlignment is similar to zConseq, but holds a multiple sequence alignment
of informants to the target.

\******************************************************************************/


struct zAlignment {
    int num_seqs;  /* number of sequences in this alignment */
	coor_t length; /* length of each sequence */
	char  *def;    /* header */
	char  **seq;    /* raw sequences */
	char  **s7;     /* sequences converted to numbers */
};
typedef struct zAlignment zAlignment;

void zFreeAlignment (zAlignment*);
void zInitAlignment (zAlignment*, coor_t, const char*, const char*);
void zFastaToAlignment (const zFastaFile*, zAlignment*);
void zCopyAlignment (const zAlignment*, zAlignment*);
void zReverseAlignment (zAlignment*);
void zComplementAlignment (zAlignment*);

#endif
