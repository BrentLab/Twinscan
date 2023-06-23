/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zProtein.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_PROTEIN_C
#define ZOE_PROTEIN_C

#include "zProtein.h"

void zInitProtein (zProtein *pro, coor_t length, const char *def, const char *seq) {
	pro->length = length;
	pro->def    = zMalloc(strlen(def) +1, "zInitProtein def");
	pro->seq    = zMalloc(pro->length,    "zInitProtein seq");
	(void)strcpy(pro->def, def);
	(void)memcpy(pro->seq, seq, pro->length);
}

void zFastaToProtein (const zFastaFile *fasta, zProtein *pro) {
	(void)zInitProtein(pro, fasta->length, fasta->def, fasta->seq);
}

void zFreeProtein(zProtein *pro) {
	zFree(pro->def);
	zFree(pro->seq);
}


#endif
