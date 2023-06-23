/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zProtein.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_PROTEIN_H
#define ZOE_PROTEIN_H

#include <stdlib.h>
#include <string.h>

#include "zFastaFile.h"
#include "zTools.h"

/******************************************************************************\
 zProtein

This is pretty minimal right now since proteins aren't a major focus of this
library. Tere are three ways to create proteins: manually, via zFastaFile, and
via zTranslateDNA.

	zProtein pro;
	zInitProtein(&pro, length,  definition, sequence);
	zFastaToProtein(&fasta, &protein);
	zTranslateDNA(&dna, &pro, frame);
	zFreeProtein(&pro);

The zFastaFile and zPtoein structs are somewhat compatible, so you can cast
zProtein to zFastaFile for printing.

	zWriteFastaFile(stdout, (zFastaFile*)&pro);

\******************************************************************************/


struct zProtein {
	coor_t length;
	char *def;
	char *seq;
};
typedef struct zProtein zProtein;

void zFreeProtein (zProtein*);
void zInitProtein (zProtein*, coor_t, const char*, const char*);
void zFastaToProtein (const zFastaFile*, zProtein*);

#endif
