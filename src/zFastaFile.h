/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zFastaFile.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_FASTA_FILE_H
#define ZOE_FASTA_FILE_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zTools.h"

/******************************************************************************\
 zFastaFile

zFastaFile reads and writes FASTA format files. Whitespaces will be skipped,
but it does no other sequence checking.

	FILE      *stream;
	zFastaFile fasta;
	
	if ((stream = fopen(filename, "r")) == NULL) error_handler();
	if (!zReadFastaFile(stream, &fasta)) error_handler();
	do_something(&fasta)
	zFreeFastaFile(&fasta);

You can read FASTA databases too.

	while (zReadFastaFile(stream, &fasta) {
		do_something(&fasta);
		zFreeFastaFile(&fasta);
	}
	

\******************************************************************************/

struct zFastaFile {
	coor_t length;
	char   *def;
	char   *seq;
};
typedef struct zFastaFile zFastaFile;

void zFreeFastaFile (zFastaFile*);
int  zReadFastaFile (FILE*, zFastaFile*);
int  zReadMultiFastaFile(FILE* stream, zVec *fasta);
void zSetFastaLineLength (unsigned int);
void zWriteFastaFile (FILE*, const zFastaFile*);

#endif
