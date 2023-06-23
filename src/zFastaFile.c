/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zFastaFile.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_FASTA_FILE_C
#define ZOE_FASTA_FILE_C

#include "zFastaFile.h"

void zFreeFastaFile(zFastaFile *entry) {
	zFree(entry->def);
	zFree(entry->seq);
}

/******************************************************\
  Reads one fasta entry in a fasta file. If there are
  multiple entries, it returns after reading the first
  one without affecting the file stream (it puts the
  '>' character back if it finds another one).

  Returns:
  0 - '>' not the first character, 
      EOF was the first character,
      0 length sequence,
      file stream error.
  1 - found a valid fasta entry and read it 
\******************************************************/

int zReadFastaFile (FILE *stream, zFastaFile* entry){

	/*  added a line for testing */

	char c;
	char defline[4096];
	size_t buffer_size = 1024; /* most sequences are smaller than this */
	size_t buffer_index;       /* used to check if buffer is full */
	char *buffer;           /* the actual buffer, it grows as necessary */

	/* initial check for fasta format */
	c = fgetc(stream);
	if (c == EOF) return 0;
	if (c != '>') {
		zWarn("zReadFastaFile > not found");
		return 0;
	}
	(void)ungetc(c, stream);
	
    /* read the def line */
	(void)fgets(defline, sizeof(defline), stream);
	entry->def = zMalloc(strlen(defline) +1, "zReadFastaFile defline");
	(void)strcpy(entry->def, defline);
	entry->def[strlen(entry->def) -1] = '\0'; /* remove newline */
    
    /* read the sequence */
	buffer_index = 0;
	buffer = zMalloc(buffer_size, "zReadFastaFile buffer"); /* grows */
    
	while ((c = (char) fgetc(stream)) != EOF) {
		if (c == '>') {
			(void)ungetc(c, stream);
			break; /* next record found */
		}
		if (isspace((int)c)) continue; /* skip spaces */
      
		buffer[buffer_index] = c;
		buffer_index++;
      
		if (buffer_index == buffer_size) {
			buffer_size += buffer_size;
			buffer = zRealloc(buffer, buffer_size, "zReadFastaFile buffer");
		}
	}
	
	entry->length = buffer_index;
	entry->seq = buffer;
    
    /* last check for errors */
	if (ferror(stream)) {
		zWarn("zReadFastaFile ferror");
		return 0;
	}
	if (entry->length == 0) {
		zWarn("zReadFastaFile no sequence");
		return 0;
	}
	
    return 1;
}


/******************************************************\
  Reads a multi-fasta file. Since we dont know how 
  many entries there are, it works with a zVec object.
  Calling functions should typecast the vector elements
  into zFastaFile objects.
  Returns: 
     number of entries or -1 if there was an error
\******************************************************/

int zReadMultiFastaFile(FILE* stream, zVec *fasta) {
	zFastaFile *entry = (zFastaFile*) zMalloc(sizeof(zFastaFile), "zReadMultiFastaFile: entry");
	if (fasta == NULL) {
		fasta = (zVec*) zMalloc(sizeof(zVec), "zReadMultiFastaFile: fasta");
		zInitVec(fasta, 2);
	}
	while(zReadFastaFile(stream, entry)) {
		zPushVec(fasta, (void*)entry);
		entry = (zFastaFile*) zMalloc(sizeof(zFastaFile), "zReadMultiFastaFile: entry");
	}
	zFree(entry);
	return fasta->size;
}

static unsigned int zFastaLineLength = 50;

void zSetFastaLineLength (unsigned int length) {
	zFastaLineLength = length;
}

void zWriteFastaFile(FILE *stream, const zFastaFile *entry) {
	coor_t i;

	(void)fprintf(stream, "%s\n", entry->def);
	
	for (i=1; i <= entry->length; i++) {
		(void)fputc((entry->seq)[i-1], stream);
		if ((i % zFastaLineLength) == 0) (void)fprintf(stream, "\n");
	}
	if ((i % zFastaLineLength) != 1) (void)fprintf(stream, "\n");
}

#endif
