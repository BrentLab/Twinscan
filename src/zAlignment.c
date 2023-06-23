/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zAlignment.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_Alignment_C
#define ZOE_Alignment_C

#include <stdio.h>
#include "zAlignment.h"


void zInitAlignment (zAlignment *alignment, coor_t length, const char *def, const char *seq) {
	int i;
	coor_t j;
	int num_seqs = 0;

	/* set definition */
	alignment->def = zMalloc(strlen(def) +1, "zInitAlignment def");
	(void)strcpy(alignment->def, def);

	/* figure out the number of sequences in this alignment */
	strtok ((char *)def, "> ");
	num_seqs = 1;
	while (strtok (NULL, "> ") != NULL) {
		num_seqs++;
	}

	alignment->num_seqs = num_seqs;
	alignment->length = length / num_seqs;

	/* allocate storage for zAlignment */
	alignment->seq = zMalloc(num_seqs * sizeof(char *), "zInitAlignment seq");
	for (i=0; i<num_seqs; i++)
		alignment->seq[i] = zMalloc(alignment->length, "zIinitAlignment seq");

	alignment->s7 = zMalloc(num_seqs * sizeof(char *), "zInitAlignment seq");
	for (i=0; i<num_seqs; i++)
		alignment->s7[i] = zMalloc(alignment->length, "zIinitAlignment seq");
	
	/* set sequences */
	for (i=0; i<num_seqs; i++)
		(void)memcpy(alignment->seq[i], seq + i*alignment->length, alignment->length);
	
	/* create s7 sequences */
	for(i=0; i < num_seqs; i++) {
		for (j=0; j<alignment->length; j++) {
			switch (alignment->seq[i][j]) {
			case 'A':
				alignment->s7[i][j] = 0;
				break;
			case 'C':
				alignment->s7[i][j] = 1;
				break;
			case 'G':
				alignment->s7[i][j] = 2;
				break;
			case 'T':
				alignment->s7[i][j] = 3;
				break;
			case '_':
				alignment->s7[i][j] = 4;
				break;
			case '.':
				alignment->s7[i][j] = 5;
				break;
			case 'N':
				alignment->s7[i][j] = 6;
				break;
			default:
				zDie ("Unrecognized character in alignment: %c\n", 
					  alignment->seq[i][j]);
			}
		}
	}
}

void zFastaToAlignment (const zFastaFile *fasta, zAlignment *alignment) {
	zInitAlignment(alignment, fasta->length, fasta->def, fasta->seq);
}

void zFreeAlignment(zAlignment *alignment) {
	zFree(alignment->def);
	zFree(alignment->seq);
	zFree(alignment->s7);
}

void zCopyAlignment(const zAlignment *alignment, zAlignment *copy) {
	int i;

	/* allocate storage */
	copy->def = zMalloc(strlen(alignment->def) + 1, "zCopyAlignment def");
	copy->seq = zMalloc(alignment->num_seqs * sizeof (char *),"zCopyAlignment seq");
	for (i=0; i<alignment->num_seqs; i++)
		copy->seq[i] = zMalloc(alignment->length, "zCopyAlignment seq");

	copy->s7 = zMalloc(alignment->num_seqs * sizeof (char *), "zCopyAlignment s7");
	for (i=0; i<alignment->num_seqs; i++)
		copy->s7[i] = zMalloc(alignment->length, "zCopyAlignment s7");
	
	/* fill */
	copy->num_seqs = alignment->num_seqs;
	copy->length = alignment->length;

	(void)strcpy(copy->def, alignment->def);
	for (i=0; i<alignment->num_seqs; i++) {
		(void)memcpy(copy->seq[i], alignment->seq[i], alignment->length);
		(void)memcpy(copy->s7[i], alignment->s7[i], alignment->length);
	}
}

void zReverseAlignment (zAlignment *alignment) {

	int i;
	coor_t j;
	char tmp;
	coor_t half_way = alignment->length / 2;

	for (i = 0; i < alignment->num_seqs; i++){
		for (j = 0; j < half_way; j++) {
	
			/* swap seq */
			tmp = alignment->seq[i][j];
			alignment->seq[i][j] = alignment->seq[i][alignment->length - j -1];
			alignment->seq[i][alignment->length -j -1] = tmp;
			
			/* swap s10 */
			tmp = alignment->s7[i][j];
			alignment->s7[i][j] = alignment->s7[i][alignment->length - j -1];
			alignment->s7[i][alignment->length -j -1] = tmp;
		}
	}
}

void zComplementAlignment (zAlignment *alignment) {
	int i;
	coor_t j;

	for (i=0; i<alignment->num_seqs; i++) {
		for (j=0; j<alignment->length; j++) {

			/* complement seq */
			switch (alignment->seq[i][j]) {
			case 'A': alignment->seq[i][j] = 'T'; break;
			case 'C': alignment->seq[i][j] = 'G'; break;
			case 'G': alignment->seq[i][j] = 'C'; break;
			case 'T': alignment->seq[i][j] = 'A'; break;
			}

			/* complement s7 */
			switch (alignment->s7[i][j]) {
			case 0: alignment->s7[i][j] = 3; break;
			case 1: alignment->s7[i][j] = 2; break;
			case 2: alignment->s7[i][j] = 1; break;
			case 3: alignment->s7[i][j] = 0; break;
			}
		}
	}
}

#endif
