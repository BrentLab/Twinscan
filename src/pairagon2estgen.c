/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */  
/*****************************************************************************\
pairagon2estgen.c

Convert a zoe list of features to a GTF file;
\*****************************************************************************/

#include <stdio.h>
#include "ZOE.h"

static const char* usage =
    "Usage: \n\n pairagon2estgen pairagon_file -cdna=cdna_file -genomic=genomic_file\n";

extern int FORWARD, REVERSE, BOTH; /* from zAlnFeature.c */
extern int zChopString(char*, int);
int main(int argc, char** argv) {
	char       *pairagon_file;
	FILE       *pairagon_stream;
	FILE       *stream;
	zAFVec     *afv;
	zDNA       *genomic, **cdna;
	zVec       *cdna_fasta;
	int cdna_entries;
	int i;
	int cdna_orientation;
	int splice_orientation;

	zChar2StrIdx("pairagon2estgen");
	zSetProgramName(argv[0]);
	zParseOptions(&argc, argv);

	if (argc != 2 ||
	    zOption("cdna") == NULL ||
	    zOption("genomic") == NULL)
		zDie(usage);

	pairagon_file = argv[1];

	/* Read in a cDNA fasta file with multiple sequences */
	if ((stream = fopen(zOption("cdna"), "r")) == NULL) {
		zDie("fasta file error (%s)", zOption("cdna"));
	}
	fclose(stream);
	cdna_fasta = (zVec*) zMalloc(sizeof(zVec), "main: cdna_fasta");
	zInitVec(cdna_fasta, 1);
	cdna_entries = zLoadMultiDNAFromMultiFasta(cdna_fasta, zOption("cdna"), NULL);

	cdna = (zDNA**) zMalloc(cdna_entries*sizeof(zDNA*), "main: cdna");
	for (i = 0; i < cdna_entries; i++) {
		cdna[i] = (zDNA*) cdna_fasta->elem[i];
		zSetDNAPadding(cdna[i], PADDING);
	}

	/* Read in genomic fasta file */
	if ((stream = fopen(zOption("genomic"), "r")) == NULL) {
		zDie("fasta file error (%s)", zOption("genomic"));
	}
	fclose(stream);
	genomic = zMalloc(sizeof(zDNA), "main: genomic");
	zInitDNA(genomic);
	zLoadDNAFromFasta(genomic,zOption("genomic"),NULL);
	zSetDNAPadding(genomic, PADDING);

	if ((pairagon_stream = fopen(pairagon_file, "r")) == NULL) {
		zDie("Couldn't open pairagon file %s", pairagon_stream);
	}
	while ((afv = zReadAFVec(pairagon_stream, genomic, NULL, &cdna_orientation, &splice_orientation)) != NULL) {
/* Pass NULL so that afv->elem[n].cdna->def will be set to what's on the output file */
		zDNA *padded_cdna = zMalloc(sizeof(zDNA), "main:cdna");
		char cdna_def[256];
		zInitDNA(padded_cdna);
		for (i = 0; i < cdna_entries; i++) {
			strcpy(cdna_def, cdna[i]->def);
			zChopString(cdna_def, 21); /* Do the exact same thing done to the def when printing the state sequence out and compare */
						   /* Add 1 since there is a '>' in the header, but comparison starts after the '>' - see next line */
			zCopyDNA(cdna[i], padded_cdna);
			if (zChar2StrIdx(cdna_def + 1) == afv->elem[0].cdna_def) {
				if (cdna_orientation == REVERSE) zAntiDNA(padded_cdna);
				break;
			}
		}

		if (padded_cdna == NULL) {
			zDie("Matching cDNA sequence cannot be found for '%s' in %s", zStrIdx2Char(afv->elem[0].cdna_def), zOption("cdna"));
		}
		for (i = 0; i < afv->size; i++) {
			afv->elem[i].cdna = padded_cdna;
			afv->elem[i].genomic = genomic;
		}
			
		fprintf(stdout, "Note Best alignment is between %s est and forward genome, %s splice sites imply %s\n", 
				(cdna_orientation == REVERSE)?"reversed":"forward", 
				(splice_orientation == REVERSE)?"but":"and", 
				(splice_orientation == REVERSE)?"REVERSED GENE":"forward gene");
		if (afv->size == 0) {
			zWarn("Did not read any features");
		} else {
			zWriteAlignment(stdout, afv, 0);
		}
		zFreeDNA(padded_cdna);
	}
	fclose(pairagon_stream);

	zFreeDNA(genomic);
	for (i = 0; i < cdna_entries; i++) {
		zFreeDNA(cdna[i]);
	}
	zFreeVec(cdna_fasta);
	if (afv != NULL) zFreeAFVec(afv);
	zFree(afv);
	zStringPoolFree();
	return 0;
}
