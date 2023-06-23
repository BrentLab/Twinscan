/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zDNA.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_DNA_H
#define ZOE_DNA_H

#include <stdio.h>
#include <string.h>

#include "zTools.h"
#include "zSequence.h"

/******************************************************************************\
 zDNA

A zDNA contains a length, a definition, and a sequence . The sequence is
simultaneously held in 3 alphabets. zDNA.seq is straight ASCII. zDNA.s5 is a 5
symbol numeric alphabet. zDNA.s16 is a 15 symbol numeric alphabet.

Num  s5    s16
===  ===   ===
 0    A    reserved but not used
 1    C    T  (0001)
 2    G    G  (0010)
 3    T    K  (0011)
 4    N    C  (0100)
 5         Y  (0101)
 6         S  (0110)
 7         B  (0111)
 8         A  (1000)
 9         W  (1001)
10         R  (1010)
11         D  (1011)
12         M  (1100)
13         H  (1101)
14         V  (1110)
15         N  (1111)

There are two ways to create a zDNA object. You can call zInitDNA directly or
you can convert a zFastaFile into a zDNA.

	zDNA dna;
	zInitDNA(&dna1, length, definition, sequence);
	zFastaToDNA(&fasta, &dna)

Several common operations are supported. The reverse, complement, and
anti-parallel operations work on the sequence itself. If you do not want to
mutate your object, copy it first. There are two kinds of subsequence
operations that differ in how they allocate memory. zSelectDNA creates a new
sequence while zSubseqDNA merely moves some pointers around. zSubseq is
therefore much faster, but you shouldn't mutate the parent or child objects.
subsequence definitions are pointers to their parents while selections get a
new string "unnamed selection". zTranslateDNA will convert a zDNA into a
zProtein. You must supply the reading frame (0..2). The translation is done
with the s5 alphabet, so some ambiguous codons with unambiguous translations
will be lost unless they translate from 'N' unambiguously.

	zReverseDNA(&dna);
	zComplementDNA(&dna);
	zAntiDNA(&dna);
	zCopyDNA(&original, &copy);
	zSelectDNA(&dna, &subseq, from, length);
	zSubseqDNA(&dna, &subseq, from, length);
	zTranslateDNA(&dna, &protein, frame);

When you're done with a zDNA object, you should free it. You can call this on
subsequences but it doesn't do anything.

	zFreeDNA(&dna);

The zFastaFile and zDNA structs are somewhat compatible, so you can cast zDNA
to zFastaFile for printing.

	zWriteFastaFile(stdout, (zFastaFile*)&dna);

\******************************************************************************/

struct zDNA {
	coor_t        length;
	char         *def;
	zSequence    *seq;
	float         gc;     /* total gc level */
	float         *gcs;     /* windowed gc level */
	bool          complement;
	zHash *companion; /* Sequences associated with this one, such as a  *
					   * positional GC level, or the unmasked sequence. *
					   * This is loosely typed, and functions that      *
					   * put things in here have to agree ahead of time *
					   * with functions that pull things out.           *
					   */

	int id;
};
typedef struct zDNA zDNA;

void zInitDNA (zDNA*);
void zFreeDNA (zDNA*);
void zCopyDNA (zDNA*, zDNA*);
void zReverseDNA (zDNA*);
void zComplementDNA (zDNA*);
void zAntiDNA (zDNA*);
void zSetDNAPadding(zDNA*,coor_t);

void zLoadDNAFromFasta(zDNA* dna, char* filename,char* snp_filename);
int zLoadMultiDNAFromMultiFasta(zVec *multi_dna, char* filename, char* snp_filename);
char zGetDNAUCSeq(zDNA* dna, coor_t pos);
char zGetDNASeq(zDNA* dna, coor_t pos);
char* zGetDNASeqRange(zDNA* dna, coor_t from, coor_t to);
char zGetDNAS5(zDNA* dna, coor_t pos);
char zGetDNAS16(zDNA* dna, coor_t pos);
float zGetDNAGC(zDNA* dna);
float zGetDNAWindowedGC(zDNA* dna, coor_t pos);

int zSetDNABlockSize(coor_t);
int zSetDNABlockCount(int);

bool zDNACheckA(zDNA* dna, coor_t pos);
bool zDNACheckC(zDNA* dna, coor_t pos);
bool zDNACheckG(zDNA* dna, coor_t pos);
bool zDNACheckT(zDNA* dna, coor_t pos);

void zDNASetSNPs(zDNA* dna);

float* zCalcSlidingGCLevel(zDNA* dna, size_t win);
void zAddDNACompanion(zDNA* dna, const char* key, void* ptr);
void *zGetDNACompanion(zDNA* dna, const char* key);
#endif
