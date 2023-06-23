#include "zHMM.h"

/***************************************************************************
	This is an attempt to avoid hard-coding in the regular
	algorithm. Although the objective was to develop a
	general GPAIRHMM decoder, it becomes necessary to use the
	domain specific knowledge for efficient decoding, and
	sometimes even regular decoding. The functions that
	need the domain specific knowledge don't use that 
	knowledge directly, since that would mean tracking all
	those places everytime that knowledge changes, or every
	time the model changes. They make a functional call to
	get the relevant information, that comes from this file.
	So you are assured that the model/domain changes can be
	done easily by fixing this file. This will only work
	to some extent. For example, you cannot fix this file
	and expect the GPAIRHMM to predict genes in both sequences,
	for example, for that is beyond the scope of this
	generalization. The domain specific knowledge here has to
	do with cDNA-to-genome alignment, splice junctions, U2 and
	U12 intron types, and the mechanism, etc.

	See also: zHardCoding.c
****************************************************************************/

/* Hard coded names of durations and suffix */
#define GENOMIC_NULL_DURATION_NAME ("RandomGenomic")
#define CDNA_NULL_DURATION_NAME    ("RandomCDna")
#define NULL_MODEL_SUFFIX          ("NULL")

/* Functions that test if some state fall into a category*/

int zIsGenomicOnly(zStrIdx state);
int zIsCDnaOnly(zStrIdx state);

int zIsFivePrimeOverhang(zStrIdx state);
int zIsThreePrimeOverhang(zStrIdx state);

int zIsGenomicInsertion(zStrIdx state);
int zIsCDnaInsertion(zStrIdx state);

int zIsIntronEntry(zStrIdx state, strand_t s);
int zIsBranchPoint(zStrIdx state);
int zIsInsideIntron(zStrIdx state);
int zIsIntronExit(zStrIdx state, strand_t s);

int zIsInsideU2Intron(zStrIdx state);
int zIsInsideU12Intron(zStrIdx state);

int zIsU2(zStrIdx state);
int zIsU12(zStrIdx state);

int zIsMatch(zStrIdx state);

/* Derived functions */

int zIsIntron(zStrIdx state);

int zIsExon(zStrIdx state);

int zIsOverhang(zStrIdx state);

int zIsGenomicOverhang(zStrIdx state);
int zIsCDnaOverhang(zStrIdx state);

/* Functions that iterate through the states of the HMM and get the index of the right state */

int zGetU2Donor(zHMM* hmm);
int zGetU2Acceptor(zHMM* hmm);
int zGetU12Donor(zHMM* hmm);
int zGetU12Acceptor(zHMM* hmm);
int zGetFivePrimeCDna(zHMM* hmm);
int zGetFivePrimeGenomic(zHMM* hmm);
int zGetThreePrimeCDna(zHMM* hmm);
int zGetThreePrimeGenomic(zHMM* hmm);
int zGetMatch(zHMM* hmm);
int zGetU12BranchPoint(zHMM* hmm);

int zSetHMMStrand(zHMM* hmm, strand_t strand);
