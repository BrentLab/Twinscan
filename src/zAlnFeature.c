/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zAlnFeature.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_AFEATURE_C
#define ZOE_AFEATURE_C

#include "zHardCoding.h"
#include "zAlnFeature.h"

extern zDNA* zMakePaddedDNA (zDNA *real_dna); /* from zPairTrellis.c */
const int FORWARD = 1, REVERSE = 2, BOTH = 3; /* 01, 10, 11 */

void zScore2IntText (score_t val, char *s) {
	     if (val == MIN_SCORE) strcpy(s, ".");
	else if (val == MAX_SCORE) strcpy(s, "*");
	else                       sprintf(s, "%.0f", val);
}

/* zAlnFeature functions */

void zClearAlnFeature (zAlnFeature *f) {
	f->name          = -1; 
	f->state         = -1; 
	f->genomic       = NULL;
	f->genomic_def   = zChar2StrIdx("");
	f->genomic_start = UNDEFINED_COOR; 
	f->genomic_end   = UNDEFINED_COOR;
	f->cdna          = NULL;
	f->cdna_def      = zChar2StrIdx("");
	f->cdna_start    = UNDEFINED_COOR; 
	f->cdna_end      = UNDEFINED_COOR;
	f->strand        = UNDEFINED_STRAND; 
	f->score         = MIN_SCORE; 
	f->padding       = 0;
}

void zFreeAlnFeature (zAlnFeature *f) {
	zClearAlnFeature(f);
}

void zDumpAlnFeature (FILE *stream, const zAlnFeature *f) {
		(void)fprintf(stream, "%-12s\t%12f\t%-3.1f\t%12u\t%12u\t%20s\t%12u\t%12u\t%20s\t%c\t%u\n",
			zStrIdx2Char(f->name), f->score, f->percent_identity, f->genomic_start, f->genomic_end, f->genomic->def, f->cdna_start, f->cdna_end, f->cdna->def, f->strand, f->length);
}

void zWriteAlnFeature (FILE *stream, const zAlnFeature *f, int is_gtf, coor_t g_offset) {
	char genomic_start[16], genomic_end[16], cdna_start[16], cdna_end[16], strand[8], score[64], length[16];
	char *genomic_seqname, *cdna_seqname;
	coor_t base_offset = 1;  /* note: must add 1 to start and end because most people use 1-based coordinates */

	if (f->genomic != NULL)
		genomic_seqname = (char*) zMalloc(strlen(f->genomic->def)*sizeof(char), "zWriteAlnFeature genomic_name");
	else if (!zIsCDnaOnly(f->name))
		genomic_seqname = (char*) zMalloc((strlen(zStrIdx2Char(f->genomic_def)) + 2)*sizeof(char), "zWriteAlnFeature genomic_name");
	else
		genomic_seqname = (char*) zMalloc(2*sizeof(char), "zWriteAlnFeature genomic_name");

	if (f->cdna != NULL)
		cdna_seqname    = (char*) zMalloc(strlen(f->cdna->def)*sizeof(char), "zWriteAlnFeature cdna_name");
	else if (!zIsGenomicOnly(f->name))
		cdna_seqname    = (char*) zMalloc((strlen(zStrIdx2Char(f->cdna_def)) + 2)*sizeof(char), "zWriteAlnFeature cdna_name");
	else
		cdna_seqname    = (char*) zMalloc(2*sizeof(char), "zWriteAlnFeature cdna_name");

	strcpy(genomic_start, " ");
	strcpy(genomic_end, " ");
	strcpy(cdna_start, " ");
	strcpy(cdna_end, " ");
	strcpy(genomic_seqname, " ");
	strcpy(cdna_seqname, " ");

	/* Check zTraceTrellis() to see why this inequality is checked */
	if (f->genomic_start != f->genomic_end) {
		zCoor2Text(f->genomic_start - f->padding + base_offset + 1 + g_offset, genomic_start); /* add 1 to start for closed interval */
		zCoor2Text(f->genomic_end   - f->padding + base_offset + g_offset,     genomic_end);      
		if (f->genomic != NULL)
			strcpy(genomic_seqname, f->genomic->def + 1);
		else
			strcpy(genomic_seqname, zStrIdx2Char(f->genomic_def));
		zChopString(genomic_seqname, 20);
	}
	if (f->cdna_start != f->cdna_end) {
		zCoor2Text(f->cdna_start - f->padding + base_offset + 1, cdna_start);
		zCoor2Text(f->cdna_end   - f->padding + base_offset,     cdna_end);
		if (f->cdna != NULL)
			strcpy(cdna_seqname, f->cdna->def + 1);
		else
			strcpy(cdna_seqname, zStrIdx2Char(f->cdna_def));
		zChopString(cdna_seqname, 20);
	}
	zCoor2Text(f->length, length);
	zStrand2Text(f->strand, strand);
/*
MANI: It looks ugly with %f.
	zScore2Text(f->score, score);
*/
	sprintf(score, "%.0f", f->score);

	if (is_gtf) {
		if (strcmp(zStrIdx2Char(f->name), "Intron") == 0) {
			(void)fprintf(stream, "%s%-12s\t%12s\t%-3.1f\t%12s\t%12s\t%20s\t%12s\t%12s\t%20s\n",
				strand, zStrIdx2Char(f->name), score, f->percent_identity, genomic_start, genomic_end, genomic_seqname, cdna_start, cdna_end, cdna_seqname);
		} else {
			(void)fprintf(stream,   "%-12s\t%12s\t%-3.1f\t%12s\t%12s\t%20s\t%12s\t%12s\t%20s\n",
				        zStrIdx2Char(f->name), score, f->percent_identity, genomic_start, genomic_end, genomic_seqname, cdna_start, cdna_end, cdna_seqname);
		}
	} else {
		(void)fprintf(stream, "%-12s\t%12s\t%-3.1f\t%12s\t%12s\t%20s\t%12s\t%12s\t%20s\t%s\n",
			zStrIdx2Char(f->name), score, f->percent_identity, genomic_start, genomic_end, genomic_seqname, cdna_start, cdna_end, cdna_seqname, strand);
	}
	zFree(genomic_seqname);
	zFree(cdna_seqname);
}

int zReadAlnFeature (FILE *stream, zAlnFeature *f, zDNA *genomic, zDNA *cdna) {
	char line[2048], state_name[17];
	zStrIdx name;
	char genomic_start[16], genomic_end[16], cdna_start[16], cdna_end[16], strand[8], score[64], /*length[16],*/ percent_identity[8];
	char genomic_seqname[32], cdna_seqname[32];
	coor_t base_offset = 1;  /* note: must add 1 to start and end because most people use 1-based coordinates */

	if (fgets(line, sizeof(line)-1, stream) == NULL) return 0;
	if ('#' == line[0]) return 0;
	
	if (sscanf(line, "%s ", state_name) == 1) {
		name = zChar2StrIdx(state_name);
	} else {
		zWarn("zReadAlnFeature format not valid: state");
		return -1;
	}


	cdna_seqname[0] = '\0';
	genomic_seqname[0] = '\0';

	if ((zIsGenomicOnly(name) && sscanf(line, "%12s\t%12s\t%5s\t%12s\t%12s\t%20s %s",
			state_name, score, percent_identity, genomic_start, genomic_end, genomic_seqname, strand) == 7) ||
	    (zIsCDnaOnly(name)    && sscanf(line, "%12s\t%12s\t%s %12s\t%12s\t%s\t%s",
			state_name, score, percent_identity, cdna_start, cdna_end, cdna_seqname, strand) == 7) ||
	    (zIsMatch(name)       && sscanf(line, "%12s %12s %s %12s %12s\t%20s\t%12s %12s %20s %s\n",
			state_name, score, percent_identity, genomic_start, genomic_end, genomic_seqname, cdna_start, cdna_end, cdna_seqname, strand) == 10)) {
	} else {
		zWarn("zReadAlnFeature format not valid for state %s in sequence %s:%s:", zStrIdx2Char(name), genomic->def, line);
		return 0;
	}
	/* consider adding gff someday */
	
	f->padding = PADDING;

	f->genomic_start = f->padding;
	f->genomic_end   = f->padding;
	f->cdna_start    = f->padding;
	f->cdna_end      = f->padding;

	/* set name to enum */
	f->name = zChar2StrIdx(state_name);
	if (zIsGenomicOnly(name) || zIsMatch(name)) {
		f->genomic_start  = zText2Coor(genomic_start) - 1 - base_offset + f->padding;
		f->genomic_end    = zText2Coor(genomic_end) - base_offset + f->padding;
		f->length    = f->genomic_end - f->genomic_start;
	}
	if (zIsCDnaOnly(name) || zIsMatch(name)) {
		f->cdna_start  = zText2Coor(cdna_start) - 1 - base_offset + f->padding;
		f->cdna_end    = zText2Coor(cdna_end) - base_offset + f->padding;
		f->length    = f->cdna_end - f->cdna_start;
	}
	f->strand = zText2Strand(strand);
	f->score  = zText2Score(score);
	f->percent_identity  = zText2Score(percent_identity);
	
	/* subtract 1 from start and end, zoe uses 0-based coordinates internally 
	f->genomic_start--;
	f->genomic_end--;
	f->cdna_start--;
	f->cdna_end--; */

	f->genomic = genomic;
	if (genomic != NULL) {
		f->genomic_def = zChar2StrIdx(f->genomic->def + 1);
	} else {
		if (strlen(genomic_seqname) == 0)
			f->genomic_def = -1;
		else
			f->genomic_def = zChar2StrIdx(genomic_seqname);
	}

	f->cdna = cdna;
	if (cdna != NULL) {
		f->cdna_def = zChar2StrIdx(f->cdna->def + 1);
	} else {
		if (strlen(cdna_seqname) == 0)
			f->cdna_def = -1;
		else
			f->cdna_def = zChar2StrIdx(cdna_seqname);
	}
	f->state   = -1;
	
	if (zVerifyAlnFeature(f)) {
		/* zWriteAlnFeature(stdout, f, 0); */
		return 1;
	} else {
		zWarn("zReadAlnFeature parse error (%s)", line);
		zFreeAlnFeature(f);
		return 0;
	}
}

int zVerifyAlnFeature (const zAlnFeature *f) {
/*
	if (f->genomic_def == NULL) return 0;
	if (f->cdna_def    == NULL) return 0;
*/
	if (f->genomic_start == UNDEFINED_COOR) return 0;
	if (f->genomic_end   == UNDEFINED_COOR) return 0;
	if (f->genomic_start > f->genomic_end)          return 0;
	if (f->cdna_start == UNDEFINED_COOR) return 0;
	if (f->cdna_end   == UNDEFINED_COOR) return 0;
	if (f->cdna_start > f->cdna_end)          return 0;
	if (f->strand != '-' && f->strand != '+' && f->strand != '.')      return 0;
	
	/* passed the tests */
	return 1;
}

void zCopyAlnFeature (const zAlnFeature *orig, zAlnFeature *copy) {
	copy->name          = orig->name;
	copy->state         = orig->state;
	copy->strand        = orig->strand;
	copy->genomic       = orig->genomic;
	copy->genomic_def   = orig->genomic_def;
	copy->genomic_start = orig->genomic_start;
	copy->genomic_end   = orig->genomic_end;
	copy->cdna          = orig->cdna;
	copy->cdna_def      = orig->cdna_def;
	copy->cdna_start    = orig->cdna_start;
	copy->cdna_end      = orig->cdna_end;
	copy->length        = orig->length;
	copy->score         = orig->score;
	copy->percent_identity = orig->percent_identity;
	copy->padding       = orig->padding;
}

int zAntiAlnFeature (zAlnFeature *f) {
	coor_t buffer;
	zStrIdx entry     = zChar2StrIdx("Entry");
	zStrIdx exit      = zChar2StrIdx("Exit");
	zStrIdx entry2    = zChar2StrIdx("Entry2");
	zStrIdx exit2     = zChar2StrIdx("Exit2");
	zStrIdx rgenomic1 = zChar2StrIdx("RGenomic1");
	zStrIdx rcdna1    = zChar2StrIdx("RCDna1");
	zStrIdx rgenomic2 = zChar2StrIdx("RGenomic2");
	zStrIdx rcdna2    = zChar2StrIdx("RCDna2");

	if (f->name == entry) {
		f->name = exit;
	} else if (f->name == exit) {
		f->name = entry;
	} else if (f->name == entry2) {
		f->name = exit2;
	} else if (f->name == exit2) {
		f->name = entry2;
	} else if (f->name == rgenomic1) {
		f->name = rgenomic2;
	} else if (f->name == rgenomic2) {
		f->name = rgenomic1;
	} else if (f->name == rcdna1) {
		f->name = rcdna2;
	} else if (f->name == rcdna2) {
		f->name = rcdna1;
	}

	switch(f->strand) {
		case '+':
			f->strand = '-';
			break;
		case '-':
			f->strand = '+';
			break;
		case '.':
			break;
		default:
			return -1;
			break;
	}

	if (!(f->genomic_start == f->padding && f->genomic_end == f->padding)) {
		buffer = f->genomic_start;
		f->genomic_start = f->genomic->length - 1 - f->genomic_end - 1;
		f->genomic_end   = f->genomic->length - 1 - buffer - 1;
	}

	if (!(f->cdna_start == f->padding && f->cdna_end == f->padding)) {
		buffer = f->cdna_start;
		f->cdna_start    = f->cdna->length - 1 - f->cdna_end - 1;
		f->cdna_end      = f->cdna->length - 1 - buffer - 1;
	}

	return 1;
}

int zAlnFeatureCmp (const zAlnFeature *f1, const zAlnFeature *f2) {
	int    buffer;

	/* Both cdna and genomic are emitted */

	/* genomic_start */
	buffer = zCoortCmp((void*)&f1->genomic_start, (void*)&f2->genomic_start);
	if (buffer != 0) return buffer;
	
	/* genomic_end */
	buffer = zCoortCmp((void*)&f1->genomic_end, (void*)&f2->genomic_end);
	if (buffer != 0) return buffer;
		
	/* cdna_start */
	buffer = zCoortCmp((void*)&f1->cdna_start, (void*)&f2->cdna_start);
	if (buffer != 0) return buffer;
	
	/* cdna_end */
	buffer = zCoortCmp((void*)&f1->cdna_end, (void*)&f2->cdna_end);
	if (buffer != 0) return buffer;
		
	if      (f1->strand < f2->strand) return -1;
	else if (f1->strand > f2->strand) return 1;
	else return f1->name - f2->name;
}

int zAFPtrCmp (const void *v1, const void *v2) {
	return zAlnFeatureCmp( (zAlnFeature*)v1, (zAlnFeature*)v2 );
}

/* zAFVec stuff */

void zInitAFVec (zAFVec *vec, int limit) {
	int i;
	vec->size = 0;
	vec->limit = limit;
	vec->last = NULL;
	if (vec->limit > 0) vec->elem = zMalloc(limit * sizeof(zAlnFeature), "zInitAFVec");
	else                vec->elem = NULL;

	for (i=0; i < limit; i++)
	{
		vec->elem[i].score = MIN_SCORE;
	}    
}

void zPushAFVec (zAFVec *vec, const zAlnFeature *feature) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zRealloc(vec->elem, vec->limit * sizeof(zAlnFeature), "zPushAFVec");
	}
	zCopyAlnFeature(feature, &vec->elem[vec->size]);
	vec->last = &vec->elem[vec->size];
	vec->size++;
}

void zFreeAFVec (zAFVec *vec) {
	if (vec->elem != NULL)
		zFree(vec->elem);
}

/* Reverse the order of elements in the AFVec */
void zFlipAFVec (zAFVec* afv) {
	int i,other;
	zAlnFeature tmp;

	for (i = 0; i < afv->size/2; ++i) {
		other = afv->size - 1 - i;
		zCopyAlnFeature(&afv->elem[i],     &tmp);
		zCopyAlnFeature(&afv->elem[other], &afv->elem[i]);
		zCopyAlnFeature(&tmp,              &afv->elem[other]);
	}
}

void zCopyAFVec (zAFVec* original, zAFVec* copy) {
	int i;

	for (i = 0; i < original->size; ++i) {
		zPushAFVec(copy, &original->elem[i]);
	}
}


/* Dump each feature in the vector */
void zDumpAFVec (FILE* stream, const zAFVec* afv) {
	int i;
	for (i = 0; i < afv->size; ++i) {
		zDumpAlnFeature(stream, &afv->elem[i]);
	}
}

/* Transform the AFVec as if it were on the reverse complement sequences */
void zAntiAFVec (zAFVec* afv) {
	int i;
	for (i = 0; i < afv->size; i++) {
		zAntiAlnFeature(&afv->elem[i]);
	}
	zFlipAFVec(afv);
	qsort(afv->elem, afv->size, sizeof(zAlnFeature), zAFPtrCmp);
}

/* Write each feature in the vector */
void zWriteAFVec (FILE* stream, const zAFVec* afv, int is_gtf, coor_t g_offset) {
	int i;
	for (i = 0; i < afv->size; ++i) {
		zWriteAlnFeature(stream, &afv->elem[i], is_gtf, g_offset);
	}
}

/* Begin Utility functions */

/* Makes the string instring at most length bytes long. char[length-1] is set to '\0' */
int zChopString(char *instring, int length) {
	int i;
	if (instring == NULL || length < 1) {
		return 0;
	}
	for (i = 0; i < (int) strlen(instring); i++) {
		if (isspace((int)instring[i])) {
			instring[i] = '\0';
			break;
		}
	}
	if ((int)strlen(instring) > length - 1) {
		instring[length - 1] = '\0';
	}
	return 1;
}

/* End   Utility functions */

/* Write the alignment in est2genome format */
void zWriteAlignment (FILE* stream, const zAFVec* afv, coor_t g_offset) {
	int i;
	coor_t offset = 1;
	coor_t display_window = 50, max_length = 0, length = 0, current, j, alignment_length = 0;
	score_t score = 0, percent_identity = 0;
	char *genomic_string, *cdna_string, *alignment_string, *genomic_seqname, *cdna_seqname;
	zAFVec gfv;
	coor_t genomic_start = (coor_t) -1;
	coor_t cdna_start = (coor_t) -1;
	int intron_length = 2;

	if (afv->elem[0].genomic != NULL) {
		genomic_seqname = (char*) zMalloc(strlen(afv->elem[0].genomic->def)*sizeof(char), "zWriteAlignment genomic_seqname");
		strcpy(genomic_seqname, afv->elem[0].genomic->def + 1);
	} else {
		genomic_seqname = (char*) zMalloc((strlen(zStrIdx2Char(afv->elem[0].genomic_def)) + 1)*sizeof(char), "zWriteAlignment genomic_seqname");
		strcpy(genomic_seqname, zStrIdx2Char(afv->elem[0].genomic_def));
	}
	if (afv->elem[0].cdna != NULL) {
		cdna_seqname = (char*) zMalloc(strlen(afv->elem[0].cdna->def)*sizeof(char), "zWriteAlignment cdna_seqname");
		strcpy(cdna_seqname, afv->elem[0].cdna->def + 1);
	} else {
		cdna_seqname = (char*) zMalloc((strlen(zStrIdx2Char(afv->elem[0].cdna_def)) + 1)*sizeof(char), "zWriteAlignment cdna_seqname");
		strcpy(cdna_seqname, zStrIdx2Char(afv->elem[0].cdna_def));
	}

	zChopString(genomic_seqname, 20);
	zChopString(cdna_seqname,    20);
	
	/* Begin est_genome header information */

	zInitAFVec(&gfv, 4);
	zAFVec2GFF(afv, &gfv);

	/* No significant alignment */
	if (gfv.size == 0) { 
		return;
	}
	zWriteAFVec(stream, &gfv, 1, g_offset);

	for (i = 0; i < afv->size; ++i) {
		zAlnFeature *af;
		af = &afv->elem[i];
		if (!zIsOverhang(af->name) && !zIsIntron(af->name)) {
			percent_identity += (af->percent_identity*af->length);
			alignment_length += af->length;
		}
		score += af->score;
		max_length += af->length;
	}
	percent_identity /= alignment_length;

	fprintf(stream, "\n%-12s\t%12.0f\t%-4.1f\t%12u\t%12u\t%20s\t%12u\t%12u\t%20s\n\n", 
			"Span", score, percent_identity, gfv.elem[0].genomic_start + 1 + offset + g_offset, gfv.last->genomic_end + offset + g_offset, genomic_seqname, 
			gfv.elem[0].cdna_start + 1 + offset, gfv.last->cdna_end + offset, cdna_seqname);
	fprintf(stream, "\n\n");

	/* End   est_genome header information */

/* Printed the full names, now go back to short names, alright? */

	/* Begin est_genome alignment information, from the -align option */

	zChopString(genomic_seqname, 20);
	zChopString(cdna_seqname,    20);

	fprintf(stream, "%s vs %s:\n\n", genomic_seqname, cdna_seqname);

	genomic_string   = (char*) zMalloc(2*max_length*sizeof(char), "zWriteAlignment genomic");
	cdna_string      = (char*) zMalloc(2*max_length*sizeof(char), "zWriteAlignment genomic");
	alignment_string = (char*) zMalloc(2*max_length*sizeof(char), "zWriteAlignment alignment");
	for (j = 0; j < 2*max_length; j++) {
		genomic_string[j]   = '\0';
		alignment_string[j] = '\0';
		cdna_string[j]      = '\0';
	}
	for (i = 0; i < afv->size; ++i) {
		zAlnFeature *af = &afv->elem[i];
		if (zIsMatch(af->name)) {
			char *genomic = zGetDNASeqRange(af->genomic, af->genomic_start + 1, af->genomic_start + af->length);
			char *cdna    = zGetDNASeqRange(af->cdna,    af->cdna_start + 1, af->cdna_start + af->length);

			for (j = 0; j < strlen(genomic); j++) genomic[j] = toupper(genomic[j]);
			strcat(genomic_string, genomic);

			for (j = 0; j < strlen(cdna); j++) cdna[j] = toupper(cdna[j]);
			strcat(cdna_string, cdna);

			for (j = 0; j < af->length; j++) {
				if (toupper(cdna[j]) == toupper(genomic[j])) {
					strcat(alignment_string, "|");
				} else {
					strcat(alignment_string, " ");
				}
			}
			zFree(genomic);
			zFree(cdna);
			if (genomic_start == (coor_t)-1 && cdna_start == (coor_t)-1) { /* Have not seen a Match yet */
				genomic_start = af->genomic_start - af->padding + 1 + 1; /* add one for 0 based to 1 based coordinates conversion */
				cdna_start    = af->cdna_start    - af->padding + 1 + 1;
			}
			length += af->length;
		} else if (zIsGenomicInsertion(af->name)) { /* WARNING: Only genomic base is emitted */
			char *buffer = zGetDNASeqRange(af->genomic, af->genomic_start + 1, af->genomic_start + af->length);
			for (j = 0; j < strlen(buffer); j++) buffer[j] = toupper(buffer[j]);
			strcat(genomic_string, buffer);
			zFree(buffer);
			for (j = 0; j < af->length; j++) {
				strcat(cdna_string, "-");
				strcat(alignment_string, " ");
			}
			length += af->length;
		} else if (zIsCDnaInsertion(af->name)) { /* WARNING: Only cDNA base is emitted */
			char *buffer = zGetDNASeqRange(af->cdna,    af->cdna_start + 1, af->cdna_start + af->length);
			for (j = 0; j < strlen(buffer); j++) buffer[j] = toupper(buffer[j]);
			strcat(cdna_string, buffer);
			zFree(buffer);
			for (j = 0; j < af->length; j++) {
				strcat(genomic_string, "-");
				strcat(alignment_string, " ");
			}
			length += af->length;
		} else if (zIsIntronEntry(af->name, af->strand) || zIsIntronExit(af->name, af->strand)) {
			coor_t length_of_buffer;
			char *buffer = zGetDNASeqRange(af->genomic, af->genomic_start + 1, af->genomic_start + af->length);
			length_of_buffer = (coor_t) strlen(buffer); /* I didnt want to use it in the loop, since it is mutable in theory */
			for (j = 0; j < length_of_buffer; j++) {
				buffer[j] = tolower(buffer[j]);
			}
			strcat(genomic_string, buffer);
			zFree(buffer);
			for (j = 0; j < af->length; j++) {
				strcat(cdna_string, ".");
				if (af->strand == '+') {
					strcat(alignment_string, ">");
				} else {
					strcat(alignment_string, "<");
				}
			}
			if (zIsIntronEntry(af->name, af->strand)) {
				intron_length = af->length;
			} else {
				intron_length = 0;
			}
			length += af->length;
		} else if (zIsBranchPoint(af->name)) {
			char *buffer;
			coor_t length_of_buffer;
			if (zIsU2(af->name)) {
				continue;
			}
			buffer = zGetDNASeqRange(af->genomic, af->genomic_start + 1, af->genomic_start + af->length);
			length_of_buffer = (coor_t) strlen(buffer); /* I didnt want to use it in the loop, since it is mutable in theory */
			for (j = 0; j < length_of_buffer; j++) {
				buffer[j] = tolower(buffer[j]);
			}
			strcat(genomic_string, buffer);
			zFree(buffer);
			for (j = 0; j < af->length; j++) {
				strcat(cdna_string, ".");
				strcat(alignment_string, "/");
			}
			intron_length = 0;
			length += af->length;
		} else if (zIsInsideIntron(af->name)) {
			/************************************
			 * introns are displayed as follows *
			 *                                  *
			 *    ddiii....iiiaa                *
			 *    >>>>> nn >>>>>                *
			 *    ..............                *
			 *                                  *
			 * d-donor,i-intron,a-acceptor      *
			 * the number of '.'s in the genomic*
			 * sequence depends on the length of*
			 * the intron (I mean the string)   *
			 * so it is 5 arrows, a space, the  *
			 * length, a space, 5 arrows.       *
			 ************************************/
			int short_length;
			int k;
			char *buffer;

			/* Get full length of this intron */
			if (zIsU2(af->name)) {
				if (zIsBranchPoint(afv->elem[i-1].name)) {
					continue;
				} else {
					for (k = i; k < afv->size; k++) {
						intron_length += afv->elem[k].length;
						if (zIsIntronExit(afv->elem[k].name, af->strand)) {
							break;
						}
					}
				}
			} else {
				if (zIsBranchPoint(afv->elem[i+1].name)) {
					intron_length += af->length;
				} else {
					for (k = i; k < afv->size; k++) {
						intron_length += afv->elem[k].length;
						if (zIsIntronExit(afv->elem[k].name, af->strand)) {
							break;
						}
					}
				}
			}

			/**************************************************************************************************************************
			 * OK this is kind of messy. The length of the intron is normally Don+Intron+Acc. When there is a branch point, it        *
			 * becomes Don+Intron, Intron2+Acc. We need to track that since when you are just in Intron you dont know what the length *
			 * is. Can't think of a better way right now. May be because I am on hydrocodone now.                                     *
			 **************************************************************************************************************************/
			buffer  = (char*) zMalloc(256*sizeof(char), "zMalloc buffer");
			if (af->strand == '-') {
				sprintf(buffer, "<<< %u <<<%c", intron_length, '\0'); 
			} else {
				sprintf(buffer, ">>> %u >>>%c", intron_length, '\0');
			}
			short_length = strlen(buffer);
			strcat(alignment_string, buffer);
			for (k = 0; k < short_length; k++) {
				buffer[k] = '.';
			}
			for (k = 0; k < 3; k++) {
				buffer[k]                    = tolower(zGetDNASeq(af->genomic, af->genomic_start + 1 + k));
				buffer[short_length - k - 1] = tolower(zGetDNASeq(af->genomic, af->genomic_end - k));
			}
			strcat(genomic_string, buffer);
			zFree(buffer);
			for (k = 0; k < short_length; k++) {
				strcat(cdna_string, ".");
			}
			length += short_length;
		} else if (zIsThreePrimeOverhang(af->name)) {
			char *buffer  = (char*) zMalloc(2*sizeof(char), "zMalloc buffer");
			buffer[0]='\0';
			strcat(genomic_string, buffer);
			strcat(alignment_string, buffer);
			strcat(cdna_string, buffer);
			zFree(buffer);
		}
	}
	/* End   generating est_genome alignment information, from the -align option */

	/* Begin printing est_genome alignment information, from the -align option */

/* Consider generating strings for each element and finally displaying them together using display_window.*
 * That way, you wont have to sscanf() for the intron length. It's potentially dangerous. 05/21/05        *
 * MANI: No time for that. I guess someone else might have to do it. 06/13/06                             */
	for (current = 0; current < length; current+=display_window) {
		static coor_t intron_length = (coor_t) -1;
		char ac, gc, cc;
		fprintf(stream, "%20s %12u ", genomic_seqname, genomic_start + g_offset);
		for (j = 0; j < display_window; j++) { 
			gc = genomic_string[current+j];
			ac = alignment_string[current+j];
			if (gc == '\0') {
				break; 
			} else {
				fprintf(stream, "%c", gc);
				if (ac == '>' || ac == '<' || ac == '?') {
					if (intron_length == (coor_t)-1) {
						if (sscanf(&alignment_string[current+j+6], "%u", &intron_length) != 0) {
							genomic_start += intron_length;
						}
					}
				} else if (gc != '-' && gc != '.') {
					genomic_start++;
				}
			}

			if (!(ac == '>' || ac == '<' || ac == '?')) {
				intron_length = (coor_t) -1;
			}
		}
		fprintf(stream, " %12u\n", genomic_start - 1 + g_offset);
		fprintf(stream, "%20s %12s ", " ", " ");
		for (j = 0; j < display_window; j++) {
			ac = alignment_string[current+j];
			if (ac == '\0') break; else fprintf(stream, "%c", ac);
		}
		fprintf(stream, " %12s\n", " ");
		fprintf(stream, "%20s %12u ", cdna_seqname, cdna_start);
		for (j = 0; j < display_window; j++) {
			cc = cdna_string[current+j];
			if (cc == '\0') break; else fprintf(stream, "%c", cc);
			if (cc != '-' && cc != '.') {
				cdna_start++;
			}
		}
		fprintf(stream, " %12u\n\n", cdna_start - 1);
	}
	fprintf(stream, "Alignment Score: %f\n", score);

	/* End   printing est_genome alignment information, from the -align option */

	zFree(genomic_string);
	zFree(cdna_string);
	zFree(alignment_string);
	zFree(genomic_seqname);
	zFree(cdna_seqname);
	zFreeAFVec(&gfv);
}

/*******************************************************************
 * Convert af to introns and exons. It keeps a count of the states *
 * in exon_count and intron_count and merges them together to make *
 * a new feature called intron or exon                             *
 *******************************************************************/
void zAFVec2GFF (const zAFVec* afv, zAFVec* gffv) {
	int i, j, intron_flag = 0, exon_flag = 0;
	zAlnFeature gff;
	int intron_begin = -1;
	for (i = 0; i < afv->size; ++i) {
/*			zWriteAlnFeature(stdout, &afv->elem[i]); */
		zAlnFeature *af = &afv->elem[i];
		/* Entering into Intron or 3' overhang, so sign off Exon */
		if (zIsIntronEntry(af->name, af->strand)        ||
		    zIsThreePrimeOverhang(af->name)) {
			score_t percent_identity = 0;
			coor_t  block_length     = 0;
			        intron_flag      = 1;

			/* Intron begins here */
			if (zIsIntronEntry(af->name, af->strand)) {
				intron_begin = i;
			}
			/* Empty alignment */
			if (zIsThreePrimeOverhang(af->name) && exon_flag == 0) {
				break;
			}
			zClearAlnFeature(&gff);
			gff.name          = zChar2StrIdx("Exon");
			gff.state         = 99;
			gff.strand        = afv->elem[i - exon_flag].strand;
			gff.genomic       = af->genomic;
			/* Inside exon features, the only thing that can have genomic coordinates that are skewed is CDna. So account for that */
			if (afv->elem[i - exon_flag].genomic_start != afv->elem[i - exon_flag].genomic_end) {
				gff.genomic_start = afv->elem[i - exon_flag].genomic_start - afv->elem[i - exon_flag].padding;
			} else {
				gff.genomic_start = afv->elem[i - exon_flag + 1].genomic_start - afv->elem[i - exon_flag + 1].padding;
			}
			if (afv->elem[i - 1].genomic_start != afv->elem[i - 1].genomic_end) {
				gff.genomic_end   = afv->elem[i - 1].genomic_end - afv->elem[i - 1].padding;
			} else {
				gff.genomic_end   = afv->elem[i - 2].genomic_end - afv->elem[i - 2].padding;
			}
			gff.cdna          = af->cdna;
			/* Inside exon features, the only thing that can have cdna coordinates that are skewed is Genomic. So account for that */
			if (afv->elem[i - exon_flag].cdna_start != afv->elem[i - exon_flag].cdna_end) {
				gff.cdna_start    = afv->elem[i - exon_flag].cdna_start - afv->elem[i - exon_flag].padding;
			} else {
				gff.cdna_start    = afv->elem[i - exon_flag + 1].cdna_start - afv->elem[i - exon_flag + 1].padding;
			}
			if (afv->elem[i - 1].cdna_start != afv->elem[i - 1].cdna_end) {
				gff.cdna_end      = afv->elem[i - 1].cdna_end - afv->elem[i - 1].padding;
			} else {
				gff.cdna_end      = afv->elem[i - 2].cdna_end - afv->elem[i - 2].padding;
			}
			gff.length        = 0;
			gff.score         = 0;
			for (j = 0; j < exon_flag; j++) {
				gff.length += afv->elem[i - j - 1].length;
				gff.score  += afv->elem[i - j - 1].score;
				percent_identity += afv->elem[i - j - 1].percent_identity*afv->elem[i - j - 1].length;
				block_length += afv->elem[i - j - 1].length;
			}
			gff.percent_identity = percent_identity / block_length;
			zPushAFVec(gffv, &gff);
			exon_flag = 0;
			if (zIsThreePrimeOverhang(af->name)) {
				break;
			}
		}
/* Inside Intron, do pretty much nothing */
		if (zIsInsideIntron(af->name)) {
			intron_flag = 1;
			exon_flag = 0;
		} 
/* Exiting Intron, so sign off Intron */
		if (zIsIntronExit(af->name, af->strand)) {
			char intron_name[20];
			intron_flag = 1;
			exon_flag = 0;
			zClearAlnFeature(&gff);
			gff.state         = 99;
			gff.strand        = af->strand;
			zStrand2Text(gff.strand, intron_name);
			strcat(intron_name, "Intron");
			gff.name          = zChar2StrIdx(intron_name);
			gff.genomic       = af->genomic;
			gff.genomic_start = afv->elem[intron_begin].genomic_start - afv->elem[intron_begin].padding;
			gff.genomic_end   = af->genomic_end - af->padding;
			gff.cdna          = af->cdna;
			if (afv->elem[intron_begin].cdna_start != afv->elem[intron_begin].cdna_end) {
				gff.cdna_start    = afv->elem[intron_begin].cdna_start - afv->elem[intron_begin].padding;
			} else {
				gff.cdna_start    = afv->elem[intron_begin+1].cdna_start - afv->elem[intron_begin+1].padding;
			}
			if (af->cdna_start != af->cdna_end) {
				gff.cdna_end      = af->cdna_end - af->padding;
			} else {
				gff.cdna_end      = afv->elem[i-1].cdna_end - afv->elem[i-1].padding;
			}
			gff.length        = 0;
			gff.score         = 0;
			/* Calculate intron length */
			for (j = intron_begin; j <= i; j++) {
				gff.length += afv->elem[j].length; 
				gff.score  += afv->elem[j].score;
			}
			gff.percent_identity = 0.0;
			zPushAFVec(gffv, &gff);
			intron_flag = 0;
			exon_flag = 0;
/* Entering or continuing Exon, so keep count */
		} else if (zIsExon(af->name)) {
			intron_flag = 0;
			exon_flag++;
		}
	}

/* If the last portion was exonic, include that */
	if (exon_flag > 0) {
			score_t percent_identity = 0;
			coor_t block_length = 0;
			i = afv->size;
			zClearAlnFeature(&gff);
			gff.name          = zChar2StrIdx("Exon");
			gff.state         = 99;
			gff.strand        = afv->elem[i - exon_flag].strand;
			gff.genomic       = afv->elem[i - 1].genomic;
			/* Inside exon features, the only thing that can have genomic coordinates that are skewed is CDna. So account for that */
			if (afv->elem[i - exon_flag].genomic_start != afv->elem[i - exon_flag].genomic_end) {
				gff.genomic_start = afv->elem[i - exon_flag].genomic_start - afv->elem[i - exon_flag].padding;
			} else {
				gff.genomic_start = afv->elem[i - exon_flag + 1].genomic_start - afv->elem[i - exon_flag + 1].padding;
			}
			if (afv->elem[i - 1].genomic_start != afv->elem[i - 1].genomic_end) {
				gff.genomic_end   = afv->elem[i - 1].genomic_end - afv->elem[i - 1].padding;
			} else {
				gff.genomic_end   = afv->elem[i - 2].genomic_end - afv->elem[i - 2].padding;
			}
			gff.cdna          = afv->elem[i - 1].cdna;
			/* Inside exon features, the only thing that can have cdna coordinates that are skewed is Genomic. So account for that */
			if (afv->elem[i - exon_flag].cdna_start != afv->elem[i - exon_flag].cdna_end) {
				gff.cdna_start    = afv->elem[i - exon_flag].cdna_start - afv->elem[i - exon_flag].padding;
			} else {
				gff.cdna_start    = afv->elem[i - exon_flag + 1].cdna_start - afv->elem[i - exon_flag + 1].padding;
			}
			if (afv->elem[i - 1].cdna_start != afv->elem[i - 1].cdna_end) {
				gff.cdna_end      = afv->elem[i - 1].cdna_end - afv->elem[i - 1].padding;
			} else {
				gff.cdna_end      = afv->elem[i - 2].cdna_end - afv->elem[i - 2].padding;
			}
			gff.length        = 0;
			gff.score         = 0;
			for (j = 0; j < exon_flag; j++) {
				gff.length += afv->elem[i - j - 1].length;
				gff.score  += afv->elem[i - j - 1].score;
				percent_identity += afv->elem[i - j - 1].percent_identity*afv->elem[i - j - 1].length;
				block_length += afv->elem[i - j - 1].length;
			}
			gff.percent_identity = percent_identity / block_length;
			zPushAFVec(gffv, &gff);
	}
}

zAFVec* zReadAFVec (FILE* stream, zDNA *genomic, zDNA *cdna, int *cdna_orientation, int *splice_orientation) {
	zAFVec      *afv;
	zAlnFeature  f;
	char         line[2048];
	char         orientation[32];
	char         splice[32];
	char         conjunction[32];
	score_t      fscore;
	int          i, j;

	afv = zMalloc(sizeof(zAFVec), "zReadAFVec");
	zInitAFVec(afv, 1);
	do {
		if (fgets(line, sizeof(line)-1, stream) == NULL) {
			zFreeAFVec(afv);
			zFree(afv);
			afv = NULL;
			return NULL;
		}
		if (sscanf(line, "# Note Best alignment is between %s est and forward genome, %s splice sites imply %s", orientation, conjunction, splice) == 3) {
			if (strcmp(orientation, "forward") == 0) {
				*cdna_orientation = FORWARD;
			} else if (strcmp(orientation, "reversed") == 0) {
				*cdna_orientation = REVERSE;
			} else {
				*cdna_orientation = BOTH;
			}
			if (strcmp(splice, "forward") == 0) {
				*splice_orientation = FORWARD;
			} else if (strcmp(splice, "REVERSED") == 0 || strcmp(splice, "reversed") == 0) {
				*splice_orientation = REVERSE;
			} else {
				*splice_orientation = BOTH;
			}
		}
		if (sscanf(line, "# Optimal score: %lf", &fscore) == 1) break;
	} while ('#' == line[0]);

	while (zReadAlnFeature(stream, &f, genomic, cdna)) {
		if (f.cdna_start == f.padding && f.cdna_end == f.padding) {
			if (afv->last != NULL) {
				f.cdna_start = afv->last->cdna_end;
				f.cdna_end   = afv->last->cdna_end;
			}
		}
		if (f.genomic_start == f.padding && f.genomic_end == f.padding) {
			if (afv->last != NULL) {
				f.genomic_start = afv->last->genomic_end;
				f.genomic_end   = afv->last->genomic_end;
			}
		}
		zPushAFVec(afv, &f);
	}
/***************************************************************************************************
 * Before this point, accessing afv->elem[0].cdna_def might throw error. Accessing cdna_def should *
 * be done with EXTREME CARE, since you might throw SERIOUS errors. But after this, it should      *
 * not! Every element should have the right cdna_def set up. Same for genomic_def too!             *
 ***************************************************************************************************/
	for (i = 0; i < afv->size; i++) {
		if (afv->elem[i].cdna_def != -1 && afv->elem[i].cdna_def != zChar2StrIdx("")) {
			for (j = 0; j < afv->size; j++) {
				if (j != i) {
					afv->elem[j].cdna_def = afv->elem[i].cdna_def;
				}
			}
			break;
		}
	}
	return afv;
}

void zTranslateAFVec (zAFVec* afv, int dist) {
	int i;
	for (i=0; i < afv->size; ++i) {
		afv->elem[i].genomic_start += dist;
		afv->elem[i].genomic_end   += dist;
		afv->elem[i].cdna_start += dist;
		afv->elem[i].cdna_end   += dist;
	}
}


/**********************************************************************\
  zAFList - linked list of zAlnFeature objects
\**********************************************************************/

static unsigned int zAFListCount = 0; 
static zAFListNode* zAFListDeadNodes;

void zIncrementAFListPool(){
	if(zAFListCount == 0){
		zAFListDeadNodes = NULL;	
	}
	zAFListCount++;
}

void zDecrementAFListPool(){
	zAFListNode* n;
    zAFListCount--;
	if(zAFListCount == 0){
		while(zAFListDeadNodes != NULL){
			n = zAFListDeadNodes;
			zAFListDeadNodes = zAFListDeadNodes->next;    
			zFree(n);
		}
	}
}

zAFListNode* zAFListAcquireNode(){
	zAFListNode* temp;
	if(zAFListDeadNodes == NULL){
		temp = zMalloc(sizeof(zAFListNode),"zAFListAcquireNode head");
	}
	else{
		temp = zAFListDeadNodes;
		zAFListDeadNodes = zAFListDeadNodes->next;
	}		
	zClearAlnFeature(&temp->data);
	temp->next = NULL;
	temp->prev = NULL;
	return temp;
}

void zAFListReleaseNode(zAFListNode* n){
	n->next = zAFListDeadNodes;
	n->prev = NULL;
	zAFListDeadNodes = n;
}

void zInitAFList(zAFList* l){
	zIncrementAFListPool();
	l->head = zAFListAcquireNode();
	l->tail = zAFListAcquireNode();
	l->current = l->head;
	l->head->prev = NULL;
	l->head->next = l->tail;
	l->tail->prev = l->head;
	l->tail->next = NULL;
	l->size = 0;
}

void zFreeAFList(zAFList* l){
	l->current = l->head->next;
	while(l->current != l->tail){
		l->current = l->current->next;
		zAFListReleaseNode(l->current->prev);
	}
	zAFListReleaseNode(l->head);
	zAFListReleaseNode(l->tail);
	l->size = 0;
    zDecrementAFListPool();
}

void zResetAFList(zAFList* l){
	l->current = l->head->next;
	while(l->current != l->tail){
		l->current = l->current->next;
		zAFListReleaseNode(l->current->prev);
	}
	l->head->next = l->tail;
	l->tail->prev = l->head;
    l->current = l->head;
	l->size = 0;
}

zAlnFeature* zAFListMoveFirst(zAFList* l){
	l->current = l->head->next;
	if(l->current == l->tail){
		return NULL;
	}
	return &l->current->data;
}

zAlnFeature* zAFListMoveLast(zAFList* l){
	l->current = l->tail->prev;
	if(l->current == l->head){
		return NULL;
	}
	return &l->current->data;
}

zAlnFeature* zAFListMoveNext(zAFList* l){
    if(l->current != l->tail){
		l->current = l->current->next;
		if(l->current != l->tail){
			return &l->current->data;
		}
	}
	return NULL;
}

zAlnFeature* zAFListMovePrev(zAFList* l){
    if(l->current != l->head){
		l->current = l->current->prev;
		if(l->current != l->head){
			return &l->current->data;
		}
	}
	return NULL;
}

zAlnFeature* zAFListGetCurrent(zAFList* l){
	if((l->current == l->head) ||
	   (l->current == l->tail)){
		return NULL;
	}
	else{
		return &l->current->data;
	}
}

zAlnFeature* zAFListAppend(zAFList* l){
	zAFListNode* n;
	n = zAFListAcquireNode();
	n->next = l->tail;
	n->prev = l->tail->prev;
	l->tail->prev->next = n;
	l->tail->prev = n;
	l->size++;
	return &n->data;
}

zAlnFeature* zAFListPrepend(zAFList* l){
	zAFListNode* n;
	n = zAFListAcquireNode();
	n->next = l->head->next;
	n->prev = l->head;
	l->head->next->prev = n;
	l->head->next = n;
	l->size++;
	return &n->data;
}

zAlnFeature* zAFListInsert(zAFList* l,zAlnFeature* f){
	zAFListNode* n;
	n = zAFListAcquireNode();
	zCopyAlnFeature(f,&n->data);
	n->next = l->tail;
	n->prev = l->tail->prev;
	l->tail->prev->next = n;
	l->tail->prev = n;
	l->size++;
	return &n->data;
}

zAlnFeature* zAFListInsertNext(zAFList* l,zAlnFeature* f){
	zAFListNode* n;
	n = zAFListAcquireNode();
	zCopyAlnFeature(f,&n->data);
	if(l->current == l->tail){
		return zAFListInsertPrev(l,f);
	}
	n->next = l->current->next;
	n->prev = l->current;
	l->current->next->prev = n;
	l->current->next = n;
	l->size++;
	return &n->data;
}

zAlnFeature* zAFListInsertPrev(zAFList* l,zAlnFeature* f){
	zAFListNode* n;
	n = zAFListAcquireNode();
	zCopyAlnFeature(f,&n->data);
	if(l->current == l->head){
		return zAFListInsertNext(l,f);
	}
	n->next = l->current;
	n->prev = l->current->prev;
	l->current->prev->next = n;
	l->current->prev = n;
	l->size++;
	return &n->data;
}

void zAFListRemoveLast(zAFList* l){
    zAFListNode* n = l->tail->prev;
	if(n != l->head){
		n->prev->next = l->tail;
		l->tail->prev = n->prev;
		zAFListReleaseNode(n);
	}
	l->size--;
}

void zAFListRemoveFirst(zAFList* l){
    zAFListNode* n = l->head->next;
	if(n != l->tail){
		n->next->prev = l->head;
		l->head->next = n->next;
		zAFListReleaseNode(n);
	}
	l->size--;
}

void zAFListRemoveCurrent(zAFList* l){
	zAFListNode* n = l->current;
	l->current = n->next;
	if((n == l->head) || (n == l->tail)){ 
		return;
	}
	n->prev->next = n->next;
	n->next->prev = n->prev;
	zAFListReleaseNode(n);
	l->size--;
}

void zAFList2AFVec(zAFList* l,zAFVec* v){
	zAFListNode* n = l->head->next;
	while(n != l->tail){
		zPushAFVec(v,&n->data);
		n = n->next;
	}
}

#endif
