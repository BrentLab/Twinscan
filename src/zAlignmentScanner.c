/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zScanner.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_ALIGNMENTSCANNER_C
#define ZOE_ALIGNMENTSCANNER_C

#include <assert.h>
#include "zAlignmentScanner.h"

score_t zGetUScoreA(const zAlignmentScanner *scanner,coor_t i){

	if ((i < scanner->min_pos) || (i > scanner->max_pos)) {
		return MIN_SCORE;
	}

	if(scanner->uscore != NULL){
		return scanner->uscore[i];
	}
	else{
		return scanner->score(scanner,i);
	}
}

score_t zGetUScoreARange(const zAlignmentScanner *scanner, coor_t start, coor_t end){
	score_t score,s;
	coor_t i;

	if ((start < scanner->min_pos) || (end > scanner->max_pos)){
		return MIN_SCORE;
	}
	if(end < start) {
		return 0;
	}
	
	if(scanner->uscore2 != NULL){
		return scanner->uscore2[end] - scanner->uscore2[start-1];
	}

	score = 0;
	for(i = start; i <= end; i++){
		s = zGetUScoreA(scanner,i);
		if(s == MIN_SCORE) continue;
		score += s;
	}
	return score;
}

/* converts an int to a string of its bytes */
/*
static char* int_to_string (int x) {
	int i;
	int remainder;
	char *string;

	string = zMalloc(sizeof(int)+1, "int_to_string");
	remainder = x;
	for (i=sizeof(int)-1; i>=0; i--) {
		string[i] = (char)(remainder / intpow(2, 8*i));
		remainder -= intpow(2, 8*i);
	}
	assert (remainder == 0);
	string[sizeof(int)] = ' ';

	return string;
}
*/

/* Runs Bayesian inference to score a tuple at the given position.
   Models remember if they have scored a given tuple before
   and record the score to save work in the future */
static score_t zBNTreeInfer (const zBNTreeModel *model, zDNA *dna,
							 const zAlignment *alignment, coor_t pos) {
	int i, j;
	int tuple_idx;
	char tuple[256];
	double prob;
	score_t score;

	/* figure out the tuple at this position */
	tuple_idx = 0;
	for (i = -1 * model->order; i <= 0; i++) {
		if (zGetDNAS5(dna,pos + i) == 4) return -10;  /* found N */

		tuple[tuple_idx++] = zGetDNAUCSeq(dna,pos + i);
		for (j=0; j<alignment->num_seqs; j++) {
			if (alignment->seq[j][pos + i] == 'N') return -10; /* found N */
			tuple[tuple_idx++] = alignment->seq[j][pos + i];
		}
		if (i != 0)
			tuple[tuple_idx++] = ' ';
	}
	tuple[tuple_idx] = '\0';

	prob = calc_tuple_str_prob (model->bntree, tuple, 1);
	score = 10 * zLog2(prob);

	return score;
}

static score_t zScoreBNTREE (const zAlignmentScanner *scanner, coor_t pos) {

	if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) {
		return MIN_SCORE;
	}

	if(scanner->uscore != NULL){
		return scanner->uscore[pos];
	}

	return zBNTreeInfer (scanner->model, scanner->dna, scanner->alignment, pos);
}

static score_t zScoreBNTREE_ARRAY (const zAlignmentScanner *scanner, coor_t pos) {
	int i;
	coor_t mfocus;
	double score = 0;
	zBNTreeModel *this_submodel;

	if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) {
		return MIN_SCORE;
	}

	if(scanner->uscore != NULL){
		return scanner->uscore[pos];
	}
	
	mfocus = pos - scanner->model->focus;
	for (i=0; i<scanner->model->submodels; i++) {
		this_submodel = &(scanner->model->submodel[i]);
		score += zBNTreeInfer (this_submodel, scanner->dna, scanner->alignment, mfocus + i);
	}

	return score;
}

static score_t zScoreBNTREE_CDS (const zAlignmentScanner *scanner, zSfeature *f) {
	coor_t i;
	score_t score;
	coor_t start, end;
	int frame = -1;
	
	if (f->name == scanner->Einit || f->name == scanner->Esngl) {
		start = f->start +6;
		end   = f->end - 3;
	}
	else if (f->name == scanner->Exon || f->name == scanner->Eterm) {
		start = f->start + 3;
		end = f->end - 3;
	}
	else {
		start = end = 0;
		zDie("zScoreBNTREE_CDS: attempt to score non-CDS feature %s with BNTREE_CDS model", f->name);
	}
	
	/* figure out subscanner index to start with */
	switch (f->lfrag) {
	  case 0: frame = 0; break;
	  case 1: frame = 2; break;
	  case 2: frame = 1; break;
	  default: zDie("zScoreBNTREE_CDS:  invalid exon lfrag");
	}
	
	score = 0;

	if (scanner->subscanner[0].uscore != NULL) {
		for (i = start; i <= end; i++) {
			score += scanner->subscanner[frame].uscore[i];
			frame = (frame + 1) % 3;
		}
	}
	else {
		for (i = start; i <= end; i++) {
			score += scanner->subscanner[frame].score (&scanner->subscanner[frame], i);
			frame = (frame + 1) % 3;
		}
	}
	
	return score;
}

/*********\
   Range
\*********/
static score_t zAlignmentScoreFeature (const zAlignmentScanner *scanner, zSfeature *f) {
	score_t score;
	coor_t i;

	if (f->lfrag || f->rfrag || f->frame) {
		zWriteSfeature(stderr, f);
		zDie("zAlignmentScoreFeature expects lfrag, rfrag, f->frag to be zero");
	}

	if ((f->start < scanner->min_pos) || (f->end > scanner->max_pos)){
		return MIN_SCORE;
	}

	if(f->end < f->start){
		return 0;
	}

	if(scanner->uscore2 != NULL){
		return scanner->uscore2[f->end] - scanner->uscore2[f->start-1];
	}

	score = 0;
	for(i = f->start; i <= f->end; i++){
		score += zGetUScoreA(scanner,i);
	}
	return score;
}

/*********************\
   ILLEGAL FUNCTIONS - prevents the use of score and count for CDS
\*********************/
static score_t zIllegalAlignmentScore (const zAlignmentScanner *s, coor_t p) {
	zDie("illegal score in AlignmentScanner (%s), at %d", s->model->name, p);
	return 0;
}

/********\
   Init
\********/
static int ID = 0;
void zInitAlignmentScanner(zAlignmentScanner *scanner, zDNA *dna, zAlignment *alignment, zBNTreeModel *model) {
	int i;

	/* set DNA, alignment, and model, clear pointers */	
	scanner->dna       = dna;
	scanner->alignment = alignment;
	scanner->model     = model;
	scanner->uscore    = NULL;
	scanner->uscore2   = NULL;
	scanner->score     = NULL;
	scanner->scoref    = NULL;
	scanner->id        = ID++;

	/* determine scoring region */
	scanner->min_pos = model->length + model->order - 1;
	scanner->max_pos = alignment->length - model->length - model->order - 2;

	/* bind scoring and counting functions to type of model */
	switch (model->type) {
	  case BNTREE:       scanner->score = zScoreBNTREE; break;
	  case BNTREE_ARRAY: scanner->score = zScoreBNTREE_ARRAY; break;
	  default: scanner->score = zIllegalAlignmentScore;
	}

	/* bind range scoring functions */
	switch (model->type) {
	  case BNTREE_CDS:
		scanner->Einit  = zChar2StrIdx("Einit");
		scanner->Exon   = zChar2StrIdx("Exon");
		scanner->Esngl  = zChar2StrIdx("Esngl");
		scanner->Eterm  = zChar2StrIdx("Eterm");
		scanner->scoref = zScoreBNTREE_CDS; 
		break;
	  default:  scanner->scoref = zAlignmentScoreFeature;
	}

	/* set up subscanners */
	switch (model->type) {
	  case BNTREE_CDS:
		scanner->subscanners = 3;
		scanner->subscanner = zMalloc (3 * sizeof(zAlignmentScanner), "zInitAlignmentScanner subscanners");
		for (i=0; i<3; i++)
			zInitAlignmentScanner (&scanner->subscanner[i], scanner->dna,
								   scanner->alignment, &scanner->model->submodel[i]);
		break;
	  default:
		scanner->subscanners = 0;
		scanner->subscanner = NULL;
	}
}

void zFreeAlignmentScanner(zAlignmentScanner *scanner) {
	zFree(scanner->uscore);     
	scanner->uscore     = NULL;
	scanner->model  = NULL;
	scanner->dna = NULL;
	scanner->alignment = NULL;
	scanner->score  = NULL;
	scanner->scoref = NULL;
}

void zPreComputeAlignmentScanner(zAlignmentScanner* scanner) {
/* Argument use_run_sum controls whether the scanner is pre-computed as a  */
/* "running sum" (uscore2[i] = uscore2[i-1] + scanner->score(scanner, i)) or */
/* as a score array (uscore2[i] = scanner->score(scanner, i))               */

  coor_t i;
  int frame;
  score_t  *u,*u2;

  if (NULL != scanner->uscore)  {
	return;
  }
  
  if (scanner->model->type == BNTREE_CDS) {
	  /* precompute all subscanners */
	  for (frame=0; frame<3; frame++) {
		  zPreComputeAlignmentScanner(&scanner->subscanner[frame]);
	  }
  }
  else {
	  u = zMalloc(scanner->alignment->length * sizeof(score_t),
								"zPreComputeAlignmentScanner uscore");
	  u2 = zMalloc(scanner->alignment->length * sizeof(score_t),
								 "zPreComputeAlignmentScanner uscore2");
	  u[0] = scanner->score(scanner, 0);
	  if(u[0] == MIN_SCORE){
		  u2[0] = 0;
	  }
	  else{
		  u2[0] = u[0];
	  }
	  for (i = 1; i < scanner->alignment->length; i++){
		  u[i] = scanner->score(scanner, i);
		  if (u[i] == MIN_SCORE){
			  u2[i] = u2[i-1];
		  }
		  else{
			  u2[i] = u[i] + u2[i-1];
		  }
	  }
	  scanner->uscore = u;
	  scanner->uscore2 = u2;
  }
}

#endif
