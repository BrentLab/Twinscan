/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zFeatureFactory.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_FEATUREFACTORY_C
#define ZOE_FEATUREFACTORY_C

/* Definitions related to twinscan compatibility */
#define MINIMUM_TWINSCAN_SIGNAL_SCORE -100 /* -100 may be a litte bit high */
#define MIN_POLYA_SCORE -105
#define SIGNAL_PEPTIDE_LENGTH 20

/* UTR maximum feature lengths */
#define MAX_EPA_LENGTH 400
#define MAX_EP_LENGTH 400
#define MAX_ENC_LENGTH 400
#define MAX_EA_LENGTH 200

/* It has been observed that codon usage statistics in first 150 bases  */
/* of Esngl and Einit differ from the rest of the exonic sequences      */
/* Thus, first 150 bases of Einits and Esngls are scored with EINIT_CDS */
/* and the remainder (if any) with NONEINIT_CDS                         */
#define EINIT_CDS_LENGTH 150

#include "zFeatureFactory.h"

/* Definitions used for the exon class of factories */
#define EINIT_CDS     0
#define NONEINIT_CDS  1
#define START         2
#define STOP          3
#define ACC           4
#define DON           5
#define SIG_PEP       6
#define EP            7 /* content model for Ep, Epa features */
#define EA            8 /* content model for Ea, Enc features */
#define NUM_EXON_SCANNERS 9

#define CDSCONS    0
#define INTRONCONS 1
#define INITCONS   2
#define TERMCONS   3
#define ACCCONS    4
#define DONCONS    5
#define UTRCONS    6
#define UTRDONCONS 7
#define UTRACCCONS 8
#define NUM_EXON_CONSCANNERS 9

#define CDSEST    0
#define INTERGENICEST 1
#define INITEST   2
#define TERMEST   3
#define ACCEST    4
#define DONEST    5
#define EPEST     6
#define ENCEST    7
#define EAEST     8
#define EPAEST    9
#define UTRACCEST 10
#define UTRDONEST 11
#define NUM_EXON_ESTSCANNERS 12

#define CDSALIGN 0
#define INTRONALIGN 1
#define INITALIGN 2
#define TERMALIGN 3
#define ACCALIGN 4
#define DONALIGN 5
#define EPA_ALIGN 6
#define EP_ALIGN 7
#define ENC_ALIGN 8
#define EA_ALIGN 9
#define UTRACCALIGN 10
#define UTRDONALIGN 11
#define NUM_EXON_ALIGNSCANNERS 12

#define TATA              0
#define PB_CAP            1
#define NUM_PROM_SCANNERS 2

static score_t zScoreRange (zScanner* scanner, coor_t start, coor_t end){
	return scanner->scorer(scanner,start,end);
}
score_t tempfunc(zScanner* scanner, coor_t start, coor_t end){
	coor_t i;
	score_t s,score = 0;
	for(i = start; i <= end; i++){
		s = scanner->score(scanner,i);
		if(s == MIN_SCORE) continue; /* This should probably return MIN_SCORE instead
										of continuing but this is how it was before and
										I am not certain that it should change
										emk 9-26-05 */
		score += s;
	}
	return score;
}

static score_t zScoreRangeA (const zAlignmentScanner* scanner, coor_t start, coor_t end){
	/*	coor_t i;
	score_t s,score = 0;
	for(i = start; i <= end; i++){
		s = zGetUScoreA(scanner,i);
		if(s == MIN_SCORE) continue; 
		score += s;
	}
	s = zGetUScoreARange(scanner,start,end);
	if(abs(s - score) > .000001){
		zDie("wtf\n");
	}
	return score;*/
	return zGetUScoreARange(scanner,start,end);
}

static score_t zScoreCoding (const zFeatureFactory* fac, zSfeature *e) {
	int     frame = 0;
	int     length;
	score_t coding_score;

	/* convert exon frame to cds[frame] -  I know, it looks weird */
	/* plus, the funky precomputation encoding scheme is undocumented */
	switch (e->frame) {
	case 0: frame = 2; break;
	case 1: frame = 1; break;
	case 2: frame = 0; break;
	default: zDie("zScoreCoding impossible");
	}
	
	/* -3 in uscore[e->end -3] added for genscan compatibiliity */
	coding_score = 0;
	
	/* added this as part of the let-very-short-exons-live campaing */
	length = e->end - e->start + 1;
	
	if ( (((e->end - 3) - e->start + 1) > fac->offset) && (length > 2) ) {
		coding_score = zGetRangeScore(fac->scan[NONEINIT_CDS],e->end-3,e->start + fac->offset,frame); 
	}
	return coding_score;
}

/* EVAN check for overlaping scoring below */
/* EVAN change to model independant scoreing (range vs position for init, term, etc models)*/
/* EVAN change FFs to be consistant on real bases to score requirements. ie term is prohibited from having the acceptor overlap the initial paddign but the stop model can overlap the terminal padding */
static score_t zScoreEsngl (const zFeatureFactory *fac, zSfeature *exon) {
	score_t   total_score;
	score_t   initcon_score, termcon_score, cdscons_score; 
	score_t   initest_score, termest_score, cdsest_score; 
	score_t   initalign_score, termalign_score, cdsalign_score;	

	score_t   overlap_markov_score, rear_markov_score;
	score_t   overlap_coding_prob, coding_score;
	int       codon_num, exon_length, rear_length;
	int       frame = 0;
	coor_t    end_point, einit_cds_endpt;
	zSfeature signal_peptide;

	exon->cdscons_score = 0;	
	exon->cdsest_score  = 0;
	exon->align_score   = 0;	

	exon_length = exon->end - exon->start + 1;
	total_score = zGetUScore(fac->scan[START],exon->start)
		+ zGetUScore(fac->scan[STOP],exon->end - 2);
	
	/*this is an ineficent way to do this, come back and make it better later*/
	cdscons_score = 0;
	if(fac->conseq_enabled) {
		if(exon->start+6 > exon->end-3){
			cdscons_score = 0;
		}
		else{
			cdscons_score = zScoreRange(fac->conscan[CDSCONS],exon->start+6,exon->end-3);
		}
		
		if (fac->conscan[INITCONS]->model->type == LUT) { /* MC */
			initcon_score = zScoreRange(fac->conscan[INITCONS],exon->start-6,exon->start+5);
		}
		else { /* WWAM */
			initcon_score = zGetUScore(fac->conscan[INITCONS],exon->start);
		} 
		/* null model */
		initcon_score -= zScoreRange(fac->conscan[INTRONCONS],exon->start-6,exon->start+5);

		
		if (fac->conscan[TERMCONS]->model->type == LUT) { /* MC */
			termcon_score = zScoreRange(fac->conscan[TERMCONS],exon->end-2,exon->end+3);
		}
		else { /* WWAM */
			termcon_score = zGetUScore(fac->conscan[TERMCONS],exon->end-2);
		} 
		/* null model */
		termcon_score -= zScoreRange(fac->conscan[INTRONCONS],exon->end-2,exon->end+3);
		
		cdscons_score += initcon_score + termcon_score;
		
		/**! now adding cdscons_score to the exon score !**/
		exon->cdscons_score = cdscons_score;
		total_score        += cdscons_score;
	}

	cdsest_score = 0;
	if(fac->est_para_mode) {
		if(exon->start+6 > exon->end-3){
			cdsest_score = 0;
		}
		else{
			cdsest_score += zScoreRange(fac->estscan[CDSEST],exon->start+6,exon->end - 3);
		}
		
		initest_score = zScoreRange(fac->estscan[INITEST],exon->start-6,exon->start+5) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start-6,exon->start+5);
		termest_score = zScoreRange(fac->estscan[TERMEST],exon->end-2,exon->end+3) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->end-2,exon->end+3);
		cdsest_score += initest_score + termest_score;
		
		exon->cdsest_score  = cdsest_score;
		total_score        += cdsest_score;
	}

	if (fac->phylo_enabled) {
		cdsalign_score = 0;
		if (fac->alignscan[CDSALIGN]->model->type == BNTREE) {
			/* no frame dependence */
			cdsalign_score = zScoreRangeA(fac->alignscan[CDSALIGN],exon->start+6,exon->end-3);
		}
		else if (fac->alignscan[CDSALIGN]->model->type == BNTREE_CDS) {
			/* frame dependence */
			cdsalign_score = fac->alignscan[CDSALIGN]->scoref(fac->alignscan[CDSALIGN], exon);
		}
		cdsalign_score -= zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start+6,exon->end-3);
		
		initalign_score = zGetUScoreA(fac->alignscan[INITALIGN],exon->start);
		initalign_score -= zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start-6,exon->start+5);

		termalign_score = zGetUScoreA(fac->alignscan[TERMALIGN],exon->end-2);
		termalign_score -= zScoreRangeA(fac->alignscan[INTRONALIGN],exon->end-2,exon->end+3);

		cdsalign_score *= PSC;
		initalign_score *= PSC;
		termalign_score *= PSC;
		
		exon->align_score = cdsalign_score + initalign_score + termalign_score;
		
		total_score += (cdsalign_score + initalign_score + termalign_score);
	}
	  
#ifdef DEBUG
	exon->begin_score = zGetUScore(fac->scan[START],exon->start);
	exon->end_score = zGetUScore(fac->scan[STOP],exon->end-2);
#endif

	codon_num = exon_length/3-1;
	/* SIGNAL_PEPT_LENGTH see twinscan InitialExon.cpp and generalize later */
	if (codon_num > SIGNAL_PEPTIDE_LENGTH) {
		codon_num = SIGNAL_PEPTIDE_LENGTH;
	}
	end_point = exon->start + 3 + 3*codon_num - 1;
		
	if (end_point > exon->end - 3)  {
		end_point = exon->end - 3;
	}
	
	/* score overlap here */
	switch (exon->frame) {
	case 0: frame = 2; break;
	case 1: frame = 1; break;
	case 2: frame = 0; break;
	default: zDie("zScoreEsngl impossible");
	}
	
        /* Scoring 5' end of Esngl with EINIT_CDS model */
	overlap_markov_score = zGetRangeScore(fac->scan[EINIT_CDS],end_point,exon->start+5,frame); 
	
	if(exon_length <= 8) {
		overlap_markov_score = 0;
	}
	
	/* score rear here */
	
	/* pos-4 changed to pos-1 for genscan compatability */
	rear_length = ((exon->end-3) - (exon->start+2+3*codon_num));
	einit_cds_endpt = exon->start+EINIT_CDS_LENGTH-1;
	
	if (einit_cds_endpt > (exon->end-3)){
		einit_cds_endpt = exon->end-3;
	}
		
	if( rear_length  > 0) {
		/* Scoring 5' end of Esngl with EINIT_CDS model */
		rear_markov_score = zGetRangeScore(fac->scan[EINIT_CDS],einit_cds_endpt,exon->start+2+3*codon_num,frame);
		
		if (einit_cds_endpt < (exon->end-3)){
			/* Scoring the rest of Esngl with NONEINIT_CDS model */
		    rear_markov_score += zGetRangeScore(fac->scan[NONEINIT_CDS],exon->end-3,
												 einit_cds_endpt,frame);
		}
	}
	else {
		rear_markov_score = 0;
	}
	
	/* score signal peptide here */
	
	signal_peptide.start = exon->start + 3;
	signal_peptide.end   = exon->start + 3 + 3*codon_num;
	signal_peptide.score = fac->scan[SIG_PEP]->scoref(fac->scan[SIG_PEP], 
                                                          &signal_peptide);
		
	overlap_coding_prob = 0.2*pow(2.0,signal_peptide.score/10.0) 
		+ 0.8*pow(2.0,overlap_markov_score/10.0);
	
	coding_score = 10*zLog2(overlap_coding_prob)+rear_markov_score;
	
	if(exon_length > 5) {
		total_score += coding_score;				  
#ifdef DEBUG
		exon->content_score = coding_score;
	}
	else {
		exon->content_score = 0;
#endif
	}

	return total_score;
}

/* return true if estchar is seen in range from start->end, and false otherwise*/
static bool check_est_region(const zFeatureFactory* fac, coor_t start, coor_t end, char estchar){
	coor_t i = zGetEstseqNextPosBounded(fac->estsearch,estchar,start,end+1);
	if(i <= end){
		return true;
	}
	return false;
}

static coor_t find_prev_est(const zFeatureFactory* fac,coor_t pos,char estchar, coor_t min){
	return zGetEstseqPrevPosBounded(fac->estsearch,estchar,pos,min);
	/*coor_t i = pos;
	while(i > min){
		if(zGetEstseqSeq(fac->estseq,i) == estchar){
			return i;
		}
		i--;
	}
	return min;	*/
}

static coor_t find_next_est(const zFeatureFactory* fac,coor_t pos,char estchar, coor_t max){
	return zGetEstseqNextPosBounded(fac->estsearch,estchar,pos,max);
	/*coor_t i = pos;
	while(i < max){
		if(zGetEstseqSeq(fac->estseq,i) == estchar){
			return i;
		}
		i++;
	}
	return max;*/
}

int not_in_intron_check(const zFeatureFactory *fac, coor_t start, coor_t end, zStrIdx name) { 
	/* return -1 if [start, end] is not an exon according to the EST alignment sequence.
       This includes two cases: 
	      1. contain a '2', 
	      2. [start, end] can extend at lease one more base in either or both ends

	   return 1 otherwise
	*/

	int j;

	if(fac->estseq == NULL) return 1;

	if(strncmp(zStrIdx2Char(name), "Exon", 4) == 0 &&
	   (zGetEstseqSeq(fac->estseq,start - 1) == '1' || 
		zGetEstseqSeq(fac->estseq,end + 1) =='1')){ 
		return -1;
	}
	else if(strncmp(zStrIdx2Char(name), "Eterm", 5) == 0 &&
			zGetEstseqSeq(fac->estseq,start -1) == '1'){ 
		return -1;
	}
	else if(strncmp(zStrIdx2Char(name), "Einit", 5) == 0 &&
			zGetEstseqSeq(fac->estseq,end+1) =='1' ){
		return -1;
	}
		
	for (j = (int)start; j <= (int) end; j++) { 
		if(zGetEstseqSeq(fac->estseq,j) == '2') {
			return -1;
		} 
	}  
	
	return 1;
}

int overlap_exon_or_intron(const zFeatureFactory *fac, coor_t start, coor_t end) { 
	/* return 1 if [start, end] overlap exon or intron region according 
        to the EST alignment sequence.
       return -1 otherwise
	*/

	int j;

	if(fac->estseq == NULL) return 1;
	
	for (j = (int)start; j <=(int) end; j++) { 
		if(zGetEstseqSeq(fac->estseq,j) == '2' ||
		   zGetEstseqSeq(fac->estseq,j) == '1' ) {
			return 1;
		} 
	}

	return -1;
}

static void zMakeEpasFrom3 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	int i, begin, end;
	zSfeature* exon;
	
	if (zGetUScore(fac->scan[START],pos+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
		end   = pos - 4; /* prevents overlap with start codon model */
		if (pos < PADDING + MAX_EPA_LENGTH - 1){
			begin = PADDING;
		}
		else{
			begin = pos - MAX_EPA_LENGTH + 1;
		}
		
		for (i = begin; i < end; i++) { 
			exon = zSFListAppend(sfl);			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name   = fac->type;
			exon->start  = i;
			exon->end    = pos;
			exon->lfrag  = 0;
			exon->rfrag  = 0;
			exon->frame  = 0;
			
			exon->score  = fac->scoref(fac, exon);
			
		}
	}
}

static void zMakeEpasFrom5 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i;
	zSfeature* exon;
	if(pos < PADDING) return;
	/* min size 6 */
	for(i = pos + 5; i < pos + MAX_EPA_LENGTH; i++){
		if(zGetUScore(fac->scan[START],i+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE){
			exon = zSFListAppend(sfl);			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name   = fac->type;
			exon->start  = pos;
			exon->end    = i;
			exon->lfrag  = 0;
			exon->rfrag  = 0;
			exon->frame  = 0;
			
			exon->score  = fac->scoref(fac, exon);
		}
	}
}


/*EVAN all zMakeXXXFromInside functions need to be tested post code merge before use */

static void zMakeEpasFromInside (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t begin, end;
	zSfeature* exon;
	if (pos < PADDING) return;
	
	for(end = pos; end < pos + MAX_EPA_LENGTH; end++){
		if (zGetUScore(fac->scan[START],end+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
			begin = end - MAX_EPA_LENGTH + 1;
			if(end  < PADDING + MAX_EPA_LENGTH - 1){
				begin = PADDING;
			}
			else{
				begin = end - MAX_EPA_LENGTH + 1;
			}
			/* min size 6 */
			for(;begin < end - 4; begin++){ 
				exon = zSFListAppend(sfl);			
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name   = fac->type;
				exon->start  = begin;
				exon->end    = end;
				exon->lfrag  = 0;
				exon->rfrag  = 0;
				exon->frame  = 0;
				
				exon->score  = fac->scoref(fac, exon);
				
			}
		}
	}
}

static void zMakeEpsFrom3 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, begin, end;
	zSfeature* exon;
	
	if (zGetUScore(fac->scan[DON],pos+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
		end   = pos - 1;  /* allows enough space for donor model */
		if (pos < PADDING + MAX_EP_LENGTH - 1){
			begin = PADDING;
		}
		else{
			begin = pos - MAX_EP_LENGTH + 1;
		}
		
		/* min size 3 */
		for (i = begin; i < end; i++) {
			exon = zSFListAppend(sfl);			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name  = fac->type;
			exon->start = i;
			exon->end   = pos;
			exon->lfrag = 0;
			exon->rfrag = 0;
			exon->frame = 0;
				
			exon->score = fac->scoref(fac, exon);
		}
	}
	return;
}

static void zMakeEpsFrom5 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i;
	zSfeature* exon;
	if (pos < PADDING) return;

	/* min size 3 */
	for(i = pos + 2; i < pos + MAX_EP_LENGTH;i++){
		if(zGetUScore(fac->scan[DON],i+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE){
			exon = zSFListAppend(sfl);			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name   = fac->type;
			exon->start  = pos;
			exon->end    = i;
			exon->lfrag  = 0;
			exon->rfrag  = 0;
			exon->frame  = 0;
			
			exon->score  = fac->scoref(fac, exon);
		}
	}
}

static void zMakeEpsFromInside (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t begin, end;
	zSfeature* exon;
	if (pos < PADDING) return;
	
	for(end = pos; end < pos + MAX_EP_LENGTH; end++){
		if (zGetUScore(fac->scan[DON],end+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
			begin = end - MAX_EP_LENGTH + 1;
			if(end < PADDING + MAX_EP_LENGTH - 1){
				begin = PADDING;
			}
			else{
				begin = end - MAX_EP_LENGTH + 1;
			}
			/* min size 2*/
			for(;begin < end; begin++){ 
				exon = zSFListAppend(sfl);			
				exon->group  = NULL;
				exon->strand = fac->strand;
				
				exon->name   = fac->type;
				exon->start  = begin;
				exon->end    = end;
				exon->lfrag  = 0;
				exon->rfrag  = 0;
				exon->frame  = 0;
				
				exon->score  = fac->scoref(fac, exon);
			}
		}
	}
}

static void zMakeEasFrom3 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, begin, end;
	
	zSfeature* exon;
	
	if (zGetUScore(fac->scan[START],pos+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
	    /* min size 9 */
		end   = pos - 7; /* prevents overlap with start codon model    */
		                 /* and allows enough space for acceptor model */

		/* 40 here allows enough space to score acceptor */
		/* model at the 5' end of the exon->              */
		if (pos < (PADDING+40)  + MAX_EA_LENGTH - 1){
			begin = PADDING + 40;
		}
		else{
			begin = pos - MAX_EA_LENGTH + 1;			
		}
	    for (i = begin; i < end; i++) { 
			if (zGetUScore(fac->scan[ACC],i-1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
				continue;
			exon = zSFListAppend(sfl);			
			exon->group  = NULL;
			exon->strand = fac->strand;

			exon->name  = fac->type;
			exon->start = i;
			exon->end   = pos;
			exon->frame = 0;
			exon->lfrag = 0;
			exon->rfrag = 0;
		
			exon->score = fac->scoref(fac, exon);
		}
	}	
}

static void zMakeEasFrom5 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i;
	
	zSfeature* exon;

	/* 40 here allows enough space to score acceptor */
	/* model at the 5' end of the exon->              */
	if (pos < PADDING + 40) return;

	if (zGetUScore(fac->scan[ACC],pos-1) > MINIMUM_TWINSCAN_SIGNAL_SCORE){
		/* min size 9 */
		for(i = pos + 8; i < pos + MAX_EA_LENGTH; i++){
			if(zGetUScore(fac->scan[START],i+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE){
				exon = zSFListAppend(sfl);			
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name   = fac->type;
				exon->start  = pos;
				exon->end    = i;
				exon->lfrag  = 0;
				exon->rfrag  = 0;
				exon->frame  = 0;
				
				exon->score  = fac->scoref(fac, exon);
			}
		}
	}
	return;
}

static void zMakeEasFromInside (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t begin, end;
	
	zSfeature* exon;
	if (pos < PADDING) return;
	
	for(end = pos; end < pos + MAX_EA_LENGTH; end++){
		if (zGetUScore(fac->scan[START],end+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
            /* 40 here allows enough space to score acceptor */
            /* model at the 5' end of the exon->              */
			if(end  < PADDING + 40 + MAX_EA_LENGTH - 1){
				begin = PADDING + 40;
			}
			else{
				begin = end - MAX_EA_LENGTH + 1;
			}

			/* min size 8*/
			for(;begin < end - 6; begin++){ 
				if(zGetUScore(fac->scan[ACC],begin-1) <= 
				   MINIMUM_TWINSCAN_SIGNAL_SCORE)
                    continue;
				exon = zSFListAppend(sfl);			
				
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name   = fac->type;
				exon->start  = begin;
				exon->end    = end;
				exon->lfrag  = 0;
				exon->rfrag  = 0;
				exon->frame  = 0;
				
				exon->score  = fac->scoref(fac, exon);
			}
		}
	}
}

static void zMakeEncsFrom3 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
    coor_t i, begin, end;
	zSfeature* exon;
	if (pos < PADDING) return;
	
	if (zGetUScore(fac->scan[DON],pos+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {	
		/* min size 6*/
		end = pos - 4; /* allows enough space for donor and acceptor */
		               /* models and prevents overlap between them   */

		/* 40 here allows enough space to score acceptor */
		/* model at the 5' end of the exon->              */
		if (pos < (PADDING+40) + MAX_ENC_LENGTH - 1){
			begin = PADDING + 40;
		}
		else{
			begin = pos - MAX_ENC_LENGTH + 1;
		}

		for (i = begin; i < end; i++) { 
			if (zGetUScore(fac->scan[ACC],i-1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE) 
				continue;
			exon = zSFListAppend(sfl);			
			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name   = fac->type;
			exon->start  = i;
			exon->end    = pos;
			exon->lfrag  = 0;
			exon->rfrag  = 0;
			exon->frame  = 0;
			
			exon->score  = fac->scoref(fac, exon);
		}
	}
}

static void zMakeEncsFrom5 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
    coor_t i;
	zSfeature* exon;

	/* 40 here allows enough space to score acceptor */
	/* model at the 5' end of the exon->              */
	if (pos < PADDING + 40) return;

	if (zGetUScore(fac->scan[ACC],pos-1) > MINIMUM_TWINSCAN_SIGNAL_SCORE){
		/* min size 6 */
		for(i = pos + 5; i < pos + MAX_ENC_LENGTH;i++){
			if(zGetUScore(fac->scan[DON],i+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE){
				exon = zSFListAppend(sfl);			
				
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name   = fac->type;
				exon->start  = pos;
				exon->end    = i;
				exon->lfrag  = 0;
				exon->rfrag  = 0;
				exon->frame  = 0;
				
				exon->score  = fac->scoref(fac, exon);
			}
		}
	}
}

static void zMakeEncsFromInside (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
    coor_t begin, end;
	zSfeature* exon;
	if (pos < PADDING) return;
	
	for(end = pos; end < pos + MAX_ENC_LENGTH; end++){
		if (zGetUScore(fac->scan[DON],end+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
			begin = end - MAX_ENC_LENGTH + 1;
            /* 40 here allows enough space to score acceptor */
            /* model at the 5' end of the exon->              */
			if(end < PADDING + 40 + MAX_ENC_LENGTH - 1){
				begin = PADDING + 40;
			}
			else{
				begin = end - MAX_ENC_LENGTH + 1;
			}
			/* min size 6*/
			for(;begin < end - 4; begin++){ 
				if(zGetUScore(fac->scan[ACC],begin-1) <= 
				   MINIMUM_TWINSCAN_SIGNAL_SCORE)
                    continue;
				exon = zSFListAppend(sfl);			
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name   = fac->type;
				exon->start  = begin;
				exon->end    = end;
				exon->lfrag  = 0;
				exon->rfrag  = 0;
				exon->frame  = 0;
				
				exon->score  = fac->scoref(fac, exon);
			}
		}
	}
}

static score_t zScoreEpa(const zFeatureFactory *fac, zSfeature *exon) {
 	score_t total_score;
	score_t content_score, cons_score;
	score_t est_score;
	score_t align_score;

	exon->cdscons_score = 0;	
	exon->cdsest_score  = 0;
	exon->align_score   = 0;	

	content_score = zScoreRange(fac->scan[EP],exon->start, exon->end-6);

	total_score = content_score;

	if (fac->conseq_enabled) {
		cons_score = zScoreRange(fac->conscan[UTRCONS],exon->start,exon->end-6) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->start,exon->end-6);

		total_score += cons_score;
	}

	if(fac->est_para_mode) {
		est_score = zScoreRange(fac->estscan[EPAEST],exon->start,exon->end - 6) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start,exon->end - 6);		
		total_score += est_score;
	}

	if (fac->phylo_enabled) {
		align_score = zScoreRangeA(fac->alignscan[EPA_ALIGN],exon->start,exon->end-6) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start,exon->end - 6);		
		
		align_score *= PSC;
		
		exon->align_score = align_score;

		total_score += align_score;
	}

	return total_score;
}

static score_t zScoreEp(const zFeatureFactory *fac, zSfeature *exon) {
	score_t total_score;
	score_t don_score, content_score, cons_score;
	score_t est_score;
	score_t align_score;

	exon->cdscons_score = 0;	
	exon->cdsest_score  = 0;
	exon->align_score   = 0;	

	don_score = zGetUScore(fac->scan[DON],exon->end + 1);
	content_score  = zScoreRange(fac->scan[EP],exon->start, exon->end-3);
	
	total_score = don_score + content_score;
	
	if (fac->conseq_enabled) {
		cons_score = zScoreRange(fac->conscan[UTRCONS],exon->start,exon->end-3) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->start,exon->end-3);

		cons_score += zGetUScore(fac->conscan[UTRDONCONS],exon->end+1) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->end-2,exon->end+6);
				
		total_score += cons_score;
	}

	if(fac->est_para_mode) {
		
		est_score = zScoreRange(fac->estscan[EPEST],exon->start,exon->end-3) -
			zScoreRange(fac->estscan[INTERGENICEST],exon->start,exon->end-3);
		
		est_score += zGetUScore(fac->estscan[UTRDONEST],exon->end+1) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->end-2,exon->end+6);
		
		total_score += est_score;
	}
	
	if (fac->phylo_enabled) {
		align_score = zScoreRangeA(fac->alignscan[EP_ALIGN],exon->start,exon->end-3) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start,exon->end-3);
		
		align_score += zGetUScoreA(fac->alignscan[UTRDONALIGN],exon->end+1) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->end-2,exon->end+6);
		
		align_score *= PSC;
		
		exon->align_score = align_score;
		
		total_score += align_score;
	}
	
	return total_score;
}

static score_t zScoreEa(const zFeatureFactory *fac, zSfeature *exon) {
	score_t total_score;
	score_t acc_score, content_score, cons_score;
	score_t est_score;
	score_t align_score;
	
	exon->cdscons_score = 0;	
	exon->cdsest_score  = 0;
	exon->align_score   = 0;	

	acc_score  = zGetUScore(fac->scan[ACC],exon->start - 1);
	content_score  = zScoreRange(fac->scan[EA],exon->start+3,exon->end-6);

	total_score = acc_score + content_score;

	if (fac->conseq_enabled) {
		cons_score = zScoreRange(fac->conscan[UTRCONS],exon->start+3,exon->end-6) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->start+3,exon->end-6);

		cons_score += zGetUScore(fac->conscan[UTRACCCONS],exon->start-1) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->start-40,exon->start+2);
			
		total_score += cons_score;
	}
	
	if(fac->est_para_mode) {
		est_score = zScoreRange(fac->estscan[EAEST],exon->start+3,exon->end-6) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start+3,exon->end-6);

		est_score += zGetUScore(fac->estscan[UTRACCEST],exon->start-1) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start-40,exon->start+2);
		
		total_score += est_score;
	}
	
	if (fac->phylo_enabled) {
		align_score = zScoreRangeA(fac->alignscan[EA_ALIGN],exon->start+3,exon->end-6) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start+3,exon->end-6);
		
		align_score += zGetUScoreA(fac->alignscan[UTRACCALIGN],exon->start-1) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start-40,exon->start+2);
		
		align_score *= PSC;
		
		exon->align_score = align_score;
		
		total_score += align_score;
	}
	
	return total_score;
}

static score_t zScoreEnc(const zFeatureFactory *fac, zSfeature *exon) {
	score_t total_score;
	score_t acc_score, don_score, content_score, cons_score;
	score_t est_score;
	score_t align_score;

	exon->cdscons_score = 0;	
	exon->cdsest_score  = 0;
	exon->align_score   = 0;	

	acc_score  = zGetUScore(fac->scan[ACC],exon->start - 1);
	don_score = zGetUScore(fac->scan[DON],exon->end + 1);
	content_score  = zScoreRange(fac->scan[EA],exon->start+3,exon->end-3);

	total_score = acc_score + don_score + content_score;

	if (fac->conseq_enabled) {
		cons_score = zScoreRange(fac->conscan[UTRCONS],exon->start+3,exon->end-3) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->start+3,exon->end-3);
		
		cons_score += zGetUScore(fac->conscan[UTRACCCONS],exon->start-1) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->start-40,exon->start+2);
		
		cons_score += zGetUScore(fac->conscan[UTRDONCONS],exon->end+1) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->end-2,exon->end+6);
		
		total_score += cons_score;
	}

	if(fac->est_para_mode) {
		est_score = zScoreRange(fac->estscan[ENCEST],exon->start+3,exon->end - 3) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start+3,exon->end - 3);
		
		est_score += zGetUScore(fac->estscan[UTRACCEST],exon->start-1) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start-40,exon->start+2);
		
		est_score += zGetUScore(fac->estscan[UTRDONEST],exon->end+1) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->end-2,exon->end+6);
		
		total_score += est_score;
	}

	if (fac->phylo_enabled) {
		align_score = zScoreRangeA(fac->alignscan[ENC_ALIGN],exon->start+3,exon->end-3) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start+3,exon->end-3);
		
		align_score += zGetUScoreA(fac->alignscan[UTRACCALIGN],exon->start-1) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start-40,exon->start+2);
		
		align_score += zGetUScoreA(fac->alignscan[UTRDONALIGN],exon->end+1) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->end-2,exon->end+6);
		
		align_score *= PSC;

		exon->align_score = align_score;

		total_score += align_score;
	}
	
	return total_score;
}

static void zMakeEsnglsFrom3 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, begin, end;
	int frame;
	zSfeature* exon;
	
	if(pos < PADDING) return;
  
	if (zGetUScore(fac->scan[STOP],pos-2) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
		frame = (pos+1) % 3; 
		/* Esngl is in the same frame as the stop codon */
		/* min size 9 (start + stop + 1 codon) */
		end   = pos - 8;
		begin = zGetStopSeqPrevPos(fac->fstop,pos-5)+3;
		if(fac->estseq != NULL && !fac->est_para_mode) {
			/*EVAN checking all but the last 2 bases for compatability.
			  This can't be right since eterm checks all but the last 3 */
			begin = find_prev_est(fac,pos-2,'2',begin-1)+1;
		}
		for (i = end; i >= begin; i -= 3) {
			if (zGetUScore(fac->scan[START],i) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
				continue;

			/*			if(fac->estseq != NULL && !fac->est_para_mode) {
				flag = 1;
				flag = not_in_intron_check(fac, i, pos-2, fac->type);
				if(flag <  0) continue;
				}*/

			exon = zSFListAppend(sfl);			

			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name  = fac->type;
			exon->start = i;
			exon->end   = pos;
			exon->lfrag = 0;
			exon->rfrag = 0;
			exon->frame = frame;
			
			exon->score = fac->scoref(fac, exon);
		}
	}
}	

static void zMakeEsnglsFrom5 (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	int frame;
	coor_t end;
	zSfeature* exon;
 
	if (pos < PADDING) return;
  	if (zGetUScore(fac->scan[START],pos) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
		frame = pos % 3; /* Esngl is in the same frame as the start codon */
		end = zGetStopSeqNextPos(fac->fstop,pos)+2;

		/* minimum of one non-start/stop codon */
		if(end - pos <= 5) return;

		if ((zGetUScore(fac->scan[STOP],end-2) <= MINIMUM_TWINSCAN_SIGNAL_SCORE))
			return;
		
		if(fac->estseq != NULL && !fac->est_para_mode) {
			/*EVAN checking all but the last 2 bases for compatability.
			  This can't be right since eterm checks all but the last 3 */
			if(check_est_region(fac,pos,end-2,'2')){
				return;
			}
			/*flag = 1;
			  flag = not_in_intron_check(fac, pos, end-2, fac->type);
			  if(flag <  0) return;*/
		}
  
		exon = zSFListAppend(sfl);			
		exon->group  = NULL;
		exon->strand = fac->strand;
		
		exon->name  = fac->type;
		exon->start = pos;
		exon->end   = end;
		exon->lfrag = 0;
		exon->rfrag = 0;
		exon->frame = frame;
		
		exon->score = fac->scoref(fac, exon);
	}
}	

static void zMakeEsnglsFromInside (const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, frame, begin, end, max_begin;
	zSfeature* exon;
 
	if (pos < PADDING) return;
  
	/* for each possible frame (relative to start of sequence)*/
	for(frame = 0; frame <= 2; frame++){
		end = zGetStopSeqNextPos(fac->fstop,pos-frame)+2;
		if ((zGetUScore(fac->scan[STOP],end-2) <= MINIMUM_TWINSCAN_SIGNAL_SCORE))
			continue;
		begin = zGetStopSeqPrevPos(fac->fstop,end-2-3)+3;
		max_begin = pos;
		/* must be at least 3 codons (start - interal - stop) long */
		if(end-8 < max_begin){
			max_begin = end-8;
		}
		
		if(fac->estseq != NULL && !fac->est_para_mode) {
			/*EVAN checking all but the last 2 bases for compatability.
			  This can't be right since eterm checks all but the last 3 */
			if(check_est_region(fac,max_begin-1,end-2,'2')){
				/*check next frame*/
				continue;
			}
			/*EVAN stupidly checking same region 3 times*/
			begin = find_prev_est(fac,max_begin,'2',begin-1)+1;
		}
		
		for(i = max_begin;i >= begin; i-=3){
			if (zGetUScore(fac->scan[START],i) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
				continue;
			
			/*if(fac->estseq != NULL && !fac->est_para_mode) {
			  flag = 1;
			  flag = not_in_intron_check(fac, i, end-2, fac->type);
			  if(flag <  0) continue;
			  }*/
		
			exon = zSFListAppend(sfl);			
			
			
			exon->group  = NULL;
			exon->strand = fac->strand;

 			exon->name  = fac->type;
			exon->start = i;
			exon->end   = end;
			exon->lfrag = 0;
			exon->rfrag = 0;
			exon->frame = i % 3;
			
			exon->score = fac->scoref(fac, exon);
		}
	}
}	

static score_t zScoreEterm(const zFeatureFactory *fac, zSfeature *exon) {
	score_t total_score;
	score_t cdscons_score, termcon_score, acccon_score;
	score_t cdsest_score, termest_score, accest_score;
	score_t cdsalign_score, termalign_score, accalign_score;
	score_t acc_score, stop_score, cds_score;

	exon->cdscons_score = 0;	
	exon->cdsest_score  = 0;
	exon->align_score   = 0;	

	acc_score  = zGetUScore(fac->scan[ACC],exon->start - 1);
	stop_score = zGetUScore(fac->scan[STOP],exon->end - 2);
	cds_score  = zScoreCoding(fac, exon);
	
	total_score = stop_score + acc_score + cds_score;
	
#ifdef DEBUG
	exon->begin_score   = acc_score;
	exon->end_score     = stop_score;
	exon->content_score = cds_score;
#endif
	
	if(fac->conseq_enabled) {
		/* this is an ineficent way to do this, come back and make it better later */
		acccon_score = zGetUScore(fac->conscan[ACCCONS],exon->start - 1) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->start-40,exon->start + 2);

		if (fac->conscan[TERMCONS]->model->type == LUT) { /* MC */
			termcon_score  = zScoreRange(fac->conscan[TERMCONS],exon->end-2,exon->end+3);
		}
		else { /* WWAM */
			termcon_score = zGetUScore(fac->conscan[TERMCONS],exon->end-2);
		} /* null model */
		termcon_score -= zScoreRange(fac->conscan[INTRONCONS],exon->end-2,exon->end+3);
		
		/* this is an ineficent way to do this, come back and make it better later */
		if(exon->start+3 > exon->end-3){
			cdscons_score = 0;
		}
		else{
			cdscons_score = zScoreRange(fac->conscan[CDSCONS],exon->start+3,exon->end-3);
		}

		cdscons_score += termcon_score + acccon_score;
		
		/**! now adding cdscons_score to the exon score !**/
		exon->cdscons_score = cdscons_score;
		total_score        += cdscons_score;
	}
	
	if(fac->est_para_mode) {
		accest_score = zGetUScore(fac->estscan[ACCEST],exon->start - 1) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start-40,exon->start + 2);
		
		termest_score = zScoreRange(fac->estscan[TERMEST],exon->end-2,exon->end+3) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->end-2,exon->end+3);
		
		if(exon->start+3 > exon->end-3){
			cdsest_score = 0;
		}
		else{
			cdsest_score = zScoreRange(fac->estscan[CDSEST],exon->start+3,exon->end-3);
		}
		
		cdsest_score += termest_score + accest_score;
		exon->cdsest_score = cdsest_score;
		total_score        += cdsest_score;
	}

	if (fac->phylo_enabled) {
		cdsalign_score = 0;
		if (fac->alignscan[CDSALIGN]->model->type == BNTREE) {
			/* no frame dependence */
			cdsalign_score = zScoreRangeA(fac->alignscan[CDSALIGN],exon->start+3,exon->end-3);
		}
		else if (fac->alignscan[CDSALIGN]->model->type == BNTREE_CDS) {
			/* frame dependence */
			cdsalign_score = fac->alignscan[CDSALIGN]->scoref(fac->alignscan[CDSALIGN], exon);
		}
		cdsalign_score -= zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start+3,exon->end-3);
		
		accalign_score = zGetUScoreA(fac->alignscan[ACCALIGN],exon->start-1) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start-40,exon->start+2);
		
		termalign_score = zGetUScoreA(fac->alignscan[TERMALIGN],exon->end-2) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->end-2,exon->end+3);
		
		cdsalign_score *= PSC;
		accalign_score *= PSC;
		termalign_score *= PSC;
		
		exon->align_score = cdsalign_score + accalign_score + termalign_score;
		
		total_score += (cdsalign_score + accalign_score + termalign_score);
	}
	
	return total_score;
}

static void zMakeEtermsFrom3(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, frame, begin, end;
	
	zSfeature* exon;
	if (pos < PADDING) return;
	
	if ((zGetUScore(fac->scan[STOP],pos-2) >  MINIMUM_TWINSCAN_SIGNAL_SCORE) &&
		(pos-2 < fac->scan[ACC]->seq->length-PADDING)) {
		frame = (pos+1) % 3; /* Eterm frame is the same as stop */
		end = pos - 5;
		begin = zGetStopSeqPrevPos(fac->fstop,end)+1;
		
		if(fac->estseq != NULL && !fac->est_para_mode){
			begin = find_prev_est(fac,pos-3,'2',begin-1)+1;
		}
		for (i = end; i >= begin; i--) { /* to allow 3 bp terminal exons */
			if ((zGetUScore(fac->scan[ACC],i-1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE) || 
				(i <= PADDING+43)) continue;
			/*if(fac->estseq != NULL && !fac->est_para_mode){
				flag = 1;
				flag = not_in_intron_check(fac, i, pos - 3 , fac->type);
				if(flag < 0) continue;
			}*/
			if(fac->estseq != NULL && !fac->est_para_mode){
				if(zGetEstseqSeq(fac->estseq,i-1) == '1'){
					continue;
				}
			}
			
			exon = zSFListAppend(sfl);			
			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name  = fac->type;
			exon->start = i;
			exon->end   = pos;
			exon->frame = frame % 3;
			exon->lfrag = (exon->end - exon->start + 1) % 3;
			exon->rfrag = 0;

			exon->score = fac->scoref(fac, exon);
		}
	}
}

static void zMakeEtermsFrom5(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t stop[3];
	coor_t end,max_end;
	int loh,frame;
	
	zSfeature* exon;
	if (pos < PADDING) return;

	if ((zGetUScore(fac->scan[ACC],pos-1) > MINIMUM_TWINSCAN_SIGNAL_SCORE)
		&& (pos > PADDING+43)){
		if(fac->estseq != NULL && !fac->est_para_mode){
			if(zGetEstseqSeq(fac->estseq,pos-1) == '1'){
				return;
			}
		}
		
		/*find the highest coored we may use for any loh */
		max_end = 0;
		for(loh = 0; loh < 3; loh++){
			stop[loh] = zGetStopSeqNextPos(fac->fstop,pos+loh)+2;
			if(stop[loh] > max_end){
				max_end = stop[loh];
			}
		}
		
		if(fac->estseq != NULL && !fac->est_para_mode){
			/* find the first non exon est pos in the range we are looking at */
			max_end = find_next_est(fac,pos,'2',max_end+1)-1+3;
		}
		
		for(loh = 0; loh < 3; loh++){
			end = stop[loh];
			if(end > max_end) continue; /* non exon est in range */
				
			if(end - pos + 1 < 6) continue; /* Minimum 3bp term exon */
			frame = (end-2) % 3; /* Eterm frame is the same as stop */
			
			if ((zGetUScore(fac->scan[STOP],end-2) <= MINIMUM_TWINSCAN_SIGNAL_SCORE) || 
				(end-2 >= fac->scan[ACC]->seq->length-PADDING)) continue;
			
			/*if(fac->estseq != NULL && !fac->est_para_mode){
			  flag = 1;
			  flag = not_in_intron_check(fac, pos, end - 3 , fac->type);
			  if(flag < 0) continue;
			  }*/
			
			exon = zSFListAppend(sfl);			
			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name  = fac->type;
			exon->start = pos;
			exon->end   = end;
			exon->frame = frame;
			exon->lfrag = loh;
			exon->rfrag = 0;
			
			exon->score = fac->scoref(fac, exon);
		}
	}
}

static void zMakeEtermsFromInside(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, begin, end, max_begin,max_end;
	coor_t stop[3];
	int frame;
	
	zSfeature* exon;
	if (pos < PADDING) return;

	max_end = 0;
	for(frame = 0; frame <= 2; frame++){
		stop[frame] = zGetStopSeqNextPos(fac->fstop,pos-frame)+2;
		if(stop[frame] > max_end){
			max_end = stop[frame];
		}
	}
	if(fac->estseq != NULL && !fac->est_para_mode){
		max_end = find_next_est(fac,pos,'2',max_end+1)-1+3;
	}

	for(frame = 0; frame <= 2; frame++){
		end = stop[frame];
		if(end > max_end) continue; /* non exon est in range */
		if ((zGetUScore(fac->scan[STOP],end-2) <= MINIMUM_TWINSCAN_SIGNAL_SCORE) || 
			(end-2 >= fac->scan[ACC]->seq->length-PADDING-10)) continue;
		begin = zGetStopSeqPrevPos(fac->fstop,end-2-3)+1;
		
		if(begin < PADDING+43) begin = PADDING+43;

		max_begin = pos;
		/* must be at least 6 bp long */
		if(end-5 < max_begin) max_begin = end-5;

		if(fac->estseq != NULL && !fac->est_para_mode){
			/*EVAN stupidly checking same region 3 times (once pre frame)*/
			begin = find_prev_est(fac,pos,'2',begin-1)+1;
		}

		for (i = max_begin; i >= begin; i--) {
			if ((zGetUScore(fac->scan[ACC],i-1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
				|| (i < PADDING+43)) continue;
			/*if(fac->estseq != NULL && !fac->est_para_mode){
				flag = 1;
				flag = not_in_intron_check(fac, i, end-3 , fac->type);
				if(flag < 0) continue;
			}*/
			if(fac->estseq != NULL && !fac->est_para_mode){
				if(zGetEstseqSeq(fac->estseq,i-1) == '1'){
					continue;
				}
			}

			exon = zSFListAppend(sfl);			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name  = fac->type;
			exon->start = i;
			exon->end   = end;
			exon->frame = (end-2) % 3; /* Eterm frame is the same as stop */
			exon->lfrag = (exon->end - exon->start + 1) % 3;
			exon->rfrag = 0;
			
			exon->score = fac->scoref(fac, exon);
		}
	}
}

static score_t zScoreEinit(const zFeatureFactory *fac, zSfeature *exon) {
	int       exon_length, codon_num, frame, rear_length;
	coor_t    end_point, einit_cds_endpt;
	score_t   total_score;
	score_t   start_score, don_score;
	score_t   initcon_score, doncon_score;
	score_t   initest_score, donest_score;
	score_t   initalign_score, donalign_score, cdsalign_score;
	score_t   overlap_markov_score, rear_markov_score;
	score_t   overlap_coding_prob, coding_score;
	zSfeature signal_peptide;
	
	exon->cdscons_score = 0;	
	exon->cdsest_score  = 0;
	exon->align_score   = 0;	

	start_score = zGetUScore(fac->scan[START],exon->start);
	don_score   = zGetUScore(fac->scan[DON],exon->end+1);

	total_score = start_score + don_score;

#ifdef DEBUG
	exon->begin_score = start_score;
	exon->end_score = don_score;
#endif

	exon_length = exon->end - exon->start + 1;

	exon->cdscons_score = 0;
	if (fac->conseq_enabled) {
		if (fac->conscan[INITCONS]->model->type == LUT) { /* MC */
			initcon_score = zScoreRange(fac->conscan[INITCONS],exon->start-6,exon->start+5);
		}
		else { /* WWAM */
			initcon_score = zGetUScore(fac->conscan[INITCONS],exon->start);
		} /* null model */
		initcon_score -= zScoreRange(fac->conscan[INTRONCONS],exon->start-6,exon->start+5);
				
		/*this is an ineficent way to do this, come back and make it better later*/
		doncon_score = zGetUScore(fac->conscan[DONCONS],exon->end+1) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->end-2,exon->end+6);
		
		/* this is an ineficent way to do this, come back and make it better later */
		if(exon->start+6 > exon->end-3){
			exon->cdscons_score = 0;
		}
		else{
			exon->cdscons_score = zScoreRange(fac->conscan[CDSCONS],exon->start+6,exon->end-3);
		}

		exon->cdscons_score += doncon_score + initcon_score;
		
		total_score += exon->cdscons_score;
	}

	exon->cdsest_score = 0;
	if (fac->est_para_mode) {
		initest_score = zScoreRange(fac->estscan[INITEST],exon->start-6,exon->start+5) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start-6,exon->start+5);
		
		donest_score = zGetUScore(fac->estscan[DONEST],exon->end+1) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->end-2,exon->end+6);
		
		if(exon->start+6 > exon->end-3){
			exon->cdsest_score = 0;
		}
		else{
			exon->cdsest_score = zScoreRange(fac->estscan[CDSEST],exon->start+6,exon->end-3);
		}

		exon->cdsest_score += donest_score + initest_score;
		
		total_score += exon->cdsest_score;
	}

	if (fac->phylo_enabled) {
		cdsalign_score = 0;
		if (fac->alignscan[CDSALIGN]->model->type == BNTREE) {
			/* no frame dependence */
			cdsalign_score = zScoreRangeA(fac->alignscan[CDSALIGN],exon->start+6,exon->end-3);
		}
		else if (fac->alignscan[CDSALIGN]->model->type == BNTREE_CDS) {
			/* frame dependence */
			cdsalign_score = fac->alignscan[CDSALIGN]->scoref(fac->alignscan[CDSALIGN], exon);
		}
		cdsalign_score -= zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start+6,exon->end-3);

		initalign_score = zGetUScoreA(fac->alignscan[INITALIGN],exon->start) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start-6,exon->start+5);

		donalign_score = zGetUScoreA(fac->alignscan[DONALIGN],exon->end+1) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->end-2,exon->end+6);
		
		cdsalign_score *= PSC;
		initalign_score *= PSC;
		donalign_score *= PSC;

		exon->align_score = cdsalign_score + initalign_score + donalign_score;

		total_score += (cdsalign_score + initalign_score + donalign_score);
	}
	
	codon_num = exon_length/3-1;
	/* SIGNAL_PEPT_LENGTH see twinscan InitialExon->cpp and generalize later */
	if (codon_num > SIGNAL_PEPTIDE_LENGTH) {
		codon_num = SIGNAL_PEPTIDE_LENGTH;
	}
	end_point = exon->start + 3 + 3*codon_num - 1;
	
	if (end_point > exon->end - 3)  {
		end_point = exon->end - 3;
	}
		
	/* score overlap here */
	
	frame = (exon->start) % 3;
	
	switch (frame) {
	case 0: frame = 2; break;
	case 1: frame = 1; break;
	case 2: frame = 0; break;
	default: zDie("zMakeExons impossible");
	}

        /* Scoring 5' end of Einit with EINIT_CDS model */
	overlap_markov_score = zGetRangeScore(fac->scan[EINIT_CDS],end_point,exon->start+5,frame); 
		
	if(exon_length <= 8) {
		overlap_markov_score = 0;
	}
		
	/* score rear here */
		
	rear_length = ((exon->end+1-4) - (exon->start+2+3*codon_num));
	einit_cds_endpt = exon->start+EINIT_CDS_LENGTH-1;

	if (einit_cds_endpt > (exon->end-3)){
		einit_cds_endpt = exon->end-3;
	}
	
	if (rear_length  > 0) {
                /* Scoring 5' end of Einit with EINIT_CDS model */
		rear_markov_score = zGetRangeScore(fac->scan[EINIT_CDS],einit_cds_endpt,exon->start+2+3*codon_num,frame); 
		
		if (einit_cds_endpt < (exon->end-3)){
			/* Scoring the rest of Einit with NONEINIT_CDS model */
			rear_markov_score += zGetRangeScore(fac->scan[NONEINIT_CDS],exon->end-3,einit_cds_endpt,frame); 
		}
	}
	else {
		rear_markov_score = 0;
	}
		
	/* score signal peptide here */
	signal_peptide.start = exon->start + 3;
	signal_peptide.end   = exon->start + 3 + 3*codon_num;
		
	signal_peptide.score = fac->scan[SIG_PEP]->scoref(fac->scan[SIG_PEP], 
                                                          &signal_peptide);
		
		
	overlap_coding_prob = 0.2*pow(2.0, signal_peptide.score/10.0) 
		+ 0.8*pow(2.0,overlap_markov_score/10.0);
	
	coding_score = 10*zLog2(overlap_coding_prob) + rear_markov_score; 
		
	if(exon_length > 5) {
		total_score += coding_score;
#ifdef DEBUG
                exon->content_score = coding_score;
#endif
	}
	else if(exon_length == 2) {
		/* your guess is as good as mine!  twinscan does it though */		  
		/* this corresponds to a probability of exactly 0.8 ?? */
		total_score += -3.21928;
#ifdef DEBUG
		exon->content_score = -3.21928;
#endif
	}
#ifdef DEBUG
	else {
		exon->content_score = 0;
	}
#endif

	return total_score;
}

static void zMakeEinitsFrom3(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	int roh, frame;
	coor_t i,begin,end;

	zSfeature* exon;
	if (pos < PADDING) return;
	if (zGetUScore(fac->scan[DON],pos+1) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {
		for (roh = 0; roh < 3; roh++) {
			frame = (pos+1 - roh) % 3;
			/* twinscan dosent allow initital exons of length less */
			/* than 3 but it does allow length 3 initial exons */
			end   = pos - roh - 2;
			begin = zGetStopSeqPrevPos(fac->fstop,end)+3;

			if(fac->estseq != NULL && !fac->est_para_mode) {
				if(zGetEstseqSeq(fac->estseq,pos+1) == '1'){
					return;
				}
				begin = find_prev_est(fac,pos,'2',begin-1)+1;
			}

			for (i = end; i >= begin; i -= 3) {
				if(zGetUScore(fac->scan[START],i) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
					continue;
				/*if(fac->estseq != NULL && !fac->est_para_mode) {
					flag = 1;
					flag = not_in_intron_check(fac, i, pos, fac->type);
					if(flag <  0) { continue; } 
				}*/
				exon = zSFListAppend(sfl);			
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name  = fac->type;
				exon->start = i;
				exon->end   = pos;
				exon->lfrag = 0;
				exon->rfrag = roh;
				exon->frame = frame;

				exon->score = fac->scoref(fac, exon);
			}
		}
	}
}


static void zMakeEinitsFrom5(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, begin, end;
	int frame;

	zSfeature* exon;
	if (pos < PADDING) return;

	if (zGetUScore(fac->scan[START],pos) > MINIMUM_TWINSCAN_SIGNAL_SCORE) {	  
		begin = pos + 2;
		end   = zGetStopSeqNextPos(fac->fstop,pos)+1;
		frame = pos % 3;
		/* twinscan dosent allow initital exons of length less */
		/* than 3 but it does allow length 3 initial exons */
		
		if(fac->estseq != NULL && !fac->est_para_mode) {
			end = find_next_est(fac,pos,'2',end+1)-1;
		}
		
		for (i = begin; i <= end; i++) {
			if(zGetUScore(fac->scan[DON],i+1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
				continue;
			/*if(fac->estseq != NULL && !fac->est_para_mode) {
				flag = 1;
				flag = not_in_intron_check(fac, pos, i, fac->type);
				if(flag <  0) continue;
			}*/
			if(fac->estseq != NULL && !fac->est_para_mode) {
				if(zGetEstseqSeq(fac->estseq,i+1) == '1'){
					continue;
				}
			}
	
			exon = zSFListAppend(sfl);			
			
			exon->group  = NULL;
			exon->strand = fac->strand;
			exon->name  = fac->type;
			exon->start = pos;
			exon->end   = i;
			exon->lfrag = 0;
			exon->rfrag = (i-pos+1)%3;
			exon->frame = frame;
			
			exon->score = fac->scoref(fac, exon);
		}
	}
}

static void zMakeEinitsFromInside(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, j, begin, end, min_end, max_end, min_begin;
	coor_t start[3];
	coor_t stop[3];
	int frame;

	zSfeature* exon;
	if (pos < PADDING) return;

	max_end = 0;
	min_begin = fac->dna->length;
	for(frame = 0;frame <= 2;frame++){
		stop[frame]  = zGetStopSeqNextPos(fac->fstop,pos-frame)+2;
		start[frame] = zGetStopSeqPrevPos(fac->fstop,stop[frame]-2-3)+3;
		if(max_end < stop[frame]){
			max_end = stop[frame];
		}
		if(min_begin > start[frame]){
			min_begin = start[frame];
		}
	}

	if(fac->estseq != NULL && !fac->est_para_mode) {
		min_begin = find_prev_est(fac,pos,'2',min_begin-1)+1;
		max_end = find_next_est(fac,pos,'2',max_end+1)-1;
	}
	for(frame = 0;frame <= 2;frame++){
		begin = start[frame];
		end = stop[frame];
		
		if(begin < min_begin){
			begin = min_begin;
		}
		if(end > max_end){
			end = max_end;
		}
		
		for(i = pos;i >= begin;i-=3){
			if(zGetUScore(fac->scan[START],i) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
				continue;
			min_end = pos;
			/* must be at lest 3bp long */
			if(min_end < i+2){
				min_end = i+2;
			}
			for(j = min_end; j < end; j++){
				if (zGetUScore(fac->scan[DON],j+1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
					continue;
				/*if(fac->estseq != NULL && !fac->est_para_mode) {
					flag = 1;
					flag = not_in_intron_check(fac, i, j, fac->type);
					if(flag <  0) continue;	
					}*/
				if(fac->estseq != NULL && !fac->est_para_mode) {
					if(zGetEstseqSeq(fac->estseq,j+1) == '1'){
						continue;
					} 
				}

				exon = zSFListAppend(sfl);			
	
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name  = fac->type;
				exon->start = i;
				exon->end   = j;
				exon->lfrag = 0;
				exon->rfrag = (exon->end - exon->start + 1) % 3;
				exon->frame = exon->start % 3;
				
				exon->score = fac->scoref(fac, exon);
			}
		}
	}
}

static score_t zScoreExon(const zFeatureFactory *fac, zSfeature *exon) {
	score_t total_score;
	score_t don_score, acc_score, cds_score;
	score_t doncon_score, acccon_score;
  	score_t donest_score, accest_score;
	score_t donalign_score, accalign_score, cdsalign_score;

	exon->cdscons_score = 0;	
	exon->cdsest_score  = 0;
	exon->align_score   = 0;	

	don_score = zGetUScore(fac->scan[DON],exon->end+1);
	acc_score = zGetUScore(fac->scan[ACC],exon->start-1);
	cds_score = zScoreCoding(fac, exon);

	total_score = don_score + acc_score + cds_score;

#ifdef DEBUG
	exon->begin_score   = acc_score;
	exon->end_score     = don_score;
	exon->content_score = cds_score;
#endif

	exon->cdscons_score = 0;
				
	if(fac->conseq_enabled) {
		/* this is an ineficent way to do this, come back and make it better later */
		doncon_score = zGetUScore(fac->conscan[DONCONS],exon->end+1) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->end-2,exon->end+6);
		
		/* this is an ineficent way to do this, come back and make it better later */
		acccon_score = zGetUScore(fac->conscan[ACCCONS],exon->start-1) - 
			zScoreRange(fac->conscan[INTRONCONS],exon->start-40,exon->start+2);
		
		/* this is an ineficent way to do this, come back and make it better later */
		if(exon->start+3 > exon->end-3){
			exon->cdscons_score = 0;
		}
		else{
			exon->cdscons_score = zScoreRange(fac->conscan[CDSCONS],exon->start+3,exon->end-3);
		}

		exon->cdscons_score += doncon_score + acccon_score;
		
		total_score += exon->cdscons_score;
	}

	exon->cdsest_score = 0;
				
	if(fac->est_para_mode) {
		
		donest_score = zGetUScore(fac->estscan[DONEST],exon->end+1) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->end-2,exon->end+6);
				  
		accest_score = zGetUScore(fac->estscan[ACCEST],exon->start-1) - 
			zScoreRange(fac->estscan[INTERGENICEST],exon->start-40,exon->start+2);

		if(exon->start+3 > exon->end-3){
			exon->cdsest_score = 0;
		}
		else{
			exon->cdsest_score = zScoreRange(fac->estscan[CDSEST],exon->start+3,exon->end-3);
		}
				  
		exon->cdsest_score += donest_score + accest_score;
		total_score += exon->cdsest_score;		
	}

	if (fac->phylo_enabled) {
		cdsalign_score = 0;
		if (fac->alignscan[CDSALIGN]->model->type == BNTREE) {
			/* no frame dependence */
			cdsalign_score = zScoreRangeA(fac->alignscan[CDSALIGN],exon->start+3,exon->end-3);
		}
		else if (fac->alignscan[CDSALIGN]->model->type == BNTREE_CDS) {
			/* frame dependence */
			cdsalign_score = fac->alignscan[CDSALIGN]->scoref(fac->alignscan[CDSALIGN], exon);
		}
		cdsalign_score -= zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start+3,exon->end-3);

		accalign_score = zGetUScoreA(fac->alignscan[ACCALIGN],exon->start-1) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->start-40,exon->start+2);
		
		donalign_score = zGetUScoreA(fac->alignscan[DONALIGN],exon->end+1) - 
			zScoreRangeA(fac->alignscan[INTRONALIGN],exon->end-2,exon->end+6);
		
		cdsalign_score *= PSC;
		accalign_score *= PSC;
		donalign_score *= PSC;

		exon->align_score = cdsalign_score + accalign_score + donalign_score;

		total_score += (cdsalign_score + accalign_score + donalign_score);
	}

	return total_score;
}

static void zMakeExonsFrom3(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, loh, roh, begin, end, min_start;
	coor_t start[3];
	int frame;
	zSfeature* exon;
	if (pos < PADDING) return;
	
	/* MINIMUM_TWINSCAN_SIGNAL_SCORE is for genscan compatibility */
	
	if (zGetUScore(fac->scan[DON],pos+1)  > MINIMUM_TWINSCAN_SIGNAL_SCORE) {	
		if(fac->estseq != NULL && !fac->est_para_mode){
			if(zGetEstseqSeq(fac->estseq,pos+1) == '1'){
				return;
			}
		}
		min_start = fac->dna->length;
		for (roh = 0; roh < 3; roh++) {
			start[roh] = zGetStopSeqPrevPos(fac->fstop,pos-2-roh)+1;
			if(min_start > start[roh]){
				min_start = start[roh];
			}
		}
		if(fac->estseq != NULL && !fac->est_para_mode){
			min_start = find_prev_est(fac,pos,'2',min_start-1)+1;
		}
		for (roh = 0; roh < 3; roh++) {
			frame = (pos+1 - roh) % 3;
			end   = pos-2;
			begin = start[roh];
			if(begin < min_start){
				begin = min_start;
			}
			/* this allows 3 bp internal exons. 1 and 2 bp exons restricted for 
			   backward compatability (although the old comments say 2 is allowed) */
			for (i = end; i >= begin; i--) { 
				if((i <= PADDING+43) || 
				   (zGetUScore(fac->scan[ACC],i-1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)) 
					continue;
				if(fac->estseq != NULL && !fac->est_para_mode){
					if(zGetEstseqSeq(fac->estseq,i-1) == '1'){
						continue;
					}
				}
				
				/*see ModelTrellis.cpp line 1153 for where the 43 comes from */
				/*if(fac->estseq != NULL && !fac->est_para_mode){
					flag = 1;
					flag = not_in_intron_check(fac, i, pos , fac->type);
					if(flag < 0) continue; 
					}*/
				loh = (pos - i + 1 - roh) % 3;
				
				exon = zSFListAppend(sfl);			
				
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name   = fac->type;
				exon->start  = i;
				exon->end    = pos;
				exon->lfrag  = loh;
				exon->rfrag  = roh;
				exon->frame  = frame;
				
				exon->score  = fac->scoref(fac, exon);
			}
		}
	}
}

static void zMakeExonsFrom5(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, loh, roh, begin, end, max_end;
	coor_t stop[3];
	int frame;
	zSfeature* exon;

	if (pos < PADDING) return;
	
	/* MINIMUM_TWINSCAN_SIGNAL_SCORE is for genscan compatibility */
	if ((pos > PADDING+43) && 
		(zGetUScore(fac->scan[ACC],pos-1) > MINIMUM_TWINSCAN_SIGNAL_SCORE)){	

		if(fac->estseq != NULL && !fac->est_para_mode){
			if(zGetEstseqSeq(fac->estseq,pos-1) == '1'){
				return;
			}
		}

		max_end = 0;
		for (loh = 0; loh < 3; loh++){		
			stop[loh] = zGetStopSeqNextPos(fac->fstop,pos+loh)+1;
			if(stop[loh] > max_end){
				max_end = stop[loh];
			}
		}
		if(fac->estseq != NULL && !fac->est_para_mode){
			max_end = find_next_est(fac,pos,'2',max_end+1)-1;
		}		
		for (loh = 0; loh < 3; loh++){		
			frame = (pos + loh) % 3;
			/* this allows 3bp exons but restricts 1 and 2bp exons */
			begin = pos+2;
			end   = stop[loh];
			if(end > max_end){
				end = max_end;
			}
			for (i = begin; i <= end; i++) { 
				if (zGetUScore(fac->scan[DON],i+1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
					continue;
				if(fac->estseq != NULL && !fac->est_para_mode){
					if(zGetEstseqSeq(fac->estseq,i+1) == '1'){
						continue;
					}
				}

				roh = (i - pos + 1 - loh) % 3; 

				exon = zSFListAppend(sfl);			

				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name   = fac->type;
				exon->start  = pos;
				exon->end    = i;
				exon->lfrag  = loh;
				exon->rfrag  = roh;
				exon->frame  = frame;
				
				exon->score  = fac->scoref(fac, exon);
			}
		}
	}
}

/* This is doing a bunch of extra work (searching at least some region crossing 
   pos 3 times, once for each frame) but the function should be used very rarely 
   so I'll fix it later EVAN*/
static void zMakeExonsFromInside(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t i, j, loh, begin, end;
	int frame;
	zSfeature* exon;

	if (pos < PADDING) return;
	
	for(frame = 0;frame <= 2;frame++){
		end   = zGetStopSeqNextPos(fac->fstop,pos-frame)+2;
		begin = zGetStopSeqPrevPos(fac->fstop,end-2-3)+3;
		if(begin < PADDING+43){
			begin = PADDING+43;
		}
		if(fac->estseq != NULL && !fac->est_para_mode){
			end = find_next_est(fac,pos,'2',end+1)-1;
			begin = find_prev_est(fac,pos,'2',begin-1)+1;
		}
		for(i = begin;i <= pos;i++){
			if(zGetUScore(fac->scan[ACC],i-1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
				continue;
			if(fac->estseq != NULL && !fac->est_para_mode){
				if(zGetEstseqSeq(fac->estseq,i-1) == '1'){
					continue;
				}				
			}
			loh = (pos-frame+3 - i) % 3; /* pos-frame is the first base of a codon
											+3 insures that pos-frame > i */
			for(j = zIntMax(pos,i+1); j < end; j++){
				if (zGetUScore(fac->scan[DON],j+1) <= MINIMUM_TWINSCAN_SIGNAL_SCORE)
					continue;
				if(fac->estseq != NULL && !fac->est_para_mode){
					if(zGetEstseqSeq(fac->estseq,j+1) == '1'){
						continue;
					}				
				}
				/*if(fac->estseq != NULL && !fac->est_para_mode){
					flag = 1;
					flag = not_in_intron_check(fac, i, pos , fac->type);
					if(flag < 0) continue; 
					}*/

				exon = zSFListAppend(sfl);			
				exon->group  = NULL;
				exon->strand = fac->strand;
				exon->name  = fac->type;
				exon->start = i;
				exon->end   = j;
				exon->lfrag = loh;
				exon->rfrag = (exon->end - exon->start + 1 - loh) %3;
				exon->frame = (exon->start + loh) % 3;
				
				exon->score = fac->scoref(fac, exon);
			}
		}
	}
}

static score_t zScorePolyA(const zFeatureFactory *fac, zSfeature *polya) {
	polya->cdscons_score = 0;	
	polya->cdsest_score  = 0;
	polya->align_score   = 0;	
	return zGetUScore(fac->scanner,polya->end - 2) - 5;
}

static void zMakePolyAsFrom3(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {

	zSfeature polya;
	int flag = -1; 
	
	if (pos < PADDING) return;

	polya.group  = NULL;
	polya.strand = '+';
	polya.name   = fac->type;
	polya.start  = pos - 5;
	polya.end    = pos;	
	polya.lfrag  = 0;
	polya.rfrag  = 0;
	polya.frame  = -1;
	polya.cdscons_score = 0;	
	polya.state = -1;
	polya.score = fac->scoref(fac, &polya);

#ifdef DEBUG
    polya.content_score = polya.score;
#endif
	
	if(fac->strand == '-' && pos == fac->dna->length - 5436 - 1){
		pos = 5436;
	}

	if(fac->estseq != NULL && !fac->est_para_mode) {
		flag = overlap_exon_or_intron(fac, polya.start, polya.end);
	}

	/*EVAN This bounds check is way off (Maybe correct for promoter?)
	  change it to something normal */

	if ((polya.score > MIN_POLYA_SCORE) && (pos > PADDING+42) 
		&& (pos < fac->scanner->seq->length-PADDING-43) && (flag < 0)) {
		zSFListInsert(sfl, &polya);		
#ifdef DEBUG
		polya.content_score = polya.score;
#endif
	}
}

static void zMakePolyAsFrom5(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	zMakePolyAsFrom3(fac,pos+5,sfl);
}

static void zMakePolyAsFromInside(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t end;
	for(end = pos;end <= pos+5;end++){
		zMakePolyAsFrom3(fac,end,sfl);
	}
}

static score_t zScoreProm(const zFeatureFactory *fac, zSfeature *prom) {
	score_t total;
	score_t best_tata = -10000;
	score_t cap_score;
	coor_t  i;
	
	prom->cdscons_score = 0;	
	prom->cdsest_score  = 0;
	prom->align_score   = 0;

	for(i = 28; i < 35; i++) {
		if (prom->end >= i && zGetUScore(fac->scan[TATA],prom->end-i) > best_tata)
			best_tata = zGetUScore(fac->scan[TATA],prom->end-i);
	}

	cap_score = fac->scan[PB_CAP]->score(fac->scan[PB_CAP], prom->end+1);
	
#ifdef DEBUG
	prom->begin_score = best_tata;
	prom->end_score = cap_score;
#endif
	if((cap_score > MINIMUM_TWINSCAN_SIGNAL_SCORE) && (best_tata > MINIMUM_TWINSCAN_SIGNAL_SCORE) && (prom->end > PADDING+81) 
	   && (prom->end < fac->scan[PB_CAP]->seq->length-PADDING-41)) {  
		total = 10*zLog2(0.7*pow(2,(best_tata+cap_score)/10.0)+0.3)-27.97;
	} else {
		total = MIN_SCORE;
	}

	return total;
}

static void zMakePromsFrom3(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	zSfeature prom;
	int flag = -1;

	if (pos < PADDING) return;
	
	prom.group  = NULL;
	prom.strand = fac->strand;
	prom.name   = fac->type;
	prom.start  = pos - 39;
	prom.end    = pos;
	prom.lfrag  = 0;
	prom.rfrag  = 0;
	prom.frame  = -1;
   
	/* Boundary condition. */
	/* Promoter feature is too close to start or end of sequence */
	if ((pos < 39) || (pos > fac->scan[PB_CAP]->seq->length-1)) return;
		
	if(fac->estseq != NULL && !fac->est_para_mode) {
		flag = overlap_exon_or_intron(fac, prom.start, prom.end);
	}
	
	if (((prom.score = fac->scoref(fac, &prom)) > MIN_SCORE) && (flag < 0)) {
		zSFListInsert(sfl, &prom);
	}  
	
}

static void zMakePromsFrom5(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	zMakePromsFrom3(fac,pos+39,sfl);
}

static void zMakePromsFromInside(const zFeatureFactory *fac, coor_t pos, zSFList* sfl) {
	coor_t end;
	for(end = pos;end <= pos+39;end++){
		zMakePromsFrom3(fac,end,sfl);
	}
}
/* !!! Really ought to break this into factories for each type of scanner */
void zInitEFactory(zFeatureFactory *factory, zDNA* dna, zEstseq* estseq, zScanner** scanners, zScanner**  conseqscanners, zScanner** estseqscanners, zAlignmentScanner **alignmentscanners, int feature_count, zStrIdx type, zFFCreate create5, zFFCreate create3, zFFCreate createi, strand_t strand, bool conseq_enabled, bool est_para_mode,bool phylo_enabled, zStopSeq* stop_seq, zEstseqSearch* estsearch){
	bool utr_struct_enabled; /* are we predicting 5' UTR structure? */

	/* scanners */
	zStrIdx   scanner;
	zScanner *einit_cds_scan, *noneinit_cds_scan;
	zScanner *accpt_scan, *donor_scan;
	zScanner *start_scan, *stop_scan;
	zScanner *sig_pep_scan;
	zScanner *ep_scan, *ea_scan;

 	zScanner *cdscons_scan;
 	zScanner *introncons_scan;
 	zScanner *initcons_scan;
 	zScanner *termcons_scan;
 	zScanner *acccons_scan;
 	zScanner *doncons_scan;
	zScanner *utrcons_scan;
 	zScanner *utracccons_scan;
 	zScanner *utrdoncons_scan;

 	zScanner *cdsest_scan;
 	zScanner *nullest_scan;
 	zScanner *initest_scan;
 	zScanner *termest_scan;
 	zScanner *accest_scan;
 	zScanner *donest_scan;
 	zScanner *epest_scan;
 	zScanner *encest_scan;
 	zScanner *eaest_scan;
 	zScanner *epaest_scan;
 	zScanner *utrdonest_scan;
 	zScanner *utraccest_scan;

	zAlignmentScanner *cdsalign_scan;
	zAlignmentScanner *intronalign_scan;
	zAlignmentScanner *initalign_scan;
	zAlignmentScanner *termalign_scan;
	zAlignmentScanner *accalign_scan;
	zAlignmentScanner *donalign_scan;
	zAlignmentScanner *ep_align_scan;
	zAlignmentScanner *enc_align_scan;
	zAlignmentScanner *ea_align_scan;
	zAlignmentScanner *epa_align_scan;
	zAlignmentScanner *utraccalign_scan;
	zAlignmentScanner *utrdonalign_scan;

	einit_cds_scan = NULL;
	noneinit_cds_scan = NULL;
	accpt_scan = NULL;
	donor_scan = NULL;
	start_scan = NULL;
	stop_scan = NULL;
	sig_pep_scan = NULL;
	ep_scan = NULL;
	ea_scan = NULL;

 	cdscons_scan = NULL;
 	introncons_scan = NULL;
 	initcons_scan = NULL;
 	termcons_scan = NULL;
 	acccons_scan = NULL;
 	doncons_scan = NULL;
	utrcons_scan  = NULL;
	utrdoncons_scan = NULL;
	utracccons_scan = NULL;	

 	cdsest_scan = NULL;
 	nullest_scan = NULL;
 	initest_scan = NULL;
 	termest_scan = NULL;
 	accest_scan = NULL;
 	donest_scan = NULL;
 	epest_scan = NULL;
 	encest_scan = NULL;
 	eaest_scan = NULL;
 	epaest_scan = NULL;
 	utrdonest_scan = NULL;
 	utraccest_scan = NULL;

	cdsalign_scan = NULL;
	intronalign_scan = NULL;
	initalign_scan = NULL;
	termalign_scan = NULL;
	accalign_scan = NULL;
	donalign_scan = NULL;
	ep_align_scan = NULL;
	enc_align_scan = NULL; 
	ea_align_scan = NULL;
	epa_align_scan = NULL;
	utraccalign_scan = NULL;
	utrdonalign_scan = NULL;

	/* setup scanners */
	scanner = zChar2StrIdx("EinitCoding");
	if (scanner >= feature_count ||
		(einit_cds_scan = scanners[scanner]) == NULL) {

                /* if there are no separate EinitCoding and   */
                /* NonEinitCoding models, look for an unsplit */
                /* Coding model                               */

                scanner = zChar2StrIdx("Coding");
                if ((scanner >= feature_count) || (scanners[scanner] == NULL)) {
		     zDie("Exons require a Coding model OR an EinitCoding AND a NonEinitCoding models");
                } else {
                     einit_cds_scan = scanners[scanner];
                     noneinit_cds_scan = scanners[scanner];
                }
        } else {
                scanner = zChar2StrIdx("NonEinitCoding");
                if (scanner >= feature_count ||
                        (noneinit_cds_scan = scanners[scanner]) == NULL) {
                        zDie("Exons require a NonEinitCoding model");
                }
        }

	scanner = zChar2StrIdx("Acceptor");
	if (scanner >= feature_count ||
		(accpt_scan = scanners[scanner]) == NULL) {
		zDie("Exons require a Acceptor model");
	}
	scanner = zChar2StrIdx("Donor");
	if (scanner >= feature_count ||
		(donor_scan = scanners[scanner]) == NULL) {
		zDie("Exons require a Donor model");
	}
	scanner = zChar2StrIdx("Start");
	if (scanner >= feature_count ||
		(start_scan = scanners[scanner]) == NULL) {
		zDie("Exons require a Start model");
	}
	scanner = zChar2StrIdx("Stop");
	if (scanner >= feature_count ||
		(stop_scan = scanners[scanner]) == NULL) {
		zDie("Exons require a Stop model");
	}
	scanner = zChar2StrIdx("SIG_PEP");
	if (scanner >= feature_count ||
		(sig_pep_scan = scanners[scanner]) == NULL) {
		zDie("Exons require a SIG_PEP model");
	}
	/* Presence of Ep model will determine whether UTR structure mode is enabled */
	
	scanner = zChar2StrIdx("Ep");
	if (scanner >= feature_count ||
		(ep_scan = scanners[scanner]) == NULL) {
		utr_struct_enabled = false;
	}
	else 
		utr_struct_enabled = true;
	
	if (utr_struct_enabled) {
		scanner = zChar2StrIdx("Ea");
		if (scanner >= feature_count ||
			(ea_scan = scanners[scanner]) == NULL) {
			zDie("Non-coding exons require an Ea model");
		}
	}
 

	if (NULL != conseqscanners) {
		scanner = zChar2StrIdx("CDSCONS");
		if (scanner >= feature_count || (cdscons_scan = conseqscanners[scanner]) == NULL) {
			zDie("Exons require a CDSCONS model");
		}
		scanner = zChar2StrIdx("INTRONCONS");
		if (scanner >= feature_count || (introncons_scan = conseqscanners[scanner]) == NULL) {
			zDie("Exons require a INTRONCONS model");
		}
		scanner = zChar2StrIdx("INITCONS");
		if (scanner >= feature_count || (initcons_scan = conseqscanners[scanner]) == NULL) {
			zDie("Exons require a INITCONS model");
		}
		scanner = zChar2StrIdx("TERMCONS");
		if (scanner >= feature_count || (termcons_scan = conseqscanners[scanner]) == NULL) {
			zDie("Exons require a TERMCONS model");
		}
		scanner = zChar2StrIdx("ACCCONS");
		if (scanner >= feature_count || (acccons_scan = conseqscanners[scanner]) == NULL) {
			zDie("Exons require a ACCCONS model");
		}
		scanner = zChar2StrIdx("DONCONS");
		if (scanner >= feature_count || (doncons_scan = conseqscanners[scanner]) == NULL) {
			zDie("Exons require a DONCONS model");
		}
		if (utr_struct_enabled) {
			scanner = zChar2StrIdx("UTRCONS");
			if (scanner >= feature_count || (utrcons_scan = conseqscanners[scanner]) == NULL) {
				zDie("Non-coding exons require a UTRCONS conseq model");
			}
			scanner = zChar2StrIdx("UTRACCCONS");
			if (scanner >= feature_count || (utracccons_scan = conseqscanners[scanner]) == NULL) {
				zDie("Exons require a UTRACCCONS model");
			}
			scanner = zChar2StrIdx("UTRDONCONS");
			if (scanner >= feature_count || (utrdoncons_scan = conseqscanners[scanner]) == NULL) {
				zDie("Exons require a UTRDONCONS model");
			}
		}
	}

	if (NULL != estseqscanners) {
		scanner = zChar2StrIdx("CDSEST");
		if (scanner >= feature_count || (cdsest_scan = estseqscanners[scanner]) == NULL) {
			zDie("Exons require a CDSEST model");
		}
		scanner = zChar2StrIdx("INTERGENICEST");
		if (scanner >= feature_count || (nullest_scan = estseqscanners[scanner]) == NULL) {
			zDie("Exons require a INTERGENICEST model");
		}
		scanner = zChar2StrIdx("INITEST");
		if (scanner >= feature_count || (initest_scan = estseqscanners[scanner]) == NULL) {
			zDie("Exons require a INITEST model");
		}
		scanner = zChar2StrIdx("TERMEST");
		if (scanner >= feature_count || (termest_scan = estseqscanners[scanner]) == NULL) {
			zDie("Exons require a TERMEST model");
		}
		scanner = zChar2StrIdx("ACCEST");
		if (scanner >= feature_count || (accest_scan = estseqscanners[scanner]) == NULL) {
			zDie("Exons require a ACCEST model");
		}
		scanner = zChar2StrIdx("DONEST");
		if (scanner >= feature_count || (donest_scan = estseqscanners[scanner]) == NULL) {
			zDie("Exons require a DONEST model");
		}
		if (utr_struct_enabled) {
			scanner = zChar2StrIdx("EPEST");
			if (scanner >= feature_count || (epest_scan = estseqscanners[scanner]) == NULL) {
				zDie("Exons require a EPEST model");
			}
			scanner = zChar2StrIdx("ENCEST");
			if (scanner >= feature_count || (encest_scan = estseqscanners[scanner]) == NULL) {
				zDie("Exons require a ENCEST model");
			}
			scanner = zChar2StrIdx("EAEST");
			if (scanner >= feature_count || (eaest_scan = estseqscanners[scanner]) == NULL) {
				zDie("Exons require a EAEST model");
			}
			scanner = zChar2StrIdx("EPAEST");
			if (scanner >= feature_count || (epaest_scan = estseqscanners[scanner]) == NULL) {
				zDie("Exons require a EPAEST model");
			}
			scanner = zChar2StrIdx("UTRDONEST");
			if (scanner >= feature_count || (utrdonest_scan = estseqscanners[scanner]) == NULL) {
				zDie("Exons require a UTRDONEST model");
			}
			scanner = zChar2StrIdx("UTRACCEST");
			if (scanner >= feature_count || (utraccest_scan = estseqscanners[scanner]) == NULL) {
				zDie("Exons require a UTRACCEST model");
			}
		}
	}

	if (NULL != alignmentscanners) {
		scanner = zChar2StrIdx("CDSCons");
		if (scanner >= feature_count || (cdsalign_scan = alignmentscanners[scanner]) == NULL) {
			zDie("Exons require a CDSCons model");
		}
		scanner = zChar2StrIdx("NCCons");
		if (scanner >= feature_count || (intronalign_scan = alignmentscanners[scanner]) == NULL) {
			zDie("Exons require a NCCons model");
		}
		scanner = zChar2StrIdx("StartCodonCons");
		if (scanner >= feature_count || (initalign_scan = alignmentscanners[scanner]) == NULL) {
			zDie("Exons require a StartCodonCons model");
		}
		scanner = zChar2StrIdx("StopCodonCons");
		if (scanner >= feature_count || (termalign_scan = alignmentscanners[scanner]) == NULL) {
			zDie("Exons require a StopCodonCons model");
		}
		scanner = zChar2StrIdx("AcceptorCons");
		if (scanner >= feature_count || (accalign_scan = alignmentscanners[scanner]) == NULL) {
			zDie("Exons require an AcceptorCons model");
		}
		scanner = zChar2StrIdx("DonorCons");
		if (scanner >= feature_count || (donalign_scan = alignmentscanners[scanner]) == NULL) {
			zDie("Exons require a DonorCons model");
		}
		if (utr_struct_enabled) {
			scanner = zChar2StrIdx("EpCons");
			if (scanner >= feature_count || (ep_align_scan = alignmentscanners[scanner]) == NULL) {
				zDie("Exons require an EpCons model");
			}
			scanner = zChar2StrIdx("EncCons");
			if (scanner >= feature_count || (enc_align_scan = alignmentscanners[scanner]) == NULL) {
				zDie("Exons require an EncCons model");
			}
			scanner = zChar2StrIdx("EaCons");
			if (scanner >= feature_count || (ea_align_scan = alignmentscanners[scanner]) == NULL) {
				zDie("Exons require an EaCons model");
			}
			scanner = zChar2StrIdx("EpaCons");
			if (scanner >= feature_count || (epa_align_scan = alignmentscanners[scanner]) == NULL) {
				zDie("Exons require an EpaCons model");
			}
			scanner = zChar2StrIdx("UtrAcceptorCons");
			if (scanner >= feature_count || (utraccalign_scan = alignmentscanners[scanner]) == NULL) {
				zDie("Exons require a UtrAcceptorCons model");
			}
			scanner = zChar2StrIdx("UtrDonorCons");
			if (scanner >= feature_count || (utrdonalign_scan = alignmentscanners[scanner]) == NULL) {
				zDie("Exons require a UtrDonorCons model");
			}
		}
	}

	/* Set Coding to a single Isochore -------------------------------------- */
	if (ISO == einit_cds_scan->model[0].type) {
		einit_cds_scan = zGetScannerForIso(einit_cds_scan, zGetDNAGC(dna));
	}
	if (ISO == noneinit_cds_scan->model[0].type) {
		noneinit_cds_scan = zGetScannerForIso(noneinit_cds_scan, zGetDNAGC(dna));
	}

	/* Precompute scores ---------------------------------------------------- */
	zPreComputeScanner(accpt_scan);
	zPreComputeScanner(donor_scan);
	zPreComputeScanner(start_scan);
	zPreComputeScanner(stop_scan);	
	zPreComputeScanner(einit_cds_scan);
	zPreComputeScanner(noneinit_cds_scan);
	if (utr_struct_enabled) {
		zPreComputeScanner(ep_scan);
		zPreComputeScanner(ea_scan);
	}

	if(conseq_enabled) {
		zPreComputeScanner(cdscons_scan);
		zPreComputeScanner(introncons_scan);
		zPreComputeScanner(initcons_scan);
		zPreComputeScanner(termcons_scan);
		zPreComputeScanner(acccons_scan);
		zPreComputeScanner(doncons_scan);
		if (utr_struct_enabled) {
			zPreComputeScanner(utrcons_scan);
			zPreComputeScanner(utracccons_scan);
			zPreComputeScanner(utrdoncons_scan);
		}
	}

	if(est_para_mode) {
		zPreComputeScanner(cdsest_scan);
		zPreComputeScanner(nullest_scan);
		zPreComputeScanner(initest_scan);
		zPreComputeScanner(termest_scan);
		zPreComputeScanner(accest_scan);
		zPreComputeScanner(donest_scan);
		if (utr_struct_enabled) {
			zPreComputeScanner(epest_scan);
			zPreComputeScanner(encest_scan);
			zPreComputeScanner(eaest_scan);
			zPreComputeScanner(epaest_scan);
			zPreComputeScanner(utraccest_scan);
			zPreComputeScanner(utrdonest_scan);
		}
	}

	if(phylo_enabled) {
		zPreComputeAlignmentScanner(cdsalign_scan);
		zPreComputeAlignmentScanner(intronalign_scan);
		zPreComputeAlignmentScanner(initalign_scan);
		zPreComputeAlignmentScanner(termalign_scan);
		zPreComputeAlignmentScanner(accalign_scan);
		zPreComputeAlignmentScanner(donalign_scan);
		if (utr_struct_enabled) {
			zPreComputeAlignmentScanner(ep_align_scan);
			zPreComputeAlignmentScanner(enc_align_scan);
			zPreComputeAlignmentScanner(ea_align_scan);
			zPreComputeAlignmentScanner(epa_align_scan);
			zPreComputeAlignmentScanner(utraccalign_scan);
			zPreComputeAlignmentScanner(utrdonalign_scan);
		}
	}

	/* set factory attributes ----------------------------------------------- */
	factory->fstop   = stop_seq;
	factory->estsearch   = estsearch;
	factory->create5 = create5;
	factory->create3 = create3;
	factory->createi = createi;
	factory->type    = type;
	factory->strand  = strand;
	factory->offset  = einit_cds_scan->subscanner[0].model->focus;
	
	factory->scan    = zMalloc(NUM_EXON_SCANNERS * sizeof(zScanner*),
							   "zInitEfactory scanner table");

	factory->conscan    = zMalloc(NUM_EXON_CONSCANNERS * sizeof(zScanner*),
								  "zInitEfactory conseqscanner table");

	factory->conseq_enabled = conseq_enabled;

	factory->estscan    = zMalloc(NUM_EXON_ESTSCANNERS * sizeof(zScanner*),
								  "zInitEfactory estseqscanner table");

	factory->est_para_mode = est_para_mode;

	factory->alignscan = zMalloc (NUM_EXON_ALIGNSCANNERS * sizeof (zAlignmentScanner *),
								  "zInitEfactory alignscanner table");
	factory->phylo_enabled = phylo_enabled;

	factory->scan[EINIT_CDS]       = einit_cds_scan;
	factory->scan[NONEINIT_CDS]    = noneinit_cds_scan;
	factory->scan[START]           = start_scan;
	factory->scan[STOP]            = stop_scan;
	factory->scan[ACC]             = accpt_scan;
	factory->scan[DON]             = donor_scan;
	factory->scan[SIG_PEP]         = sig_pep_scan;
	factory->scan[EA]              = ea_scan;
	factory->scan[EP]              = ep_scan;

 	factory->conscan[CDSCONS]      = cdscons_scan;
	factory->conscan[INTRONCONS]   = introncons_scan;
	factory->conscan[INITCONS]     = initcons_scan;
	factory->conscan[TERMCONS]     = termcons_scan;
	factory->conscan[ACCCONS]      = acccons_scan;
	factory->conscan[DONCONS]      = doncons_scan;
        factory->conscan[UTRCONS]      = utrcons_scan;
        factory->conscan[UTRDONCONS]   = utrdoncons_scan;
        factory->conscan[UTRACCCONS]   = utracccons_scan;

	factory->estseq = estseq;

 	factory->estscan[CDSEST]      = cdsest_scan;
	factory->estscan[INTERGENICEST]   = nullest_scan;
	factory->estscan[INITEST]     = initest_scan;
	factory->estscan[TERMEST]     = termest_scan;
	factory->estscan[ACCEST]      = accest_scan;
	factory->estscan[DONEST]      = donest_scan;
	factory->estscan[EPEST]       = epest_scan;
	factory->estscan[ENCEST]      = encest_scan;
	factory->estscan[EAEST]       = eaest_scan;
	factory->estscan[EPAEST]      = epaest_scan;
	factory->estscan[UTRACCEST]   = utraccest_scan;
	factory->estscan[UTRDONEST]   = utrdonest_scan;

 	factory->alignscan[CDSALIGN]      = cdsalign_scan;
	factory->alignscan[INTRONALIGN]   = intronalign_scan;
	factory->alignscan[INITALIGN]     = initalign_scan;
	factory->alignscan[TERMALIGN]     = termalign_scan;
	factory->alignscan[ACCALIGN]      = accalign_scan;
	factory->alignscan[DONALIGN]      = donalign_scan;
	factory->alignscan[EP_ALIGN]      = ep_align_scan;
	factory->alignscan[ENC_ALIGN]     = enc_align_scan;
	factory->alignscan[EA_ALIGN]      = ea_align_scan;
	factory->alignscan[EPA_ALIGN]     = epa_align_scan;
	factory->alignscan[UTRACCALIGN]   = utraccalign_scan;
	factory->alignscan[UTRDONALIGN]   = utrdonalign_scan;

	/* this stuff isn't used by an EFactory */
	factory->length  = 0;
	factory->scanner = NULL;
	factory->dna     = dna;
	factory->min_len = 0;
	factory->score   = 0;
}

static void zInitEpaFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq, 
							zScanner **scanners, 
							zScanner **conseqscanners,
							zScanner **estseqscanners,
							zAlignmentScanner **alignmentscanners,
							int feature_count, strand_t strand, 
							bool conseq_enabled, bool est_para_mode,
							bool phylo_enabled, zStopSeq* stop_seq,
							zEstseqSearch* estsearch) {

	fac->end_offset = 0;
	fac->scoref = zScoreEpa;
	zInitEFactory(fac, dna, estseq, scanners, conseqscanners, estseqscanners, 
				  alignmentscanners, feature_count, zChar2StrIdx("Epa"), 
				  zMakeEpasFrom5, zMakeEpasFrom3, zMakeEpasFromInside,
				  strand, conseq_enabled, est_para_mode, phylo_enabled, stop_seq, estsearch);
}

static void zInitEpFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq, 
						   zScanner **scanners, 
						   zScanner **conseqscanners,
						   zScanner **estseqscanners,
						   zAlignmentScanner **alignmentscanners,
						   int feature_count, strand_t strand, 
						   bool conseq_enabled, bool est_para_mode,
						   bool phylo_enabled, zStopSeq* stop_seq,
						   zEstseqSearch* estsearch) {
	
	fac->end_offset = 0;
	fac->scoref = zScoreEp;
	zInitEFactory(fac, dna, estseq, scanners, conseqscanners, estseqscanners, 
				  alignmentscanners, feature_count, zChar2StrIdx("Ep"), 
				  zMakeEpsFrom5, zMakeEpsFrom3, zMakeEpsFromInside, 
				  strand, conseq_enabled, est_para_mode, phylo_enabled, stop_seq, estsearch);
}

static void zInitEncFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq, 
							zScanner **scanners, 
							zScanner **conseqscanners,
							zScanner **estseqscanners,
							zAlignmentScanner **alignmentscanners,
							int feature_count, strand_t strand, 
							bool conseq_enabled, bool est_para_mode,
							bool phylo_enabled, zStopSeq* stop_seq,
							zEstseqSearch* estsearch) {
	
	fac->end_offset = 0;
	fac->scoref = zScoreEnc;
	zInitEFactory(fac, dna, estseq, scanners, conseqscanners, estseqscanners, 
				  alignmentscanners, feature_count, zChar2StrIdx("Enc"), 
				  zMakeEncsFrom5, zMakeEncsFrom3, zMakeEncsFromInside, 
				  strand, conseq_enabled, est_para_mode, phylo_enabled, stop_seq, estsearch);
}

static void zInitEaFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq, 
						   zScanner **scanners, 
						   zScanner **conseqscanners,
						   zScanner **estseqscanners,
						   zAlignmentScanner **alignmentscanners,
						   int feature_count, strand_t strand, 
						   bool conseq_enabled, bool est_para_mode,
						   bool phylo_enabled, zStopSeq* stop_seq,
						   zEstseqSearch* estsearch) {


	fac->end_offset = 0;
	fac->scoref = zScoreEa;
	zInitEFactory(fac, dna, estseq, scanners, conseqscanners, estseqscanners, 
				  alignmentscanners, feature_count, zChar2StrIdx("Ea"), 
				  zMakeEasFrom5,zMakeEasFrom3,zMakeEasFromInside, 
				  strand, conseq_enabled, est_para_mode, phylo_enabled, stop_seq, estsearch);
}

static void zInitExonFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq, 
							 zScanner **scanners, 
							 zScanner **conseqscanners, 
							 zScanner **estseqscanners,
							 zAlignmentScanner **alignmentscanners,
							 int feature_count, strand_t strand, 
							 bool conseq_enabled, bool est_para_mode,
							 bool phylo_enabled, zStopSeq* stop_seq,
							 zEstseqSearch* estsearch) {
	
	fac->end_offset = 0;
	fac->scoref = zScoreExon;
	zInitEFactory(fac, dna, estseq, scanners, conseqscanners, estseqscanners, 
				  alignmentscanners,feature_count, zChar2StrIdx("Exon"), 
				  zMakeExonsFrom5, zMakeExonsFrom3, zMakeExonsFromInside, 
				  strand, conseq_enabled, est_para_mode, phylo_enabled, stop_seq, estsearch);
}

static void zInitEinitFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq,
							  zScanner **scanners, 
							  zScanner **conseqscanners, 
							  zScanner **estseqscanners,
							  zAlignmentScanner **alignmentscanners,
							  int feature_count, strand_t strand, 
							  bool conseq_enabled, bool est_para_mode,
							  bool phylo_enabled, zStopSeq* stop_seq,
							  zEstseqSearch* estsearch) {

	fac->end_offset = 0;
	fac->scoref = zScoreEinit;
	zInitEFactory(fac, dna, estseq, scanners, conseqscanners, estseqscanners, 
				  alignmentscanners,feature_count, zChar2StrIdx("Einit"), 
				  zMakeEinitsFrom5, zMakeEinitsFrom3, zMakeEinitsFromInside, 
				  strand, conseq_enabled, est_para_mode, phylo_enabled, stop_seq, estsearch);
}

static void zInitEtermFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq,
							  zScanner **scanners,
							  zScanner **conseqscanners, 
							  zScanner **estseqscanners,
							  zAlignmentScanner **alignmentscanners,
							  int feature_count, strand_t strand, 
							  bool conseq_enabled, bool est_para_mode,
							  bool phylo_enabled, zStopSeq* stop_seq,
							  zEstseqSearch* estsearch) {

	fac->end_offset = 0;
	fac->scoref = zScoreEterm;
	zInitEFactory(fac, dna, estseq, scanners, conseqscanners, estseqscanners, 
				  alignmentscanners,feature_count, zChar2StrIdx("Eterm"), 
				  zMakeEtermsFrom5, zMakeEtermsFrom3, zMakeEtermsFromInside, 
				  strand, conseq_enabled, est_para_mode, phylo_enabled, stop_seq, estsearch);
}

static void zInitEsnglFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq,
							  zScanner **scanners, 
							  zScanner **conseqscanners, 
							  zScanner **estseqscanners,
							  zAlignmentScanner **alignmentscanners,
							  int feature_count, strand_t strand, 
							  bool conseq_enabled, bool est_para_mode,
							  bool phylo_enabled, zStopSeq* stop_seq,
							  zEstseqSearch* estsearch) {

	fac->end_offset = 0;
	fac->scoref = zScoreEsngl;
	zInitEFactory(fac, dna, estseq, scanners, conseqscanners, estseqscanners, 
				  alignmentscanners,feature_count, zChar2StrIdx("Esngl"), 
				  zMakeEsnglsFrom5, zMakeEsnglsFrom3, zMakeEsnglsFromInside, 
				  strand, conseq_enabled, est_para_mode, phylo_enabled, stop_seq, estsearch);
}

static void zInitPolyAFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq,
							  zScanner **scanners, 
							  zScanner **conseqscanners, 
							  zScanner **estseqscanners,
							  zAlignmentScanner **alignmentscanners,
							  int feature_count, strand_t strand, 
							  bool conseq_enabled, bool est_para_mode,
							  bool phylo_enabled, zStopSeq* stop_seq,
							  zEstseqSearch* estsearch) {
	
	stop_seq = stop_seq; /* compiler hush */	
	fac->end_offset = 0;
	fac->create5 = zMakePolyAsFrom5;
	fac->create3 = zMakePolyAsFrom3;
	fac->createi = zMakePolyAsFromInside;
	fac->scoref  = zScorePolyA;
	fac->type    = zChar2StrIdx("PolyA");
	fac->strand  = strand;
	
	/* Unused by promoters ------------------------- */
	fac->length  = 0;
	fac->scanner = NULL;

	fac->dna     = dna;
	fac->min_len = 0;
	fac->score   = MIN_SCORE;

	fac->offset  = 0;
	fac->fstop   = NULL;
	fac->estsearch = estsearch;

	fac->scan    = NULL;
	fac->conscan = conseqscanners = NULL; /* compiler hush */
	fac->conseq_enabled = conseq_enabled;

	fac->estseq = estseq;

	fac->estscan = estseqscanners = NULL; /* compiler hush */
	fac->est_para_mode = est_para_mode;

	fac->alignscan = alignmentscanners = NULL; /* compiler hush */
	fac->phylo_enabled = phylo_enabled;

	/* --------------------------------------------- */

	if (fac->type >= feature_count || scanners[fac->type] == NULL) {
		zDie("PolyA's require a PolyA model");
	}
	fac->scanner = scanners[fac->type];
	zPreComputeScanner(fac->scanner);

}

static void zInitPromFactory(zFeatureFactory *fac, zDNA *dna, zEstseq *estseq,
							 zScanner **scanners, 
							 zScanner **conseqscanners, 
							 zScanner **estseqscanners,
                             zAlignmentScanner **alignmentscanners,
							 int feature_count, 
							 strand_t strand, 
							 bool conseq_enabled, bool est_para_mode,
							 bool phylo_enabled, zStopSeq* stop_seq,
							 zEstseqSearch* estsearch) {

	zStrIdx tmp_name;

	stop_seq = stop_seq; /* compiler hush */	
	fac->end_offset = 0;
	fac->create5 = zMakePromsFrom5;
	fac->create3 = zMakePromsFrom3;
	fac->createi = zMakePromsFromInside;
	fac->scoref  = zScoreProm;
	fac->type    = zChar2StrIdx("Prom");
	fac->strand  = strand;
	
	/* Unused by promoters ------------------------- */
	fac->length  = 0;
	fac->scanner = NULL;

	fac->dna     = dna;
	fac->min_len = 0;
	fac->score   = MIN_SCORE;

	fac->offset  = 0;
	fac->fstop   = NULL;
	fac->estsearch = estsearch;
	/* --------------------------------------------- */


	fac->scan    = zCalloc(NUM_PROM_SCANNERS, sizeof(zScanner*),
						   "zInitPromFactory: scanner array");
	fac->conscan = conseqscanners = NULL; /* compiler hush */
	fac->conseq_enabled = conseq_enabled;

	fac->estseq = estseq;
	fac->estscan = estseqscanners = NULL; /* compiler hush */
	fac->est_para_mode = est_para_mode;

	fac->alignscan = alignmentscanners = NULL; /* compiler hush */
	fac->phylo_enabled = phylo_enabled;

	tmp_name = zChar2StrIdx("TATA");
	if (tmp_name >= feature_count || scanners[tmp_name] == NULL) {
		zDie("Promoters require a TATA model");
	}
	fac->scan[TATA] = scanners[tmp_name];

	tmp_name = zChar2StrIdx("PB_CAP");
	if (tmp_name >= feature_count || scanners[tmp_name] == NULL) {
		zDie("Promoters require a PB_CAP model");
	}
	fac->scan[PB_CAP] = scanners[tmp_name];
	   
	zPreComputeScanner(fac->scan[TATA]);
}

/*
  void zInitSFactory(zFeatureFactory *factory, zScanner *scanner, coor_t length,
  zSfeatureName type)
  {	
  factory->create  = zMakeSfeatures;
  factory->type    = type;
  factory->length  = length;
  factory->scanner = scanner;
	
  / not used by SFactory /
  factory->scanner = NULL;
  factory->dna     = NULL;
  factory->min_len = 0;
  factory->score   = 0;
  factory->Einit   = -1;
  factory->Exon    = -1;
  factory->Eterm   = -1;
  factory->Esngl   = -1;
  factory->cds[0]  = NULL;
  factory->cds[1]  = NULL;
  factory->cds[2]  = NULL;
  factory->acc     = NULL;
  factory->don     = NULL;
  factory->start   = NULL;
  factory->stop    = NULL;
  factory->fstop = NULL;
  }

  void zInitRFactory(zFeatureFactory *factory, zDNA *dna, score_t score,
  coor_t min_len)
  {	
  factory->create  = zMakeRepeats;
  factory->type    = Repeat;
  factory->dna     = dna;
  factory->min_len = min_len;
  factory->score   = score;

  / not used by RFactory /
  factory->length  = 0;
  factory->scanner = 0;
  factory->Einit   = -1;
  factory->Exon    = -1;
  factory->Eterm   = -1;
  factory->Esngl   = -1;
  factory->cds[0]  = NULL;
  factory->cds[1]  = NULL;
  factory->cds[2]  = NULL;
  factory->acc     = NULL;
  factory->don     = NULL;
  factory->start   = NULL;
  factory->stop    = NULL;
  factory->fstop   = NULL;
  }
*/

void zFreeFeatureFactory(zFeatureFactory *f) {
	
	/* fstop is shared among multiple factories and freed by the trellis */
	f->fstop   = NULL;
	zFree(f->scan);    f->scan    = NULL;
	zFree(f->conscan); f->conscan = NULL;
	zFree(f->estscan); f->estscan = NULL;
	zFree(f->alignscan); f->alignscan = NULL;
}

/******************************************************************************\
 The Factory Registry
\******************************************************************************/

zFFInit zGetFactory(const char* feature_name) {

	     if (strcmp(feature_name, "Epa")   == 0) return zInitEpaFactory;
	else if (strcmp(feature_name, "Ep")    == 0) return zInitEpFactory;
	else if (strcmp(feature_name, "Ea")    == 0) return zInitEaFactory;
	else if (strcmp(feature_name, "Enc")   == 0) return zInitEncFactory;
	else if (strcmp(feature_name, "Exon")  == 0) return zInitExonFactory;
	else if (strcmp(feature_name, "Einit") == 0) return zInitEinitFactory;
	else if (strcmp(feature_name, "Eterm") == 0) return zInitEtermFactory;
	else if (strcmp(feature_name, "Esngl") == 0) return zInitEsnglFactory;
	else if (strcmp(feature_name, "PolyA") == 0) return zInitPolyAFactory;
	else if (strcmp(feature_name, "Prom")  == 0) return zInitPromFactory;

	return NULL;
}

zScanner* zGetScanner( const zFeatureFactory *fac, int type) {
	float gc;
	int i;
	zScanner *scanner = fac->scan[type];
	gc = zGetDNAGC(fac->dna);

	if(scanner == NULL) {
		zDie("Invalid scanner type in zGetScanner()");
	}
	if ( scanner->model->type == ISO ) {
		i = 0;
		while (gc > scanner->model->data[i] && i < scanner->model->submodels-1) i++;
		return &scanner->subscanner[i];
	} else {
		return scanner;
	}
}

#endif
