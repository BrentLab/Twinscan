/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
  zScanner.c - part of the ZOE library for genomic analysis
																			   
  Copyright (C) 2001-2002 Ian F. Korf
  
\******************************************************************************/

#ifndef ZOE_SCANNER_C
#define ZOE_SCANNER_C

#include "zScanner.h"

#define CDS_BLOCK_SIZE 50000;
#define NONCDS_BLOCK_SIZE 5000;
#define USCORE_BLOCK_COUNT 2;

static int zGetScannerVariantScoreIndex(zScannerVariant*);
static void zFillCDSUScore(zScanner*,coor_t,coor_t,score_t*,score_t*,score_t*);
static score_t* zFillUScoreVariant(zScanner*,int,int,int,int);
static void zFillUScore(zScanner*,int,int);

static char seq2sig (int c) {
	switch (c) {
	case 'A': case 'a':	return  8; /* 1000 */
	case 'C': case 'c':	return  4; /* 0100 */
	case 'G': case 'g':	return  2; /* 0010 */
	case 'T': case 't':	return  1; /* 0001 */
	case 'R': case 'r':	return 10; /* 1010 */
	case 'Y': case 'y':	return  5; /* 0101 */
	case 'M': case 'm':	return 12; /* 1100 */
	case 'K': case 'k':	return  3; /* 0011 */
	case 'W': case 'w':	return  9; /* 1001 */
	case 'S': case 's':	return  6; /* 0110 */
	case 'B': case 'b':	return  7; /* 0111 */
	case 'D': case 'd':	return 11; /* 1011 */
	case 'H': case 'h':	return 13; /* 1101 */
	case 'V': case 'v':	return 14; /* 1110 */
	case 'N': case 'n':	return 15; /* 1111 */
	default:
		zWarn("illegal symbols seq2sig");
		return 15;
	}
}

static void zFillUScore(zScanner*,int,int);

/*******\
  SIG
\*******/
static score_t zDNAScoreSIGRange (zScanner *scanner, coor_t start, coor_t end){

	int sequence_size, codons, j;
	int index=0;
	score_t score=0;
	score_t codon_score;

	sequence_size = end - start;
	codons = sequence_size / 3;

	for(j=0;j<codons;j++) {

		index = 0;

		if(j >= 5 ) {
			index += 125;
		}

		
		index += 25*zGetDNAS5(scanner->seq->parent,start + 3*j + 0)
			+ 5*zGetDNAS5(scanner->seq->parent,start + 3*j + 1)
			+ 1*zGetDNAS5(scanner->seq->parent,start + 3*j + 2);
		
		codon_score = scanner->model->data[index];

		score += codon_score;
		
	}

	return score;
}

static score_t zDNAScoreSIG (zScanner *scanner, zSfeature *signal_peptide) {
	return zDNAScoreSIGRange(scanner,signal_peptide->start,signal_peptide->end);
}

/*******\
   WMM
\*******/
static score_t zDNAScoreWMM (zScanner *scanner, coor_t pos) {
	coor_t i;
	score_t score;
	coor_t mfocus = pos - scanner->model->focus;
	int num;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos))
		return MIN_SCORE;

	/* scoring */
	score = 0;
	for (i = 0; i < scanner->model->length; i++) {
		num = (i * scanner->model->symbols) + zGetDNAS5(scanner->seq->parent,i + mfocus);

		score += scanner->model->data[num];
	}
	return score;
}

/*******\
   LUT
\*******/
static score_t zDNAScoreLUT (zScanner *scanner, coor_t pos) {
	coor_t i, p, index;
	coor_t mfocus = pos - scanner->model->focus;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos))
		return MIN_SCORE;
	
	/* get the precomputed score if available */

	/* scoring */
	index = 0;
	for (i = 0; i < scanner->model->length; i++) {
		p = zPOWER[scanner->model->symbols][scanner->model->length -i -1];
		index += (p * zGetDNAS5(scanner->seq->parent,i + mfocus));
	}

	return scanner->model->data[index];	
}

/************\
   Pair LUT
\************/
static score_t zDNAScorePairLUT (zScanner *scanner, zDNA* genomic, zDNA* cdna, coor_t gpos, coor_t cpos) {
	coor_t i, p, index, length;
	coor_t gfocus = gpos - scanner->model->focus;
	coor_t cfocus = cpos - scanner->model->focus;

	/* user defines and boundaries */
	if ((gpos < scanner->min_gpos) || (gpos > scanner->max_pos) ||
	    (cpos < scanner->min_pos) || (cpos > scanner->max_pos))
		return MIN_SCORE;
	
	/* scoring */
	if (scanner->model->length == 2) { /* Bypass for the basic version where we dont need the complex loops that take time */
		index = scanner->model->symbols*zGetDNAS5(genomic, gfocus) + zGetDNAS5(cdna, cfocus);
		return scanner->model->data[index];
	}

	length = scanner->model->length / 2;
	/* scoring */
	index = 0;

	/* Conditional probs */
	for (i = 0; i < length; i++) {
		p = zPOWER[scanner->model->symbols][i];
		index += (p * zGetDNAS5(genomic, gfocus - i));
	}
	index *= zPOWER[scanner->model->symbols][length];

	for (i = 0; i < length; i++) {
		p = zPOWER[scanner->model->symbols][i];
		index += (p * zGetDNAS5(cdna, cfocus - i));
	}

	return scanner->model->data[index];
}

static score_t zConseqScoreLUT (zScanner *scanner, coor_t pos) {
	coor_t i, p, index;
	coor_t mfocus = pos - scanner->model->focus;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos))
		return MIN_SCORE;

	/* scoring */
	index = 0;
	for (i = 0; i < scanner->model->length; i++) {
		p = zPOWER[scanner->model->symbols][scanner->model->length -i -1];
		index += (p * zGetConseqS10(scanner->seq->parent,i + mfocus));
	}

	return scanner->model->data[index];	
}

static score_t zEstseqScoreLUT (zScanner *scanner, coor_t pos) {
	coor_t i, p, index;
	coor_t mfocus = pos - scanner->model->focus;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos))
		return MIN_SCORE;
	
	/* scoring */
	index = 0;
	for (i = 0; i < scanner->model->length; i++) {
		p = zPOWER[scanner->model->symbols][scanner->model->length -i -1];
		index += (p * zGetEstseqS10(scanner->seq->parent,i + mfocus));
	}
	
	return scanner->model->data[index];	
}


/********\
   WWAM
\********/
/* EVAN modified to use context outside of scored region, previously ignored last two positions, need to verify that parameter file are still correct */
static score_t zDNAScoreWWAM (zScanner *scanner, coor_t pos) {
	coor_t  i, p, index, offset;
	int     j;
	coor_t  mfocus = pos - scanner->model->focus;
	score_t score = 0;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) {
		return MIN_SCORE;
	}
	score = 0;

	/* scoring */
	for (i = 0; i < scanner->model->length; i++) {
		index = 0;
		for (j = 0; j < scanner->model->order; j++) {
			p = zPOWER[scanner->model->symbols][scanner->model->order - j - 1];
			/* EVAN changed this back to orig for compatability, it should be fixed again 
			   index += (p * zGetDNAS5(scanner->seq->parent,i + mfocus - 
			   (scanner->model->order - j - 1)));
			*/
			index += (p * zGetDNAS5(scanner->seq->parent,i + mfocus +j));
		}
		
		offset = i*zPOWER[scanner->model->symbols][scanner->model->order];
		score += scanner->model->data[index+offset];
	}
	
	return score;
}

static score_t zConseqScoreWWAM (zScanner *scanner, coor_t pos) {
	coor_t  i, p, index, offset;
	int j;
	coor_t  mfocus = pos - scanner->model->focus;
	score_t score = 0;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) {
		return MIN_SCORE;
	}
	score = 0;

	/* scoring */
	for (i = 0; i < scanner->model->length; i++) {
		index = 0;
		for (j = 0; j < scanner->model->order; j++) {
			p = zPOWER[scanner->model->symbols][scanner->model->order - j -1];
			/*EVAN changed this back for compatability, needs to change back 
			index += (p * zGetConseqS10(scanner->seq->parent,i + j +
										mfocus - scanner->model->order + 1));
			*/
			index += (p * zGetConseqS10(scanner->seq->parent,i + j + mfocus));
		}

		offset = i*zPOWER[scanner->model->symbols][scanner->model->order];
		score += scanner->model->data[index+offset];
	}

	return score;
}

static score_t zConseqScoreWWAMRange (zScanner *scanner, coor_t start, coor_t end){

	if(end-start+1 != scanner->model->length){
		zDie("range must have same size as model length in zConseqScoreWWAMRange");
	}
	return zConseqScoreWWAM(scanner,start+scanner->model->focus);
}

static score_t zConseqScoreWWAMFeature (zScanner *scanner, zSfeature *f) {
	return zConseqScoreWWAMRange(scanner,f->start,f->end);
}

static score_t zEstseqScoreWWAM (zScanner *scanner, coor_t pos) {
	coor_t  i, p, index, offset;
	int j;
	coor_t  mfocus = pos - scanner->model->focus;
	score_t score = 0;

	if ((pos < scanner->min_pos) || (pos > scanner->max_pos)) {
		return MIN_SCORE;
	}
	score = 0;

	/* scoring */
	for (i = 0; i < scanner->model->length; i++) {
		index = 0;
		for (j = 0; j < scanner->model->order; j++) {
			p = zPOWER[scanner->model->symbols][scanner->model->order - j -1];
			/* EVAN this is different from the WWAM scoring above */
			index += (p * zGetEstseqS10(scanner->seq->parent,i + j + mfocus));

		}

		offset = i*zPOWER[scanner->model->symbols][scanner->model->order];
		score += scanner->model->data[index+offset];
	}
	
	return score;
}

/*******\
  SAM
\*******/

static score_t zDNAScoreSAM (zScanner *scanner, coor_t pos) {
	coor_t i;
	score_t score;
	coor_t mfocus = pos - scanner->model->focus;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos))
		return MIN_SCORE;
	
	/* scoring */
	score = 0;
	for (i = 0; i < scanner->model->length; i++) {
		score += scanner->subscanner[i].score(&scanner->subscanner[i],
											  i + mfocus);
	}

	return score;
}

/*******\
   SDT
\*******/
static score_t zDNAScoreSDT (zScanner *scanner, coor_t pos) {
	char   c;
	int    i, found, scanner_number;
	coor_t j;
	coor_t mfocus = pos - scanner->model->focus;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos))
		return MIN_SCORE;
	
	/* determine which subscanner to use */
	scanner_number = -1; /* impossible number, used to check error */
	for (i = 0; i < scanner->model->submodels; i++) {
		found = 1; /* will set to false and break */		
		for (j = 0; j < scanner->model->length; j++) {
			if (scanner->subscanner[i].sig[j] == 15) continue;
			
			c = zGetDNAS16(scanner->seq->parent,j + mfocus) | scanner->subscanner[i].sig[j];
			if (c != scanner->subscanner[i].sig[j]) {
				found = 0;
				break;
			}
		}
		
		if (found) {
			scanner_number = i;
			break;
		}
	}
	
	if (scanner_number == -1) zDie("no scanner found? zScoreSDT");
	return scanner->subscanner[scanner_number].score(&scanner->subscanner[scanner_number], pos);
}

/*******\
   MIX
\*******/
static score_t zDNAScoreMIX (zScanner *scanner, coor_t pos) {
	int i;
	score_t score;
	coor_t mfocus = pos - scanner->model->focus;

	/* user defines and boundaries */
	if ((pos < scanner->min_pos) || (pos > scanner->max_pos))
		return MIN_SCORE;
	
	/* scoring */
	score = MIN_SCORE;
	for (i = 0; i < scanner->model->submodels; i++) {
		score = zFloatwiseScoreAdd(score, scanner->subscanner[i].score
								   (&scanner->subscanner[i], i + mfocus));
	}

	return score;
}

/*******\
   ISO
\*******/
static score_t zDNAScoreISO (zScanner *scanner, coor_t pos) {
	int i;
	float gc;

	gc = zGetDNAWindowedGC(scanner->seq->parent,pos);

	i = zGetScannerIso(scanner, gc);
	/*EVAN why the hell was this check here? 
	if (scanner->model->submodel[i].type == CDS || scanner->model->submodel[i].type == SIG) {
		zDie("Cannot call zScoreISO with coor_t on CDS/SIG");
		return MIN_SCORE;
	} 
	else {
		return scanner->subscanner[i].score(&scanner->subscanner[i], pos);
	}
	*/
	return scanner->subscanner[i].score(&scanner->subscanner[i], pos);
}

static score_t zDNAScoreISORange (zScanner *scanner, coor_t start, coor_t end) {
	float gc;
	int i;

	gc = zGetDNAWindowedGC(scanner->seq->parent,(int)((start+end)/2));
	i = 0;
	while (gc > scanner->model->data[i] && i < scanner->model->submodels-1) i++;
	
	/*EVAN why the hell was this check here?
	if ((scanner->model->submodel[i].type == CDS || 
		 scanner->model->submodel[i].type == SIG)) {
		return scanner->subscanner[i].scorer(&scanner->subscanner[i], start,end);
	}
	else{
		zDie("Cannot call zScoreISORange on non-SIG");
		return MIN_SCORE;
	}
	*/
	return scanner->subscanner[i].scorer(&scanner->subscanner[i], start,end);
}

static score_t zDNAScoreISOFeature (zScanner *scanner, zSfeature *f) {
	float gc;
	int i;

	gc = zGetDNAWindowedGC(scanner->seq->parent,(int)((f->start+f->end)/2));
	i = 0;
	while (gc > scanner->model->data[i] && i < scanner->model->submodels-1) i++;

	/*EVAN why the hell was this check here?
	if (scanner->model->submodel[i].type == CDS || 
		scanner->model->submodel[i].type == SIG) {
		return scanner->subscanner[i].scoref(&scanner->subscanner[i], f);
	}
	else{
		zDie("Cannot call zScoreISOFeature with zSFeature on non-CDS/non-SIG");
		return MIN_SCORE;
	}*/
	return scanner->subscanner[i].scoref(&scanner->subscanner[i], f);
}

/*********\
   Range
\*********/

score_t zScoreRange (zScanner *scanner, coor_t start, coor_t end) {
	coor_t  i;
	score_t s, score = 0;
	
	if(scanner->uscore_map == NULL){
		/*EVAN not scoring last base for compatability, this will need to be fixed, incompatable with zGetRangeScore below */
		for (i = start; i < end; i++) {
			s = scanner->score(scanner, i);
			if (s == MIN_SCORE) continue; /* boundary condition */
			score += s;
		}
		return score;
	}
	else{
		/*EVAN start-1 since we want to include the start score, 
		  this should maybe be moved into zGetRangeScore? */
		return zGetRangeScore(scanner,end,start-1,-1);
	}
}

score_t zScoreFeature (zScanner *scanner, zSfeature *f) {
	if (f->lfrag || f->rfrag || f->frame) {
		zWriteSfeature(stderr, f);
		zDie("zScoreFeature expects lfrag, rfrag, frame to be zero");
	}	
	return zScoreRange(scanner,f->start,f->end);
}

/********\
   CDS
\********/

static score_t zDNAScoreCDS (zScanner *scanner, zSfeature *f) {
	score_t score = 0;
	coor_t i;
	int n;
	coor_t     start, end;
	
	if (f->name == scanner->Einit || f->name == scanner->Exon) {
		start = f->start +6; /* make room for 5th order Markov Model */
		end   = f->end;
	} 
	else if (f->name == scanner->Eterm || f->name == scanner->Esngl) {
		start = f->start +6;
		end   = f->end   -3; /* don't include stop codon */
	} 
	else {
		start = 0; end = 0;
		zDie("attempt to score CDS with non-exon");
	}
		
	if (f->lfrag < 0 || f->lfrag >= 3){
		zDie("zScoreCDS out of bounds (%d)", f->lfrag);
	}

	for (i = start; i <= end; i++) {
		n = (i - start + 4 - f->lfrag) % 3;		
		score += scanner->subscanner[n].score(&scanner->subscanner[n], i);
	}
		
	return score;
}

/*********************\
   ILLEGAL FUNCTIONS - prevents the use of score and count for CDS
\*********************/
static score_t zIllegalScore (zScanner *s, coor_t p) {
	zDie("illegal score in Scanner (%s), at %d", s->model->name, p);
	return 0;
}

static score_t zIllegalPairScore (zScanner *scanner, zDNA* genomic, zDNA* cdna, coor_t genomic_pos, coor_t cdna_pos) {
	zDie("illegal pairscore in Scanner (%s) for %s/%s, at %d, %d", scanner->model->name, genomic->def, cdna->def, genomic_pos, cdna_pos);
	return 0;
}

static score_t zIllegalScoreRange (zScanner *s, coor_t start, coor_t end) {
	zDie("illegal range score in Scanner (%s), from %d to %d", s->model->name,start,end);
	return 0;
}


/********\
  Init
\********/
int scanner_id = 0;

static void zInitScannerReal(zScanner *scanner, zSequence *seq, 
							 zModel *model, zScanner *parent) {
	int i;
	coor_t j;
	
	scanner->id = scanner_id++;
	/* set dna and model, clear pointers */
	scanner->seq         = seq;
	scanner->model       = model;
	scanner->subscanner  = NULL;
	scanner->sig         = NULL;
	scanner->score       = NULL;
	scanner->scoref      = NULL;
	scanner->pairscore   = NULL;
	scanner->uscore_map  = NULL;
	scanner->uscore_map_size = 0;
	scanner->uscore      = NULL;
	scanner->uscore_block_count = 0;
	scanner->parent      = parent;
	scanner->vars        = NULL;

	/* determine scoring region */
	scanner->min_pos = model->length;
	scanner->max_pos = scanner->seq->length - model->length - 1;
	scanner->min_gpos = scanner->min_pos; /* Only used for PAIR models to check both min_positions. All other models use min_pos only */

	if(model->seq_type == DNA || model->seq_type == GENOMIC) {
		/* bind scoring and counting functions to type of model */
		switch (model->type) {
		case WMM:  scanner->score = zDNAScoreWMM; break;
		case LUT:  scanner->score = zDNAScoreLUT; break;
		case WWAM: scanner->score = zDNAScoreWWAM; break;
		case SAM:  scanner->score = zDNAScoreSAM; break;
		case SDT:  scanner->score = zDNAScoreSDT; break;
		case MIX:  scanner->score = zDNAScoreMIX; break;
		case ISO:  scanner->score = zDNAScoreISO; break;
		default:   scanner->score = zIllegalScore;
		}
		/* bind range scoring functions */
		switch (model->type) {
		case CDS:
			scanner->Einit  = zChar2StrIdx("Einit");
			scanner->Exon   = zChar2StrIdx("Exon");
			scanner->Esngl  = zChar2StrIdx("Esngl");
			scanner->Eterm  = zChar2StrIdx("Eterm");
			scanner->scoref = zDNAScoreCDS; 
			scanner->scorer = zIllegalScoreRange; 
			break;
		case SIG: 
			scanner->scoref = zDNAScoreSIG;
			scanner->scorer = zDNAScoreSIGRange; 
			break;
		case ISO:
			scanner->scoref = zDNAScoreISOFeature; 
			scanner->scorer = zDNAScoreISORange; 
			break;
		default:  
			scanner->scoref = zScoreFeature;
			scanner->scorer = zScoreRange;
		}
	} else if (model->seq_type == PAIR) {
		/* Only LUT and MLUT (for quality values) are allowed as of now */
		scanner->min_pos = PADDING - 1 + model->length/2;
		scanner->max_pos = scanner->seq->length - PADDING - 1;	

		switch (model->type) {
		case LUT:  
			if (model->length%2 != 0) {
				zDie("PAIR LUT model length must be even");
			}
			scanner->pairscore  = zDNAScorePairLUT; 
			break;
		case WMM:
		case WWAM:
		case SAM:
		case SDT: 
		case MIX:
		case ISO:
		default:   scanner->pairscore  = zIllegalPairScore;
		}
	} else if(model->seq_type == CONSEQ) {
		/* bind scoring and counting functions to type of model */
		switch (model->type) {
		case LUT:  scanner->score = zConseqScoreLUT;  break;
		case WWAM: scanner->score = zConseqScoreWWAM; break;
		default:   scanner->score = zIllegalScore;
		}
		/* bind range scoring functions */
		switch (model->type) {
		case WWAM: 
			scanner->scoref = zConseqScoreWWAMFeature;
			scanner->scorer = zConseqScoreWWAMRange;
			break;
		default:  
			scanner->scoref = zScoreFeature;
			scanner->scorer = zScoreRange;
		}
	}
	else if(model->seq_type == ESTSEQ){
		/* bind scoring and counting functions to type of model */
		switch (model->type) {
		case LUT:  scanner->score = zEstseqScoreLUT; break;
		case WWAM: scanner->score = zEstseqScoreWWAM; break;
		default:   scanner->score = zIllegalScore;
		}
		/* bind range scoring functions */
		switch (model->type) {
		default:  
			scanner->scoref = zScoreFeature;
			scanner->scorer = zScoreRange;
		}
	}
	else {
		zDie("unrecognized model type for model %s\n", model->name);
	}

	/* create subscanners for constructed types */
	if (model->type == WMM || model->type == LUT || model->type == WWAM || 
		model->type == SIG) {
		scanner->subscanner = NULL;
	} 
	else {
		scanner->subscanner = zMalloc(model->submodels * sizeof(zScanner),
									  "zInitScanner subscanner");
		for (i = 0; i< model->submodels; i++){
			zInitScannerReal(&scanner->subscanner[i], seq, &model->submodel[i],
							 scanner);
		}
	}
	
	/* ensure CDS model is correctly used */
	if (model->type == CDS) {
		if (model->submodels != 3) zDie("CDS must have 3 submodels");
		for (i = 0; i < model->submodels; i++) {
			if (model->submodel[i].type != LUT) zDie("CDS submodels not LUT");
		}
	}
	
	/* create decision tree for SDT */
	if (model->type == SDT) {
		
		/* create sigs for submodels */
		scanner->sig = NULL; /* zMalloc(model->submodels,"zInitScanner sig");*/
		for (i = 0; i < model->submodels; i++) {
			scanner->subscanner[i].sig = zMalloc(model->length,
												 "zInitScanner sig 2");
			for (j = 0; j < model->length; j++) {
				scanner->subscanner[i].sig[j] = seq2sig(model->submodel[i].name[j]);
				/* name contains the SDT sequence, so name[j] is the j'th alphabet 
				   in the sequence */
			}
		}		
	} 
	else {
		scanner->sig = NULL;
	}
}

void zInitScanner(zScanner *scanner, zSequence *seq, zModel *model){
	zInitScannerReal(scanner,seq,model,NULL);
}

void zFreeScanner(zScanner *scanner) {
	int i;

	zDePreComputeScanner(scanner);

	for (i = 0; i < scanner->model->submodels; i++)
		zFreeScanner(&scanner->subscanner[i]);
	zFree(scanner->subscanner);	

	zFree(scanner->sig);

	scanner->seq = NULL;
	scanner->model = NULL;
	scanner->sig = NULL;
	scanner->subscanner = NULL;
	scanner->score = NULL;
	scanner->scoref = NULL;
}

static int zModelGetPreEffectRange(zModel* model){
	int i;
	int pre = model->focus;
	int	new_pre;
	for(i = 0;i < model->submodels;i++){
		new_pre = zModelGetPreEffectRange(&model->submodel[i]);
		if(new_pre > pre){
			pre = new_pre;
		}
	}
	return pre;
}

static int zModelGetPostEffectRange(zModel* model){
	int i;
	int post = model->length - model->focus - 1;
	int	new_post;
	for(i = 0;i < model->submodels;i++){
		new_post = zModelGetPostEffectRange(&model->submodel[i]);
		if(new_post > post){
			post = new_post;
		}
	}
	return post;
}

static void zPreComputeScannerMap (zScanner* scanner){
	zScannerVariant   *var;
	int                first_idx,i,v;
	zPtrListPos        list_pos;
	coor_t             pre,post;

	/* get the block size */
	if(scanner->model->type == CDS){
		scanner->uscore_block_size = CDS_BLOCK_SIZE;
	}
	else{
		scanner->uscore_block_size = NONCDS_BLOCK_SIZE;
	}
	scanner->uscore_block_count = USCORE_BLOCK_COUNT;
	scanner->uscore_map_size = 
		zIntMax(scanner->uscore_block_count,
				scanner->seq->length/scanner->uscore_block_size + 1);
	scanner->uscore_map = 
		zMalloc(sizeof(zUScoreMap)*scanner->uscore_map_size,
				"zInitScanner uscore_map");
	/* fill each map */
	for(i = 0;i < scanner->uscore_map_size;i++){
		scanner->uscore_map[i].pos = i*scanner->uscore_block_size;
		scanner->uscore_map[i].vars = NULL;
		scanner->uscore_map[i].var_count = 0;
		scanner->uscore_map[i].block_idx = -1;
	}
	/* get pointers to sequecnce variants in each uscore_map */
	if(scanner->seq->var_count > 0){
		scanner->vars = zMalloc(sizeof(zPtrList),
								"zPreComputeScanner vars");
		zInitPtrList(scanner->vars);
		/* get model range */
		pre = zModelGetPreEffectRange(scanner->model);
		post = zModelGetPostEffectRange(scanner->model);
		/* create new variant */
		var = zMalloc(sizeof(zScannerVariant),
					  "zPreComputeScanner var");
		zPtrListAddLast(scanner->vars,var);
		var->var_count = 1;
		first_idx = 0;
		var->first_pos = scanner->seq->variants[0]->pos - post;
		var->last_pos = scanner->seq->variants[0]->pos + pre;
		for(v = 1;v < scanner->seq->var_count;v++){
			/* if the next one overlaps the current block */
			if(var->last_pos >= scanner->seq->variants[v]->pos - post){
				var->last_pos = scanner->seq->variants[v]->pos + pre;
				var->var_count++;
			}
			else{
				/* store the seq variants in this scanner variant */
				var->vars = zMalloc(sizeof(zSeqVariant*)*var->var_count,
									"zPreComputeScanner var->vars");
				var->var_vals = 1;
				for(i = 0; i < var->var_count; i++){
					var->vars[i] = scanner->seq->variants[first_idx+i];
					var->var_vals *= var->vars[i]->variants; 
				}
				/* allocate score array */
				var->score = zMalloc(sizeof(score_t*)*var->var_vals,
									 "zPreComputeScanner var->score");
				var->score_sum = zMalloc(sizeof(score_t*)*var->var_vals,
										 "zPreComputeScanner var->score_sum");
				/* allocate and fill in on first request */ 
				for(i = 0; i < var->var_vals; i++){ 
					var->score[i] = NULL;
					var->score_sum[i] = NULL;
				}
				/* create new variant */
				var = zMalloc(sizeof(zScannerVariant),
							  "zPreComputeScanner var");
				zPtrListAddLast(scanner->vars,var);
				var->var_count = 1;
				first_idx = v;
				var->first_pos = scanner->seq->variants[v]->pos - post;
				var->last_pos = scanner->seq->variants[v]->pos + pre;
			}
		}
		/* finish the final var we created */
		/* store the seq variants in this scanner variant */
		var->vars = zMalloc(sizeof(zSeqVariant*)*var->var_count,
							"zPreComputeScanner var->vars");
		var->var_vals = 1;
		for(i = 0; i < var->var_count; i++){
			var->vars[i] = scanner->seq->variants[first_idx+i];
			var->var_vals *= var->vars[i]->variants; 
		}
		/* allocate score array */
		var->score = zMalloc(sizeof(score_t*)*var->var_vals,
							 "zPreComputeScanner var->score");
		var->score_sum = zMalloc(sizeof(score_t*)*var->var_vals,
								 "zPreComputeScanner var->score_sum");
		/* allocate and fill in on first request */ 
		for(i = 0; i < var->var_vals; i++){ 
			var->score[i] = NULL;
			var->score_sum[i] = NULL;
		}

		/* store var pointers in appropriate map block */
		for(i = 0;i < scanner->uscore_map_size;i++){
			var = zPtrListMoveFirst(scanner->vars);
			scanner->uscore_map[i].var_count = 0;
			while(var != NULL && var->last_pos < scanner->uscore_map[i].pos){
				var = zPtrListMoveNext(scanner->vars);
			}
			list_pos = zPtrListGetPos(scanner->vars);
			while(var != NULL && 
				  var->first_pos < 
				  scanner->uscore_map[i].pos + scanner->uscore_block_size){
				scanner->uscore_map[i].var_count++;
				var = zPtrListMoveNext(scanner->vars);
			}
			if(scanner->uscore_map[i].var_count > 0){
				scanner->uscore_map[i].vars = 
					zMalloc(sizeof(zScannerVariant*)*scanner->uscore_map[i].var_count,
							"zPreCoputeScanner uscore_map[i].vairants");
				zPtrListSetPos(scanner->vars,list_pos);
				v = 0;
				var = zPtrListGetCurrent(scanner->vars);
				while(v < scanner->uscore_map[i].var_count){
					scanner->uscore_map[i].vars[v++] = var;
					var = zPtrListMoveNext(scanner->vars);
				}
			}
		}
	}
}

void zPreComputeScanner(zScanner* scanner) {
	coor_t    i,k;
	int       j;

	if(scanner->uscore_map != NULL) return;

	/* create uscore map or link to parents map */	
	/* ISO submodels will rarely need blocks to cover the same sequence 
	   so we give each its own map */
	if((scanner->parent == NULL) || (scanner->parent->model->type == ISO)){
		zPreComputeScannerMap(scanner);
	}
	else{
		/* point to the parents map/info */
		scanner->uscore_map = scanner->parent->uscore_map;
		scanner->uscore_map_size = scanner->parent->uscore_map_size;
		scanner->uscore_block_count = scanner->parent->uscore_block_count;
		scanner->uscore_block_size = scanner->parent->uscore_block_size;
	}
	

	if (ISO == scanner->model->type) {
		for (i = 0; (int)i < scanner->model->submodels; i++) {
			zPreComputeScanner (&scanner->subscanner[i]);
		}
	} 
	else if (CDS == scanner->model->type) {
		if (NULL != scanner->subscanner[0].uscore) return;

		zInitPtrList(&scanner->uscore_list);
		scanner->uscore = 
			zMalloc(sizeof(zUScoreBlock)*scanner->uscore_block_count,
					"zPreComputeScanner scanner->uscore (CDS)");
		for(i = 0;(int)i < scanner->uscore_block_count;i++){
			scanner->uscore[i].score = NULL; 
			scanner->uscore[i].score_sum = NULL; 
			scanner->uscore[i].var_map = 
				zMalloc(sizeof(int)*scanner->uscore_block_size,
						"zPrecomputeScanner scanner->uscore[i].var_map (CDS)");
			for(k = 0;k < scanner->uscore_block_size;k++){
				scanner->uscore[i].var_map[k] = 0;
			}
			scanner->uscore[i].map_idx = -1;
			scanner->uscore[i].idx = i;
		}

		for(j = 0; j < 3; j++){		
			scanner->subscanner[j].uscore_map = scanner->uscore_map;
			scanner->subscanner[j].uscore = 
				zMalloc(sizeof(zUScoreBlock)*scanner->uscore_block_count,
						"zPreComputeScanner scanner->subscanner[j]->uscore");
			for(i = 0;(int)i < scanner->uscore_block_count;i++){
				scanner->subscanner[j].uscore[i].score = NULL;
				scanner->subscanner[j].uscore[i].var_map = NULL;
				scanner->subscanner[j].uscore[i].score_sum = 
					zMalloc(sizeof(score_t)*scanner->uscore_block_size,
							"zPrecomputeScanner scanner->uscore[i].score_sum (CDS)");
				scanner->subscanner[j].uscore[i].map_idx = -1;
				scanner->subscanner[j].uscore[i].idx = i;
			}
		}
		for(j = 0; j < zIntMin(scanner->uscore_block_count,scanner->uscore_map_size); j++){		
			zFillUScore(scanner,j,j);
			zPtrListAddLast(&scanner->uscore_list,&scanner->uscore[j]);
		}
	} 
	else {

		if (NULL != scanner->uscore) return;

		zInitPtrList(&scanner->uscore_list);
		scanner->uscore = zMalloc(sizeof(zUScoreBlock)*scanner->uscore_block_count,
								  "zPreComputeScanner scanner->uscore");
		for(i = 0;(int)i < scanner->uscore_block_count;i++){
			scanner->uscore[i].score = 
				zMalloc(sizeof(score_t)*scanner->uscore_block_size,
						"zPrecomputeScanner scanner->uscore[i].score");
			scanner->uscore[i].score_sum = 
				zMalloc(sizeof(score_t)*scanner->uscore_block_size,
						"zPrecomputeScanner scanner->uscore[i].score_sum");
			scanner->uscore[i].var_map = 
				zMalloc(sizeof(int)*scanner->uscore_block_size,
						"zPrecomputeScanner scanner->uscore[i].var_map");
			for(k = 0;k < scanner->uscore_block_size;k++){
				scanner->uscore[i].var_map[k] = 0;
			}
			/* handle seq variant */
			scanner->uscore[i].map_idx = -1;
			scanner->uscore[i].idx = i;
			zFillUScore(scanner,i,i);
			zPtrListAddLast(&scanner->uscore_list,&scanner->uscore[i]);
		}
	}
}

static void zDePreComputeScannerMap(zScanner* scanner) {
	int i;
	zScannerVariant *var;
	if(scanner->vars != NULL){
		var = zPtrListMoveFirst(scanner->vars);
		while(var != NULL){
			for(i = 0; i < var->var_vals; i++){ 
				if(var->score[i] != NULL){
					zFree(var->score[i]);
					var->score[i] = NULL;
				}
			}
			zFree(var->score);
			var->score = NULL;
			zFree(var->vars);
			var->vars = NULL;
			zFree(var);
			var = zPtrListMoveNext(scanner->vars);
		}
		zFreePtrList(scanner->vars);
		zFree(scanner->vars);
		scanner->vars = NULL;
		for (i = 0; i < scanner->uscore_map_size; i++) {
			if(scanner->uscore_map[i].var_count > 0){
				zFree(scanner->uscore_map[i].vars);
				scanner->uscore_map[i].vars = NULL;
			}
		}
	}
	zFree(scanner->uscore_map);
	scanner->uscore_map = NULL;
}

void zDePreComputeScanner(zScanner* scanner) {
	int i,j;

	if((scanner->parent == NULL) || (scanner->parent->model->type == ISO)){
		zDePreComputeScannerMap(scanner);
	}

	if (ISO == scanner->model->type) {
		for (i = 0; i < scanner->model->submodels; i++) {
			zDePreComputeScanner(&scanner->subscanner[i]);
		}
	} 
	else if (CDS == scanner->model->type) {
		for (i = 0; i < 3; i++) {
			if(scanner->subscanner[i].uscore != NULL){
				for(j = 0;j < scanner->uscore_block_count;j++){
					zFree(scanner->subscanner[i].uscore[j].score_sum);
				}
				zFree(scanner->subscanner[i].uscore);
				scanner->subscanner[i].uscore = NULL;
			}
		}
	} 
	
	if(scanner->uscore != NULL){
		for(i = 0;i < scanner->uscore_block_count;i++){
			if(scanner->uscore[i].score != NULL){
				zFree(scanner->uscore[i].score);
			}
			if(scanner->uscore[i].score_sum != NULL){
				zFree(scanner->uscore[i].score_sum);
			}
			if(scanner->uscore[i].var_map != NULL){
				zFree(scanner->uscore[i].var_map);
			}
		}
		zFree(scanner->uscore);
		scanner->uscore = NULL;
	}
}

int zGetScannerIso(zScanner *scanner, float gc) {
	int i = 0;
	while (gc > scanner->model->data[i] && i < scanner->model->submodels) i++;
	return i;
}

zScanner* zGetScannerForIso(zScanner *scanner, float gc) {
	return &scanner->subscanner[zGetScannerIso(scanner, gc)];
}

static score_t zScoreRangeBlock(zScanner* scanner,zUScoreBlock* block, coor_t end, 
								 coor_t start, int frame){
	zUScoreMap      *map;
	zScannerVariant *var;
	score_t          score;
	score_t         *var_scores;
	int              v;

	if(frame == -1){
		frame = 1;
	}

	if(end >= block->pos + scanner->uscore_block_size){
		end = block->pos + scanner->uscore_block_size - 1;
	}

	if(scanner->seq->var_count == 0){	
		score = block->score_sum[end-block->pos];
		if(start >= block->pos){
			score -= block->score_sum[start-block->pos];
		}
	}
	else{
		map = &scanner->uscore_map[block->map_idx];
		score = 0;
		while(end > start){
			/* get CDS scanner block rather than LUT scanner block */
			v = scanner->uscore[map->block_idx].var_map[end-block->pos];
			if(v == 0){
				score += block->score_sum[end-block->pos];
				if(start > block->pos){
					score -= block->score_sum[start-block->pos];
				}
				break;
			}
			else if(v > 0){
				v--;/* -1 readjust to actual index */
				var = map->vars[v]; 
				var_scores = zFillUScoreVariant(scanner,block->map_idx,block->idx,v,1);
				var_scores = &(var_scores[frame*(var->last_pos-var->first_pos+1)]);
				score += var_scores[end-var->first_pos];
				if(start > var->first_pos){
					score -= var_scores[start-var->first_pos];
				}
				end = var->first_pos-1;
			}
			else if(v < 0){
				score += block->score_sum[end-block->pos];
				end += v;
				if(end < start){
					score -= block->score_sum[start-block->pos];
				}
				else{
					score -= block->score_sum[end-block->pos];
				}
			}
		}
	}
	return score;
}

score_t zGetRangeScore(zScanner *scanner,coor_t end, coor_t start, int frame){
	zScanner *s;
	zUScoreBlock* block;
	coor_t pos;
	int first_map,last_map,i;
	score_t score = 0;

	if(end < start){
		return MIN_SCORE;
	}

	if(scanner->model->type == ISO){
		pos = (end+start)/2;
		scanner = zGetScannerForIso(scanner,zGetDNAWindowedGC(scanner->seq->parent,pos));
	}
	if(frame == -1){
		s = scanner;
	}
	else{
		s = &scanner->subscanner[frame];
	}

	first_map = (int)(start/scanner->uscore_block_size);
	last_map = (int)(end/scanner->uscore_block_size);
	if(scanner->uscore_map[first_map].pos > start ||
	   scanner->uscore_map[first_map].pos + scanner->uscore_block_size <= start){
		zDie("Damn1");
	}
	if(scanner->uscore_map[last_map].pos > end ||
	   scanner->uscore_map[last_map].pos + scanner->uscore_block_size <= end){
		zDie("Damn2");
	}

	for(i = first_map; i <= last_map; i++){
		scanner->uscore_map[i].tag = 0;
		if(scanner->uscore_map[i].block_idx == -1){
			scanner->uscore_map[i].tag = 1;
		}
		else{
			/* get the block for this frame */
			block = &s->uscore[scanner->uscore_map[i].block_idx];
			score += zScoreRangeBlock(scanner,block,end,start,frame);
		}
	}
	for(i = last_map; i >= first_map; i--){
		if(scanner->uscore_map[i].tag == 1){
			block = zPtrListMoveLast(&scanner->uscore_list);
			zPtrListMoveCurrentToFirst(&scanner->uscore_list);
			/* get the score for this frame */
			zFillUScore(scanner,i,block->idx);
			block = &s->uscore[block->idx];
			score += zScoreRangeBlock(scanner,block,end,start,frame);
		}
	}
	return score;
}

static zUScoreBlock* zGetUScoreBlock(zScanner* scanner, coor_t pos){
	zUScoreBlock* block;
	int i;

	block = zPtrListMoveFirst(&scanner->uscore_list);
	while(block != NULL){
		if((block->pos <= pos) && (block->pos + scanner->uscore_block_size > pos)){
			zPtrListMoveCurrentToFirst(&scanner->uscore_list);
			return block;
		}
		block = zPtrListMoveNext(&scanner->uscore_list);
	}
	/* load new block */
	block = zPtrListMoveLast(&scanner->uscore_list);
	zPtrListMoveCurrentToFirst(&scanner->uscore_list);
	i = (int)pos/scanner->uscore_block_size;

	zFillUScore(scanner,i,block->idx);
	return block;
}

score_t zGetUScore(zScanner *scanner,coor_t i){
	zUScoreBlock *block;
	zUScoreMap   *map;
	coor_t        idx;
	int           v;
	score_t      *scores;

	if(i < scanner->min_pos || i > scanner->max_pos){
		return MIN_SCORE;
	}

	block = zGetUScoreBlock(scanner,i);
	map = &scanner->uscore_map[block->map_idx];
	idx = i - block->pos;
	v = block->var_map[idx];
	if(v == 0){
		return block->score[idx];
	}
	else{
		v--; /* -1 readjust to actual index */
		scores = map->vars[v]->score[zGetScannerVariantScoreIndex(map->vars[v])];
		if(scores == NULL){
			scores = zFillUScoreVariant(scanner,block->map_idx,block->idx,v,0);
		}
		idx = i - map->vars[v]->first_pos;
		return scores[idx];
	}
}

static int zGetScannerVariantScoreIndex(zScannerVariant* var){
	int i,var_val;
	var_val = 0;
	for(i = 0; i < var->var_count; i++){
		var_val += (var->vars[i]->cur * (i+1));
	}
	return var_val;
}

static void zFillCDSUScore(zScanner* scanner,coor_t start,coor_t size,
						   score_t *pre1,score_t  *pre2,score_t  *pre3){
	int       f0, f1, f2;
	coor_t    i,j;
	zScanner* s;
	
	i = start;
	
	f0 = (i - 1) % 3;
	f1 = (i + 0) % 3;
	f2 = (i + 1) % 3;
	
	s = &scanner->subscanner[f0];
	pre1[0] = s->score(s, i);
	s = &scanner->subscanner[f1];
	pre2[0] = s->score(s, i);
	s = &scanner->subscanner[f2];
	pre3[0] = s->score(s, i);
	
	if (pre1[0] == MIN_SCORE) pre1[0]  = 0;
	if (pre2[0] == MIN_SCORE) pre2[0]  = 0;
	if (pre3[0] == MIN_SCORE) pre3[0]  = 0;
	
	for(j = 1;j < size;j++){
		
		i++;
		
		f0 = (i - 1) % 3;
		f1 = (i + 0) % 3;
		f2 = (i + 1) % 3;
		
		s = &scanner->subscanner[f0];
		pre1[j] = s->score(s, i);
		s = &scanner->subscanner[f1];
		pre2[j] = s->score(s, i);
		s = &scanner->subscanner[f2];
		pre3[j] = s->score(s, i);
		
		if (pre1[j] == MIN_SCORE) pre1[j]  = 0;
		else pre1[j] += pre1[j-1];
		if (pre2[j] == MIN_SCORE) pre2[j]  = 0;
		else pre2[j] += pre2[j-1];
		if (pre3[j] == MIN_SCORE) pre3[j]  = 0;
		else pre3[j] += pre3[j-1];
	}
}

static score_t* zFillUScoreVariant(zScanner* scanner,int map_idx,int block_idx,int var_idx,
								   int sum){
	zUScoreBlock*    block;
	zUScoreMap*      map;
	zScannerVariant* var;
	int              v;
	coor_t           size,i;

	map = &scanner->uscore_map[map_idx];
	block = &scanner->uscore[block_idx];
	var = map->vars[var_idx];
	
	v = zGetScannerVariantScoreIndex(var);
	if(sum){
		if(var->score_sum[v] != NULL) return var->score_sum[v];
	}
	else{
		if(var->score[v] != NULL) return var->score[v];
	}
	size = var->last_pos - var->first_pos + 1;

	/* fill in the scores */
	if (CDS == scanner->model->type) {
		var->score_sum[v] = zMalloc(sizeof(score_t)*3*size,
									"zPreComputeScanner var->score[i] (CDS)");
		zFillCDSUScore(scanner,var->first_pos, size,
					   &(var->score_sum[v][0]),&(var->score_sum[v][size]),
					   &(var->score_sum[v][2*size]));
	}
	else{
		var->score[v] = zMalloc(sizeof(score_t)*size,
								"zPreComputeScanner var->score[i]");
		var->score_sum[v] = zMalloc(sizeof(score_t)*size,
									"zPreComputeScanner var->score_sum[i]");
		var->score[v][0] = scanner->score(scanner,var->first_pos);
		if(var->score[v][0] == 0){
			var->score_sum[v][0] = 0;
		}
		else{
			var->score_sum[v][0] = zGetUScore(scanner,var->first_pos-1) + 
				var->score[v][0];
		}
		for(i = var->first_pos+1;i <= var->last_pos; i++){
			var->score[v][i - var->first_pos] = scanner->score(scanner,i);
			/* EVAN included for consistency.  see comment in zFillUScore below */
			if(var->score[v][i-var->first_pos] == MIN_SCORE){
				var->score_sum[v][i - var->first_pos] = 0;
			}
			else{
				var->score_sum[v][i - var->first_pos] = var->score_sum[v][i-1-var->first_pos] +  
					var->score[v][i - var->first_pos];
			}
		}
	}
	if(sum){
		return var->score_sum[v];
	}
	else{
		return var->score[v];
	}
}

/* fill scanner->uscore[block_idx] with the block described by 
   scanner->uscore_map[map_idx] */
static void zFillUScore(zScanner *scanner,int map_idx,int block_idx){
	zUScoreBlock  *block;
	zUScoreMap    *map;
	score_t       *pre1,*pre2,*pre3;
	int            v;
	coor_t         i,j,last,last_pos;
	coor_t         first_idx,last_idx;

	if (ISO == scanner->model->type) {
		/* kind of stupid to call on an ISO scanner since each blocks will 
           rarely be in two or more ISO catagories and this fills a block for
		   all isochores */
		for (i = 0; (int)i < scanner->model->submodels; i++) {
			zFillUScore(&scanner->subscanner[i],map_idx,block_idx);
		}
	} 
	else {
		map = &scanner->uscore_map[map_idx];
		block = &scanner->uscore[block_idx];

		/* update block/map indices */
		if(block->map_idx != -1){
			scanner->uscore_map[block->map_idx].block_idx = -1;
		}
		map->block_idx = block_idx;
		block->map_idx = map_idx;
		block->pos = map->pos;

		if (CDS == scanner->model->type) {
			for(j = 0;j < 3; j++){
				scanner->subscanner[j].uscore[block_idx].map_idx = map_idx;
				scanner->subscanner[j].uscore[block_idx].pos = map->pos;
				scanner->subscanner[j].uscore_map[map_idx].block_idx = 
					block_idx;
			}

			/* get pointers to the score arrays */
			pre1 = scanner->subscanner[0].uscore[block_idx].score_sum;
			pre2 = scanner->subscanner[1].uscore[block_idx].score_sum;
			pre3 = scanner->subscanner[2].uscore[block_idx].score_sum;

			zFillCDSUScore(scanner,map->pos,scanner->uscore_block_size,
						   pre1, pre2, pre3);
		
			if(scanner->seq->var_count > 0){
				/* set variant map */
				if(map->var_count > 0){
					for(i = block->pos; i < map->vars[0]->first_pos; i++){
						block->var_map[i-block->pos] = 0;
					}
					last_pos = zIntMax(map->vars[0]->first_pos-1,block->pos-1);
					for(v = 0;v < map->var_count; v++){
						last = -1;
						for(i = last_pos+1;i < map->vars[v]->first_pos; i++){
							block->var_map[i-block->pos] = last--;
						}
						first_idx = zIntMax(map->vars[v]->first_pos,block->pos);
						last_idx = zIntMin(map->vars[v]->last_pos,
										   block->pos+scanner->uscore_block_size-1);
						for(i = first_idx;i <= last_idx;i++){
							/* v+1 because 0 is special "no more vars" symbol */
							block->var_map[i-block->pos] = v+1;
						}
						last_pos = map->vars[v]->last_pos;
					}
					last = -1;
					for(i = last_pos+1-block->pos; i < scanner->uscore_block_size; i++){
						block->var_map[i] = last--;
					}
				}
				else{
					for(i = 0; i < scanner->uscore_block_size; i++){
						block->var_map[i] = 0;
					}
				}
			}
		}
		else {
			/* score block */
			block->score[0] = scanner->score(scanner,block->pos);				
			if(block->score[0] == MIN_SCORE){
				block->score_sum[0] = 0;
			}
			else{
				block->score_sum[0] = block->score[0];				
			}
			for(i = 1;i < scanner->uscore_block_size;i++){
				block->score[i] = scanner->score(scanner,i+block->pos);
				/*EVAN included for consistency but this is wierd
				 without this we get no introns though */
				if(block->score[i] == MIN_SCORE){
					block->score_sum[i] = 0;
				}
				else{
					block->score_sum[i] = block->score_sum[i-1] + block->score[i];
				}
			}
			if(scanner->seq->var_count > 0){
				/* reset var values */
				for(i = 0;i < scanner->uscore_block_size;i++){
					block->var_map[i] = 0;
				}			
				/* set variant map */
				for(v = 0;v < map->var_count; v++){
					first_idx = zIntMax(map->vars[v]->first_pos,block->pos);
					last_idx = zIntMin(map->vars[v]->last_pos,
									   block->pos+scanner->uscore_block_size-1);
					for(i = first_idx; i <= last_idx; i++){
						block->var_map[i-block->pos] = v+1;
					}
				}
			}
		}
	}
}
#endif
