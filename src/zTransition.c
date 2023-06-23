/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zTransition.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_TRANSITION_C
#define ZOE_TRANSITION_C

#include "zTransition.h"
#include "zTrellis.h"
#include "zFeatureFactory.h" /* for PSC definition */

void zSetScanners(zScanner**, zScanner**, zScanner**,
                  zScanner**,
                  zAlignmentScanner **, zAlignmentScanner **,
                  zTrellis*, zStrIdx, bool);
score_t zScoreScanners(zTrellis *, zStrIdx, const char *, coor_t,
                       coor_t, bool, int *, int *);
void zSubInternalTrans(zTrellis *, int, int, coor_t, int, int, int,
                       score_t, zDistribution *);
void zSubExplicitTrans(zTrellis *, int, int, coor_t, int, coor_t,
                       coor_t, zDistribution *);

int zReadTransition (FILE *stream, zTransition *t, int iso_groups) {
	char          from[17], to[17];
	float         *prob;
	int i;

	prob = zCalloc(iso_groups, sizeof(float), "zReadTransition Isochores prob");
	if (fscanf(stream, "%16s %16s %f", from, to, prob) != 3) {
		zWarn("fscanf error in zReadTransition");
		return 0;
	}
	for (i = 1; i < iso_groups; i++){
		if(fscanf(stream, "%f", &prob[i]) != 1) {
			zWarn("fscanf error in zReadTransition Isochores");
			return 0;
		}
	}

	t->prob = prob;
	t->score = zCalloc(iso_groups, sizeof(score_t), "zReadTransition Isochores score");
	for (i = 0; i < iso_groups; i++){
		t->score[i] = zFloat2Score(prob[i]);
	}

	/* from */
	if (!zStrIdxExists(from)) {
		zWarn("zReadTransition from unrecognized state (%s)", from);
		return 0;
	}
	t->from = zChar2StrIdx(from);

	/* to */
	if (!zStrIdxExists(to)) {
		zWarn("zReadTransition from unrecognized state (%s)", to);
		return 0;
	}
	t->to = zChar2StrIdx(to);

	return 1;
}

void zWriteTransition (FILE *stream, const zTransition *t, int iso_groups) {
	int i;
	(void)fprintf(stream, "%s\t%s\t%f", zStrIdx2Char(t->from),
				  zStrIdx2Char(t->to), t->prob[0]);
	for(i = 1; i < iso_groups; i++){
		(void) fprintf(stream, "\t%f", t->prob[i]);
	}
	(void) fprintf(stream, "\n");
}

/******************************************************************************\
 Transition Functions
 
\******************************************************************************/

int zCompatibleJumpFromExon (zDNA *dna, zPhase_t phase, zSfeature *exon) {
	/* this function determines if an exon is compatible with an intron */
	int terminal_base;
	coor_t end;
	
	end = ('-' == exon->strand) ? dna->length - 1 - exon->start : exon->end;
		
	switch (exon->rfrag) {
	case 1:  terminal_base = zGetDNAS5(dna,end) +1; break;
	case 2:  terminal_base = (zGetDNAS5(dna,end -1) +1) * 5
				 +  zGetDNAS5(dna,end) +1; break;
	default: terminal_base = 0;
	}

	switch (phase) {
	case Phase0: /* must have no terminal bases */
		switch (terminal_base) {
		case 0:  return 1;
		default: return 0;
		}
	case Phase1: /* must have one terminal, but not T */
		switch (terminal_base) {
		case 1: case 2: case 3: case 5: return 1;
		default:                        return 0;
		}
	case Phase1T:
		switch (terminal_base) {
		case 4:  return 1;
		default: return 0;
		}
	case Phase2: /* must have two terminals, but not TA or TG  */
		if (terminal_base > 5
			&& terminal_base != 21
			&& terminal_base != 23) return 1;
		else                        return 0;
	case Phase2TA: /* must have terminal TA */
		switch (terminal_base) {
		case 21: return 1;
		default: return 0;
		}
	case Phase2TG: /* must have terminal TG */
		switch (terminal_base) {
		case 23: return 1;
		default: return 0;
		}
	default: return 1;
	}
}

int zCompatibleJumpToExon (zDNA *dna, zPhase_t phase, zSfeature *exon) {
	int initial_base;
	coor_t start;
	
	start = ('-' == exon->strand) ? dna->length - 1 - exon->end : exon->start;
	
	switch (exon->lfrag) {
	case 1:  initial_base =  zGetDNAS5(dna,start)   + 1; break;
	case 2:  initial_base = (zGetDNAS5(dna,start)   + 1) *5
				 +  zGetDNAS5(dna,start+1) + 1; break;
	default: initial_base = 0;
	}
	
	switch (phase) {
	case Phase0:   return (0 == initial_base);
	case Phase1:   return (2 == exon->lfrag);
	case Phase1T:  return ((2 == exon->lfrag)  &  /* T + ... */
						   (6 != initial_base) &  /*     AA  */
						   (9 != initial_base) &  /*     AG  */
						   (16 != initial_base)); /*     GA  */
	case Phase2:   return (1 == exon->lfrag);
	case Phase2TA: return ((1 == exon->lfrag)  & /* TA + . */
						   (1 != initial_base) & /*      A */
						   (3 != initial_base)); /*      G */
		
	case Phase2TG: return ((1 == exon->lfrag) &
						   (1 != initial_base)); /* TG + A */
		/*		default: return 1; */
	}
	return 1; /* shush */
}

static void zExplicitTrans (zTrellis *trellis, int from_state, int state, coor_t i) 
{
	/* Transition function for a state with arbitrary length distribution. */
	/* The length distribution can be defined as a single function (or a   */
	/* single range of score values) or as a number of consecutive ranges, */
	/* each defined as a function or a set of score values. INTERNAL state */
	/* (single geometric distribition) is also dealt with here (see below) */
 
	int di; /* indexes ranges of length distribution */
	int iso_group, submax_from_state, submax_length;
	zPhase_t int_phase, ext_phase;
	score_t submax_score;
	zDurationGroup *group;
	zDistribution *d;
	coor_t min, max;

	if (EXTERNAL == trellis->hmm->state[from_state].type) 
		{ /*  check phase */
			int_phase = trellis->hmm->state[state].phase;
			ext_phase = trellis->hmm->state[from_state].phase;
			if (ext_phase != int_phase) return;
		}

	group     = trellis->hmm->dmap[trellis->hmm->state[state].duration];
	iso_group = zGetDurationIsochoreGroup(group, zGetDNAGC(trellis->dna));

	/* Iterate over all ranges of length distribution which are compatible */
	/* with position i. The highest-scoring feature for each range will be */
	/* stored in trellis->cell[state][i].submax->elem[di].                  */

	di = 0;
	while ((di < group->duration[iso_group].distributions) 
		   && ((i-PADDING+1) >= group->duration[iso_group].distribution[di].start))
		{
			d = &group->duration[iso_group].distribution[di];

			/* Set min and max start coordinates for the range */
			max = i - (d->start);
			if (d->end > (i-PADDING+1))
				{
					min = PADDING-1;
				}
			else
				{
					min = i - (d->end);
				}

			if (d->type == GEOMETRIC) /* the range has geometric distribution */
				{
					/* submax_from_state, submax_score and submax_length refer to */ 
					/* the highest-scoring feature in this range at position i-1  */
					
					submax_score = trellis->cell[state][i-1].submax->elem[di].intrinsic_score;
					submax_from_state = trellis->cell[state][i-1].submax->elem[di].from_state;
					submax_length = trellis->cell[state][i-1].submax->elem[di].end 
						- trellis->cell[state][i-1].submax->elem[di].start
						+ 1;
					
					if ((d->end > trellis->dna->length) || (from_state != submax_from_state))
						{
							/* the range extends into infinity OR highest-scoring feature  */
							/* in this range at position i-1 connects to a different state */
							/* INTERNAL states are included in this category               */
							
							zSubInternalTrans(trellis, from_state, state, i, di,
											  submax_from_state, submax_length, submax_score, d);
						}
					else
						{
							/* the range is bounded from above AND optimal feature at */
							/* position i-1 connects to the same state                */

							if (submax_length < ((int) d->end))
								{
									zSubInternalTrans(trellis, from_state, state, i, di,
													  submax_from_state, submax_length, submax_score, d);
								}
							else
								{
									/* At submax_length == d->end the range has to be  */
									/* handled explicitly (can't argue with the math!) */

									zSubExplicitTrans(trellis, from_state, state, i, di, min, 
													  max, d);
								}
						}
				}
			else /* all non-geometric distributions are handled explicitly */
				{
					zSubExplicitTrans(trellis, from_state, state, i, di, min, max, d);
				}  

			/* Set overall highest-scoring feature, if necessary */
			if (trellis->cell[state][i].submax->elem[di].score > trellis->cell[state][i].score)
				{
					trellis->cell[state][i].score = 
						trellis->cell[state][i].submax->elem[di].score;
					trellis->cell[state][i].trace = 
						trellis->cell[state][i].submax->elem[di].from_state;
					trellis->cell[state][i].length = 
						trellis->cell[state][i].submax->elem[di].end
						- trellis->cell[state][i].submax->elem[di].start + 1;
				}

			di++;
		}
}

void zSubInternalTrans (zTrellis *trellis, int from_state, int state, coor_t i, 
						int di, int subtrace, int sublen, score_t subscore, zDistribution *d)
{
	int ps;
	zStrIdx name;
	score_t tscore, len_score, scan_score, score;
	zSfeature maxf;

	int flag = 1; /*  1 stands for estseq structure consistent    */
	/* -1 stands for estseq structure inconsistent  */
	int do_return = 0;

	const char *state_name = zStrIdx2Char(trellis->hmm->state[state].name);

	tscore = zGetTransitionScore(trellis->hmm, from_state, state,trellis->tiso_group);
	ps = ('-' != trellis->hmm->state[state].strand);
	name = trellis->hmm->state[state].model;

	/* Consider boundary (min. length in range) separately */
	maxf.end = i;
	maxf.start = i - (d->start) + 1;
	maxf.from_state = from_state;

	len_score = zScoreDistribution(d, d->start);
	scan_score = zScoreScanners(trellis, name, state_name, maxf.start,
								maxf.end, ps, &do_return, &flag);

	maxf.score = trellis->cell[from_state][i-(d->start)].score
		+ tscore + len_score + scan_score;

	maxf.intrinsic_score = scan_score;

	if (do_return) return; /* do_return == TRUE means the feature */
	/* is inconsistent with EST sequence   */
	/* (see zScoreScanners function)       */

	/* Now deal with the rest of length values in range       */

	score = MIN_SCORE; 

	if (from_state == subtrace)
		{
			len_score = zScoreDistribution(d, (sublen+1));
			scan_score = zScoreScanners(trellis, name, state_name, i, i, ps,  
										&do_return, &flag);

			tscore = zGetTransitionScore(trellis->hmm, subtrace, state,trellis->tiso_group);
			score = trellis->cell[subtrace][i-sublen-1].score 
				+ subscore + len_score + scan_score + tscore;

			if (do_return) return; /* do_return == TRUE means the feature */
			/* is inconsistent with EST sequence   */
			/* (see zScoreScanners function)       */
		}

	/* Find highest-scoring feature in this range at position i */
	/* and store it in trellis->cell[state][i].submax->elem[di]  */

	if (score > maxf.score)
		{
			maxf.start = i - sublen; /* = i - (sublen+1) + 1 */
			maxf.score = score;
			maxf.intrinsic_score = subscore + scan_score;
		}

	if (maxf.score > trellis->cell[state][i].submax->elem[di].score)
		{
			trellis->cell[state][i].submax->elem[di].score = maxf.score;
			trellis->cell[state][i].submax->elem[di].intrinsic_score = maxf.intrinsic_score;
			trellis->cell[state][i].submax->elem[di].start = maxf.start;
			trellis->cell[state][i].submax->elem[di].end = maxf.end;
			trellis->cell[state][i].submax->elem[di].from_state = maxf.from_state;
		}
}


void zSubExplicitTrans (zTrellis *trellis, int from_state, int state, coor_t i,
                        int di, coor_t min, coor_t max, zDistribution *d)
{
	zStrIdx    name;
	zPhase_t   ext_phase, int_phase;
	score_t    tscore, scan_score, total_score, best_score;
	coor_t     best_length = 0, dna_len, length, prev_start;
	coor_t     start; /* start coord of external feature that preceeds */
	/* the current explicit feature                  */
	int        extind;
            
	int flag = 1; /*  1 stands for estseq structure consistent    */
	/* -1 stands for estseq structure inconsistent  */
	int do_return;  
        
	const char *state_name = zStrIdx2Char(trellis->hmm->state[state].name);

	int        ps = ('-' != trellis->hmm->state[state].strand);
	zDNA       *dna = (ps) ? trellis->dna : trellis->rdna;
           
	dna_len   = dna->length;
	name      = trellis->hmm->state[state].model;
                        
	if (EXTERNAL == trellis->hmm->state[from_state].type) { /* check phases */
		int_phase = trellis->hmm->state[state].phase;
		ext_phase = trellis->hmm->state[from_state].phase;
		if (ext_phase != int_phase) return;
	}          

	extind = trellis->extpos[from_state].size - 1;
	if (extind < 0) return; /* Empty extpos array               */
	/* No external states to connect to */
               
	while ((trellis->extpos[from_state].elem[extind] > ((int) max))
		   && (extind >= 0))
		{
            extind--;
		}
            
	start = trellis->extpos[from_state].elem[extind];
	prev_start = i;
	scan_score = 0.;
	do_return = 0;
                                  
	best_score  = MIN_SCORE;
	tscore  = zGetTransitionScore(trellis->hmm,from_state,state,trellis->tiso_group);
	while ((start >= min) && (extind >= 0))
		{
            length = i - start;
            total_score = tscore
				+ zScoreDistribution(d, length)
				+ trellis->cell[from_state][start].score;
			
            scan_score += zScoreScanners(trellis, name, state_name, (start+1),
                                         prev_start, ps, &do_return, &flag);
                   
			/* do_return == TRUE means the feature */
			/* is inconsistent with EST sequence   */
			/* (see zScoreScanners function)       */
            if (do_return)
				{
					extind = -1; /* force exit from loop */
				}
            else
				{
					total_score += scan_score;
            
					if (total_score >= best_score) {
						best_score = total_score;
						best_length = length;
					}
                                  
					extind--;
					prev_start = start;

					/* This is to avoid using negative index (-1) which */
					/* can happen in some cases at the last iteration.  */
					/* Insure++ complains if you don't do this!         */
					if (extind < 0)
						{
							start = 0;
						}
					else
						{
							start = trellis->extpos[from_state].elem[extind];
						}
				}
		}
               
	if (best_score > trellis->cell[state][i].submax->elem[di].score)
		{
            trellis->cell[state][i].submax->elem[di].score = best_score;
            trellis->cell[state][i].submax->elem[di].start = i - best_length + 1;
            trellis->cell[state][i].submax->elem[di].end = i;
            trellis->cell[state][i].submax->elem[di].from_state = from_state;
		}
}

/* this function scores a single base of an */
score_t zScoreInternalState(zTrellis *trellis, int state, 
							coor_t pos, bool first_base){
	int        ps; /* positive strand */
	int        do_return,flag;
	zStrIdx    name;
	zScanner*  scanner;
	score_t    total_score;
	coor_t     scored_coor;
	const char *state_name = zStrIdx2Char(trellis->hmm->state[state].name);

	ps = ('-' != trellis->hmm->state[state].strand);
	name = trellis->hmm->state[state].model;
	scanner = (ps) ? trellis->scanner[name] : trellis->rscanner[name];
	total_score = 0;
	scored_coor = (ps) ? pos : trellis->dna->length-pos-1;

	/* score content models */
	do_return = 0;
	flag = 0;
	
	total_score = zScoreScanners(trellis, name, state_name, pos, pos,
								 ps, &do_return, &flag);
	if(do_return){
		/* do_return == TRUE means the feature is inconsistent with EST sequence */
		return MIN_SCORE;
	}
	
	/* Give the duration score the correct geometric dist., if starting a new duration */
	if (first_base) {
		zDurationGroup *dg = 
			trellis->hmm->dmap[trellis->hmm->state[state].duration];
		total_score += 
			zScoreDurationGroup(dg, 1, zGetDNAGC(trellis->dna));
		/*EVAN if this changes to windowed gc we will get different results for 
		  running viterbi in the forward direction and in the backward direction 
		  because the duration score will be applied and the lower or higher coor
		  which may have different gc content*/
	}
	else{
		if(flag == -1){
			total_score += 2000;
		}
	}
	return total_score;
}


void zInternalTransHelper (zTrellis *trellis, int from_state, 
						   int state, coor_t i, score_t trans) {
	int        ps; /* positive strand */
	score_t    total_score;

	ps = ('-' != trellis->hmm->state[state].strand);
	
	/* check phase */
	if (EXTERNAL == trellis->hmm->state[from_state].type) {
		zPhase_t int_phase, ext_phase;
		int_phase = trellis->hmm->state[state].phase;
		ext_phase = trellis->hmm->state[from_state].phase;
		if (ext_phase != int_phase) return;
	}
	
	total_score = zScoreInternalState(trellis,state,i, 
									(trellis->cell[state][i-1].length == 0 ||
									 from_state != state));
	total_score += trans;
	
	if (trellis->forward) {
		trellis->forward[state][i] = 
			zFloatwiseScoreAdd(total_score + trellis->forward[state][i-1], 
							   trellis->forward[state][i]);
	} 
	if (trellis->cell) {
		total_score += trellis->cell[from_state][i-1].score;
		if (total_score > trellis->cell[state][i].score) {
			trellis->cell[state][i].score = total_score;
			if (from_state == state) {
				trellis->cell[state][i].length = 
					trellis->cell[state][i-1].length + 1;
				trellis->cell[state][i].trace = trellis->cell[state][i-1].trace;
			} else {
				trellis->cell[state][i].trace = from_state;
				trellis->cell[state][i].length = 1;
			}
		}
	}
}

/*EVAN this ignores the actual transition score for the first base */
static void zGInternalTrans (zTrellis *trellis, int from_state, 
							 int state, coor_t i) {
	zInternalTransHelper(trellis, from_state, state, i, zGetIntergenicContinueScore(trellis) );
}

static void zInternalTrans (zTrellis *trellis, int from_state, 
							int state, coor_t i) {
	zInternalTransHelper(trellis, from_state, state, i, zGetTransitionScore(trellis->hmm,from_state,state,trellis->tiso_group));
}

/* Score an external state.  pos is the last base of the state */
score_t zScoreExternalState(zTrellis *trellis, int ext_state, 
							coor_t pos, coor_t length, score_t cscore){
	/* pos - pos of higher coordinate state */
	int           ps;
	score_t       dscore, total_score;
	zHMM*         hmm = trellis->hmm;
	zDNA*         dna;

	pos = 0;/*EVAN shut up compiler*/ 

	ps = ('-' != trellis->hmm->state[ext_state].strand);
	dna = (ps) ? trellis->dna : trellis->rdna;
	
	dscore = zScoreDurationGroup(hmm->dmap[hmm->state[ext_state].duration], 
								 length, zGetDNAGC(dna));	
	
	total_score = cscore + dscore;
	/*EVAN removed for compatability with current code
	  - length*zGetIntergenicContinueScore(trellis);*/
	
	return total_score;
}

static void zExternalTrans (zTrellis *trellis, int int_state, int ext_state, coor_t pos){

	/* int_state is the state number of the internal state before the
	   external state that you are considering, ext_state is the state
	   number of the external state you are considering, pos is
	   (roughly) the high index (end on + strand begining on - strand)
	   of the external state being considered.  Back in the trellis
	   you get to this call by looking at each position in the
	   sequence, for each poistion looking at each external state, and
	   for each external state looking at each internal state that can
	   transistion into it as specified in the jmap vector for that
	   external state. */

	/*	    int_state          ext_state        
			-----------------|||||||||||||||||||||------------
			intron               exon           intron
			->                   ->
			transition1          transition2
	
			...[......content......]...
			<------duration----->
	*/

	int           length, pre_state;
	score_t       cscore, pre_score, total_score;
	zPhase_t      pre_phase, ext_phase, left_ext_phase;
	zSfeature    *f;
	zHMM         *hmm = trellis->hmm;
	zSFList      *sfl;
	char          c;
	int           flag = 1;
	const char    *state_name = zStrIdx2Char(trellis->hmm->state[ext_state].name);
	bool          ps = ('-' != trellis->hmm->state[ext_state].strand);
	zDNA          *dna = (ps) ? trellis->dna : trellis->rdna;
	bool          estseq_enabled = trellis->estseq != NULL;

	if(pos > trellis->dna->length - PADDING){
		return;
	}
	
	if (ps){
		sfl = trellis->external[hmm->state[ext_state].model];
	}
	else{
		sfl = trellis->rexternal[hmm->state[ext_state].model];
	}

	f = zSFListMoveFirst(sfl);
	while(f != NULL){
		if(ps){
			ext_phase = zSFEndPhase(dna, f);
			left_ext_phase = zSFBeginPhase(dna, f);
		}
		else{
			ext_phase = zSFBeginPhase(dna, f);
			left_ext_phase = zSFEndPhase(dna, f);
		}
		if (ext_phase != trellis->hmm->state[ext_state].phase){
			f = zSFListMoveNext(sfl);
			continue;
		}
		
		pre_state = int_state;
		pre_phase = hmm->state[pre_state].phase;
		
		if (!zCompatPhases(ps, pre_phase, left_ext_phase)){
			f = zSFListMoveNext(sfl);
			continue;
		}
		
		if(estseq_enabled  && !trellis->est_para_mode) {
			c = zGetEstseqSeq(trellis->estseq,f->end); 
			if (c == '1' && state_name[0] != 'E'){
				/* polyA or promoter */
				flag = - 1;
				f = zSFListMoveNext(sfl);
				continue;
			}
		}
		
		/*fprintf(stderr,"%d -> %d %s (%d%c) %f\n",f->start,f->end,zStrIdx2Char(hmm->state[ext_state].model),hmm->state[ext_state].model,hmm->state[ext_state].strand,f->score);*/
		
		length = f->end - f->start + 1;
		cscore = f->score;

		total_score = zScoreExternalState(trellis,ext_state,pos,length,cscore);
		total_score += zGetTransitionScore(trellis->hmm, pre_state, ext_state,trellis->tiso_group);

		if(zGetVerbosityLevel() > 5) {
			printf("%u ", f->start-PADDING+1);
			printf("%u ", f->end-PADDING+1);
			printf("%s ", zStrIdx2Char(f->name));
			printf("%c ", f->strand);
#ifdef DEBUG
			printf("%f ", f->content_score);
			printf("%f ", f->begin_score);
			printf("%f ", f->end_score);
#endif /* DEBUG */
			printf("%f ", f->cdscons_score);
			printf("%f ", total_score);				
			printf("%s ", zStrIdx2Char(trellis->hmm->reverse_somap[pre_state]));
			printf("\n");
		}
		
		pre_score = trellis->cell[pre_state][f->start-1].score;
		total_score += pre_score;
		
		if (total_score > trellis->cell[ext_state][f->end].score) {
			trellis->cell[ext_state][f->end].score = total_score;
			trellis->cell[ext_state][f->end].length = length;
			trellis->cell[ext_state][f->end].trace = pre_state;
			trellis->cell[ext_state][f->end].frag_frame =
                           zFragFrame2Char(f->lfrag, f->rfrag, f->frame);
		}
		f = zSFListMoveNext(sfl);
	}
}

void zInternalTransBack(zTrellis *trellis, int from_state, 
						int to_state, coor_t i) {
    zScanner* scanner;
	zStrIdx   feature;
	score_t   total_score;
	int       ps; /* positive strand */

	feature = trellis->hmm->state[from_state].model;
	ps = ('-' != trellis->hmm->state[from_state].strand);

	scanner = (ps) ? trellis->scanner[feature] :
					 trellis->rscanner[feature];

	total_score = scanner->score(scanner, (ps) ? i : trellis->dna->length-i-1)
	  + zGetTransitionScore(trellis->hmm, from_state, to_state,trellis->tiso_group)
	  + trellis->backward[to_state][i+1];

	trellis->backward[from_state][i] = 
	  zFloatwiseScoreAdd(total_score, trellis->backward[from_state][i]);
}

void zExternalTransBack (zTrellis *trellis, int int_state, 
							int ext_state, coor_t pos) {
	int           j, length, pre_state, ps;
	score_t       cscore, dscore, t1score, t2score,
	              gscore, pre_score, total_score;
	zPhase_t      int_phase, pre_phase;
	zSfeature    *f;
	zHMM         *hmm = trellis->hmm;
	zDNA         *dna = trellis->dna;
	zDNA         *rdna = trellis->rdna;
	zSFVec       *sfv;

/*	    int_state          ext_state        pre_state
	-----------------|||||||||||||||||||||------------
	     intron               exon           intron
	                ->                   ->
	            transition1          transition2
	
	              ...[......content......]...
	                 <------duration----->
*/
	int_phase = hmm->state[int_state].phase;
	t1score   = zGetTransitionScore(trellis->hmm, int_state, ext_state, trellis->tiso_group);
	if (MIN_SCORE == t1score) {return;}

	sfv = &trellis->fexternal[ext_state][pos+1];
	for (j = 0; j < sfv->size; j++) {
	  f = &sfv->elem[j];

	  ps = ('-' != f->strand);
	  if ((ps) ? !zCompatibleJumpToExon(dna, int_phase, f) :
		         !zCompatibleJumpFromExon(rdna, int_phase, f)) continue;

	  for (pre_state = 0; pre_state < hmm->states; pre_state++) {
		t2score = zGetTransitionScore(trellis->hmm, ext_state, pre_state,trellis->tiso_group);
		if (MIN_SCORE == t2score) continue;

		pre_phase = hmm->state[pre_state].phase;
		if ((ps) ? !zCompatibleJumpFromExon(dna, pre_phase, f) :
			       !zCompatibleJumpToExon(rdna, pre_phase, f)) continue;

		length     = f->end = f->start + 1;
		cscore     = f->score;
		dscore     = zScoreDurationGroup(hmm->dmap[hmm->state[ext_state].duration], length, zGetDNAGC(dna));
		gscore     = zGetTransitionScore(trellis->hmm, int_state, int_state,trellis->tiso_group);
		pre_score   = trellis->backward[pre_state][f->end+1];
		total_score = cscore + dscore + t1score + t2score 
			+ gscore + pre_score;

		trellis->backward[int_state][pos] = 
		  zFloatwiseScoreAdd(total_score, trellis->backward[int_state][pos]);
	  }
	}
}

#define STATE_TYPES 4
static zTransFunc zTransLookup[STATE_TYPES] = {
	zInternalTrans, /* INTERNAL */
	zGInternalTrans, /* GINTERNAL */
	zExternalTrans, /* EXTERNAL */
	zExplicitTrans, /* EXPLICIT */
};

zTransFunc zGetTransFunc(int type) {
 	if (type >= STATE_TYPES || type < 0) {
		zWarn("Attempt to get invalid state type");
		return NULL;
	}

	return zTransLookup[type]; 
}

static zTransFunc zBackTransLookup[STATE_TYPES] = {
	zInternalTransBack, /* INTERNAL */
	zInternalTransBack, /* GINTERNAL */
	zExternalTransBack, /* EXTERNAL */
	zExternalTransBack,  /* SHUTTLE */
};

zTransFunc zGetBackTransFunc(int type) {
	if (type >= STATE_TYPES || type < 0) {
		zWarn("Attempt to get invalid state type");
		return NULL;
	}
	return zBackTransLookup[type];
}

score_t zScoreScanners(zTrellis *trellis, zStrIdx name, const char *state_name, 
					   coor_t start, coor_t end, bool ps, int *do_return, int *flag){
	/* This function returns the sum of all scanner scores (seq, conseq, etc.) */
	/* for given region (start to end) and model type (name). It is called     */
	/* only for states of type EXPLICIT.                                       */
	/* At present (Aug. 18, 2004), there are no models for introns. If such    */
	/* models are developed later, or this function is called for a state      */
	/* other than intron, double-check that the scanners are working properly. */
	/* Nov. 16, 2004: Now works for explicit UTR states, too.                  */

  zSfeature f;
  zScanner *scanner;
  zScanner *c_scanner;
  zScanner *nc_scanner;
  zScanner *e_scanner;
  zAlignmentScanner *a_scanner;
  zAlignmentScanner *na_scanner;
  coor_t dna_len, pos;
  score_t score;

  zSetScanners(&scanner, &c_scanner, &nc_scanner, &e_scanner, &a_scanner, 
               &na_scanner, trellis, name, ps);

  dna_len = trellis->dna->length;

  /* Required by the scoring functions  */
  /* See zScoreFeature in zScanner.c,   */
  /* zScanner.c and zEstScanner.c */
  f.name = name;
  f.lfrag = 0;
  f.rfrag = 0;
  f.frame = 0;

  f.start = (ps) ? start : (dna_len-end-1);
  f.end = (ps) ? end : (dna_len-start-1);
  score = 0.;

  if (scanner != NULL) score += scanner->scoref(scanner, &f);

  /*printf("scan score %s %d->%d %f\n",state_name,start,end,score);*/
  if (c_scanner != NULL) score += c_scanner->scoref(c_scanner, &f);
  if (a_scanner != NULL)
	  {
		if (na_scanner == NULL)
          zDie("Could not find NCCons alignment scanner at zScoreScanners");
	  
      /* PSC is a scaling constant defined in zFeatureFactory.h */
      score += (PSC*(a_scanner->scoref(a_scanner, &f)
                     - na_scanner->scoref(na_scanner, &f)));
    }

  *do_return = 0;
  if (trellis->estseq != NULL){
      if (trellis->est_para_mode){
          /* EST parametric mode. Use EST scanner */
          if (e_scanner != NULL)  score += e_scanner->scoref(e_scanner, &f);
	  }
      else if(*flag != -1){
		  /* This is a mess but is backward compatable with Mikhails old explicit
			 state code*/
		  coor_t two_pos;
		  coor_t one_pos;
		  /* EST non-parametric mode. Enforce */
          /* consistency with EST sequence    */		  
		  two_pos = zGetEstseqPrevPosBounded(trellis->estsearch,'2',end,start-1);
		  /* A feature containing "I" can be non-intronic 
			 with an extremely low probability           */
		  if (two_pos >= start && strncmp(state_name, "Intron", 6) != 0){
			  pos = two_pos;
			  *flag = -1;
			  score -= 1000;
		  }
		  else{
			  pos = start;
		  }
		  one_pos = zGetEstseqPrevPosBounded(trellis->estsearch,'1',end,pos-1);
		  /* Only UTR can contain an "E"                  
			 Warning: if EXPLICIT and EXTERNAL states are 
			 ever merged in any way, make sure exons are  
			 excluded here as well!!!                     */
		  if (one_pos >= pos && strncmp(state_name, "Utr", 3) != 0){
			  *flag = -1;
			  *do_return = 1;
			  score = MIN_SCORE;
		  }
	  }
  }

  return score;
}

static bool set_nc_scanners = false;
static zScanner* ncs_plus = NULL;
static zScanner* ncs_minus = NULL;
static zAlignmentScanner* nas_plus = NULL;
static zAlignmentScanner* nas_minus = NULL;

void zSetScanners(zScanner **s, zScanner **cs, zScanner **ncs,
                  zScanner **es, zAlignmentScanner **as, zAlignmentScanner **nas, 
                  zTrellis *trellis, zStrIdx name, bool ps)
{
  /* Set scanner pointers                                          */
  /* This will assign a valid pointer only if both the scanner     */
  /* and its associated sequence are present. This is necessary    */
  /* in order to avoid dangling pointers and possible segm. faults */
  /* The only exception is the DNA sequence scanner - if DNA seq   */ 
  /* isn't there, you wouldn't even get to this point.             */

  *s = NULL;
  *cs = NULL;
  *ncs = NULL;
  *es = NULL;
  *as = NULL;
  *nas = NULL;

  if(set_nc_scanners == false){
	  if(trellis->conseqscanner != NULL){
		  ncs_plus = trellis->conseqscanner[zChar2StrIdx("INTRONCONS")];
		  ncs_minus = trellis->rconseqscanner[zChar2StrIdx("INTRONCONS")];
	  }
	  if(trellis->alignmentscanner != NULL){
		  nas_plus = trellis->alignmentscanner[zChar2StrIdx("NCCons")];
		  nas_minus = trellis->ralignmentscanner[zChar2StrIdx("NCCons")];
	  }
	  set_nc_scanners = true;
  }

  if (ps)
    {
      if (trellis->scanner != NULL)
        { 
          *s = trellis->scanner[name];
        }
      if ((trellis->conseq != NULL) && (trellis->conseqscanner != NULL))
        {
          *cs = trellis->conseqscanner[name];
          *ncs = ncs_plus;
        }
      if ((trellis->estseq != NULL) && (trellis->estseqscanner != NULL))
        {
          *es = trellis->estseqscanner[name];
        }
      if ((trellis->alignment != NULL) && (trellis->alignmentscanner != NULL))
        {
          *as = trellis->alignmentscanner[name];
          *nas = nas_plus;
        }
    }
  else
    {
      if (trellis->rscanner != NULL)
        { 
          *s = trellis->rscanner[name];
        }
      if ((trellis->conseq != NULL) && (trellis->rconseqscanner != NULL))
        {
          *cs = trellis->rconseqscanner[name];
          *ncs = ncs_minus;
        }
      if ((trellis->estseq != NULL) && (trellis->restseqscanner != NULL))
        {
          *es = trellis->restseqscanner[name];
        }
      if ((trellis->alignment != NULL) && (trellis->ralignmentscanner != NULL))
        {
          *as = trellis->ralignmentscanner[name];
          *nas = nas_minus;
        }
    } 
}

#endif
