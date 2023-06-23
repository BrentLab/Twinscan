/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
  zPairTrellis.c - part of the ZOE library for genomic analysis
  
  Copyright (C) 2001-2002 Ian F. Korf
\******************************************************************************/

#ifndef ZOE_PAIR_TRELLIS_C
#define ZOE_PAIR_TRELLIS_C

#include <assert.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "zPairTrellis.h"
#include "zHardCoding.h"

static const score_t MIN_INIT_SCORE = -10000;
extern coor_t any_new_cpoint;

void zShowPairTrellisCell(zPairTrellis *trellis, coor_t i, coor_t j, int k); 
static void zQuickTracePartialTrellis(zPairTrellis *trellis, coor_t g_begin, coor_t g_end, coor_t c_begin, coor_t c_end); 
static void zTracePartialPairTrellis(zPairTrellis *trellis, coor_t g_begin, coor_t g_end, coor_t c_begin, coor_t c_end, int state, zAFVec *afv); 

/*********************************************\
 Regular Viterbi Variables
\*********************************************/

static void zAllocPartialViterbiVars(zPairTrellis *trellis, coor_t genomic_start, coor_t genomic_end, coor_t cdna_start, coor_t cdna_end) {
	coor_t i, j;
	coor_t length = (cdna_end - cdna_start + 1);
	size_t cell_size = sizeof(zPairTrellisCell)*trellis->hmm->states;
	size_t row_size = length*sizeof(zPairTrellisCell*);

	for (i = genomic_start; i <= genomic_end; i++) {
		trellis->cell[i] = (zPairTrellisCell**) zSliceAlloc(row_size, "zAllocViterbi cells[i]");
		trellis->cell[i] -= cdna_start;
		for (j = cdna_start; j <= cdna_end; j++) {
			int k;
			char strtemp[64];
			sprintf(strtemp, "zAllocViterbi cells[%u][%u]", i, j);
			trellis->cell[i][j] = (zPairTrellisCell*) zSliceAlloc(cell_size, strtemp);
			for (k = 0; k < trellis->hmm->states; k++) {
				trellis->cell[i][j][k].length = -1;
				trellis->cell[i][j][k].trace  = -1;
				trellis->cell[i][j][k].keep   = 0;
				trellis->cell[i][j][k].score  = MIN_SCORE;
			}
		}
	}
}

static void zFreePartialViterbiVars(zPairTrellis *trellis, coor_t genomic_start, coor_t genomic_end, coor_t cdna_start, coor_t cdna_end) {
	coor_t i, j;
	coor_t length = (cdna_end - cdna_start + 1);
	size_t cell_size = sizeof(zPairTrellisCell)*trellis->hmm->states;
	size_t row_size = length*sizeof(zPairTrellisCell*);

	for (i = genomic_start; i <= genomic_end; i++) {
		trellis->cell[i] += cdna_start;
		if (trellis->cell[i] != NULL) {
			for (j = 0; j < length; j++) {
				zSliceFree(cell_size, trellis->cell[i][j]);
			}
			zSliceFree(row_size, trellis->cell[i]);
		}
	}
}

static void zAllocPartialVars(zPairTrellis* trellis, score_t ***vars, coor_t genomic_start, coor_t genomic_end, coor_t cdna_start, coor_t cdna_end) {
	coor_t i, j;
	coor_t length = (cdna_end - cdna_start + 1);
	size_t cell_size = sizeof(score_t)*trellis->hmm->states;
	size_t row_size = length*sizeof(score_t*);

	for (i = genomic_start; i <= genomic_end; i++) {
		vars[i] = (score_t**) zSliceAlloc(row_size, "zAllocViterbi vars[i]");
		vars[i] -= cdna_start;
		for (j = cdna_start; j <= cdna_end; j++) {
			int k;
			char strtemp[64];
			sprintf(strtemp, "zAllocViterbi cells[%u][%u]", i, j);
			vars[i][j] = (score_t*) zSliceAlloc(cell_size, strtemp);
			for (k = 0; k < trellis->hmm->states; k++) {
				vars[i][j][k]     = MIN_SCORE;
			}
		}
	}
}

static void zFreePartialVars(zPairTrellis *trellis, score_t ***vars, coor_t genomic_start, coor_t genomic_end, coor_t cdna_start, coor_t cdna_end) {
	coor_t i, j;
	coor_t length = (cdna_end - cdna_start + 1);
	size_t cell_size = sizeof(score_t)*trellis->hmm->states;
	size_t row_size = length*sizeof(score_t*);

	for (i = genomic_start; i <= genomic_end; i++) {
		vars[i] += cdna_start;
		if (vars[i] != NULL) {
			for (j = 0; j < length; j++) {
				zSliceFree(cell_size, vars[i][j]);
			}
			zSliceFree(row_size, vars[i]);
		}
	}
}

static void zAllocViterbiVars(zPairTrellis *trellis) {
	int       i;
	if (NULL != trellis->cell) return;

	zTrace2("initializing Viterbi vars");

	trellis->cell = zCalloc(trellis->genomic->length, sizeof(zPairTrellisCell**), "zAllocViterbi cell");

	/* Allocating backward link arrays for all states                 */
	/* But only those referring to external states will be used later */
	trellis->extpos = zMalloc((trellis->hmm->states*sizeof(zIVec)), "zAllocViterbi extpos");

	for (i = 0; i < trellis->hmm->states; i++) {
		zInitIVec(&trellis->extpos[i], 1);
	}
}

static void zAllocForwardVars(zPairTrellis *trellis) {
	zTrace2("initializing forward vars");
	if (NULL != trellis->forward) return;
	trellis->forward = zCalloc(trellis->genomic->length, sizeof(score_t**), "zAllocForwardVars forward");
}

static void zAllocBackwardVars(zPairTrellis *trellis) {
	zTrace2("initializing backward vars");
	if (NULL != trellis->backward) return;
	trellis->backward = zCalloc(trellis->genomic->length, sizeof(score_t**), "zAllocBackwardVars backward");
}

static void zFreeViterbiVars(zPairTrellis *trellis) {
	int p;

	if (NULL == trellis->cell) return;

	for (p = 0; p < trellis->hmm->states; p++) {
		zFreeIVec(&trellis->extpos[p]);
	}
	zFree(trellis->extpos);

	for (p = 0; p < trellis->mem_blocks->hsps; p++) {
		zHSP *r = &trellis->mem_blocks->hsp[p];
		zFreePartialViterbiVars(trellis, r->g_start, r->g_end, r->c_start, r->c_end);
	}

	zFree(trellis->cell);
	trellis->cell = NULL;
}

static void zFreeForwardVars(zPairTrellis *trellis) {
	int p;

	if (NULL == trellis->forward) return;

	for (p = 0; p < trellis->mem_blocks->hsps; p++) {
		zHSP *r = &trellis->mem_blocks->hsp[p];
		zFreePartialVars(trellis, trellis->forward, r->g_start, r->g_end, r->c_start, r->c_end);
	}
	zFree(trellis->forward);
	trellis->forward = NULL;
}

static void zFreeBackwardVars(zPairTrellis *trellis) {
	int p;

	if (NULL == trellis->backward) return;

	for (p = 0; p < trellis->mem_blocks->hsps; p++) {
		zHSP *r = &trellis->mem_blocks->hsp[p];
		zFreePartialVars(trellis, trellis->backward, r->g_start, r->g_end, r->c_start, r->c_end);
	}
	zFree(trellis->backward);
	trellis->backward = NULL;
}

static void zGarbageCollectPartialTrellis(zPairTrellis *trellis, coor_t genomic, coor_t cdna) {
	coor_t i, j;
	int k, limit;
	zSeedAlignment *mem_blocks = trellis->mem_blocks;
	size_t cell_size = sizeof(zPairTrellisCell)*trellis->hmm->states;
	zTrace2("Garbage Collecting (%u, %u)", genomic, cdna);

	for (k = 0; k < mem_blocks->hsps; k++) {
		zHSP* r = &mem_blocks->hsp[k];
		if (genomic <= r->g_end && cdna <= r->c_end) {
			break;
		}
	}
	limit = k;

	for (k = 0; k < limit; k++) {
		zHSP* r = &mem_blocks->hsp[k];
		for (i = r->g_start; i <= r->g_end; i++) {
			for (j = r->c_start; j <= r->c_end; j++) {
				if(trellis->cell[i][j] != NULL && trellis->cell[i][j][0].keep == 0) {
					zSliceFree(cell_size, trellis->cell[i][j]);
					trellis->cell[i][j] = NULL;
					if (trellis->forward != NULL) {
						zSliceFree(sizeof(score_t)*trellis->hmm->states, trellis->forward[i][j]);
						trellis->forward[i][j] = NULL;
					}
				}
			}
		}
	}
}

static void zCheckViterbiVariables(zPairTrellis *trellis, coor_t gpos, coor_t cpos) {
	int k;
	int required_block = zGetAlignmentBlock(trellis->mem_blocks, gpos, cpos);
	if (trellis->allocated_mem_blocks >= required_block) return; /* Thumbs up! */

	for (k = trellis->allocated_mem_blocks + 1; k <= required_block; k++) {
		zHSP* r = &trellis->mem_blocks->hsp[k];
		zAllocPartialViterbiVars(trellis, r->g_start, r->g_end, r->c_start, r->c_end);
	}
	trellis->allocated_mem_blocks = required_block;
}

static void zCheckForwardVariables(zPairTrellis *trellis, coor_t gpos, coor_t cpos) {
	int k;
	int required_block = zGetAlignmentBlock(trellis->mem_blocks, gpos, cpos);
	if (trellis->forward == NULL) return;
	if (trellis->allocated_fwd_blocks >= required_block) return; /* Thumbs up! */
	for (k = trellis->allocated_fwd_blocks + 1; k <= required_block; k++) {
		zHSP* r = &trellis->mem_blocks->hsp[k];
		zAllocPartialVars(trellis, trellis->forward, r->g_start, r->g_end, r->c_start, r->c_end);
	}
	trellis->allocated_fwd_blocks = required_block;
}

static void zCheckBackwardVariables(zPairTrellis *trellis, coor_t gpos, coor_t cpos) {
	int k;
	int required_block = zGetAlignmentBlock(trellis->mem_blocks, gpos, cpos);
	if (trellis->allocated_bak_blocks <= required_block) return; /* Thumbs up! */
	for (k = required_block; k < trellis->allocated_bak_blocks; k++) {
		zHSP* r = &trellis->mem_blocks->hsp[k];
		zAllocPartialVars(trellis, trellis->backward, r->g_start, r->g_end, r->c_start, r->c_end);
	}
	trellis->allocated_bak_blocks = required_block;
}

/*********************************************\
 Regular Viterbi Factories
\*********************************************/

static void zAllocFactories(zPairTrellis *trellis) {
	zFeatureFactory **ffactory;
	zHMM_State        *state;
	zIVec           **fmap5,**fmap3,*jumps;

	zHMM    *hmm = trellis->hmm;
	int      i,j,k,insert;

	/* create factories for non-internal states */
	trellis->factory = zCalloc(hmm->feature_count, sizeof(zFeatureFactory*),
							"zAllocFactories: factory[]");
	trellis->fmap5   = zCalloc(hmm->feature_count, sizeof(zIVec*),
							"zAllocFactories: fmap5[]");
	trellis->fmap3   = zCalloc(hmm->feature_count, sizeof(zIVec*),
							"zAllocFactories: fmap3[]");
	for (i = 0; i <  hmm->states; i++) {
		state = &hmm->state[i];
		if (NULL == state->ffactory) continue;

		ffactory = &trellis->factory[state->model];
		fmap5 = &trellis->fmap5[state->model];
		fmap3 = &trellis->fmap3[state->model];
		if (NULL != *ffactory){
			zPushIVec(*fmap3,i);
			jumps = trellis->hmm->jmap[i];
			for(j = 0;j < jumps->size; j++){
				insert = 1;
				for(k = 0;k < (*fmap5)->size; k++){
					if((*fmap5)->elem[k] == jumps->elem[j]){
						insert = 0;
						break;
					}
				}
				if(insert == 1){
					zPushIVec(*fmap5,jumps->elem[j]);
				}
			}
			continue;
		}
		*fmap5 = zMalloc(sizeof(zIVec),"zAllocFacotries: fmap5");
		*fmap3 = zMalloc(sizeof(zIVec),"zAllocFacotries: fmap3");
		zInitIVec(*fmap5,3);
		zInitIVec(*fmap3,3);
		zPushIVec(*fmap3,i);
		jumps = trellis->hmm->jmap[i];
		for(j = 0;j < jumps->size; j++){
			zPushIVec(*fmap5,jumps->elem[j]);
		}
		*ffactory = zMalloc(sizeof(zFeatureFactory),"zAllocFac: factory[i]");
		state->ffactory(*ffactory,
						trellis->genomic, 
						NULL,
						trellis->scanner, 
						NULL,
						NULL,
						NULL,
						hmm->feature_count, 
						state->strand, 
						false,
						false,
						false,
						NULL,
						NULL
						);
	}
	/* allcate external arrays */
	trellis->external = zMalloc(hmm->feature_count*sizeof(zSFList*),
								"zAllocFactories: external[]");
	for (i = 0; i <  hmm->feature_count; i++) {
		if(trellis->factory[i] != NULL){
			trellis->external[i] = zMalloc(sizeof(zSFList),"zAllocFactories: external[i]");
			zInitSFList(trellis->external[i]);
		}
	}
}

static void zFreeFactories(zPairTrellis *trellis) {
	int i;

	for (i = 0; i < trellis->hmm->feature_count; i++) {
		if (NULL != trellis->factory[i]) {
			zFreeFeatureFactory(trellis->factory[i]);
			zFree(trellis->factory[i]);
			zFreeSFList(trellis->external[i]);
			zFree(trellis->external[i]);
			zFreeIVec(trellis->fmap5[i]);
			zFreeIVec(trellis->fmap3[i]);
			zFree(trellis->fmap5[i]);
			zFree(trellis->fmap3[i]);
		}
	}
	zFree(trellis->fmap5);
	zFree(trellis->fmap3);
	zFree(trellis->factory);
	zFree(trellis->external);
	trellis->external = NULL;
	trellis->factory = NULL;
}

/*********************************************\
 zPairTrellis Utilities
\*********************************************/

void zInitPairTrellis (zPairTrellis *trellis, zSeedAlignment *seed, zDNA *genomic, zDNA *cdna, zHMM *hmm) {
	int         i;

	/* clear out pointers */
	trellis->cell      = NULL;
	trellis->forward   = NULL;
	trellis->backward  = NULL;
	trellis->genomic   = NULL;
	trellis->cdna      = NULL;
	trellis->fcdna     = NULL;
	trellis->rcdna     = NULL;
	trellis->hmm       = NULL;
	trellis->fexternal = NULL;
	trellis->scanner   = NULL;
	trellis->factory   = NULL;
	trellis->extpos    = NULL;
	trellis->external  = NULL;
	trellis->seed      = NULL;
	trellis->blocks    = NULL;
	trellis->mem_blocks= NULL;
	
	trellis->padding   = PADDING;

	/* initial setup */
	trellis->fcdna   = zMalloc(sizeof(zDNA), "zInitPairTrellis dna");
	trellis->rcdna   = zMalloc(sizeof(zDNA), "zInitPairTrellis rdna");
	trellis->genomic = zMalloc(sizeof(zDNA), "zInitPairTrellis rdna");
	zInitDNA(trellis->fcdna);
	zInitDNA(trellis->rcdna);
	zInitDNA(trellis->genomic);

	zCopyDNA(genomic, trellis->genomic);
	zSetDNAPadding(trellis->genomic,trellis->padding);

	zCopyDNA(cdna, trellis->fcdna);
	zSetDNAPadding(trellis->fcdna,trellis->padding);
	zCopyDNA(trellis->fcdna, trellis->rcdna);
	zAntiDNA(trellis->rcdna);
	trellis->cdna = trellis->fcdna;

	/* Seed alignment */

	if (seed != NULL) {
		trellis->seed       = zMalloc(sizeof(zSeedAlignment), "zInitPairTrellis: trellis->seed");
		zInitSeedAlignment(trellis->seed);
		zCopySeedAlignment(seed, trellis->seed); 
	} else {
		trellis->seed = NULL;
	}

	trellis->blocks     = zMalloc(sizeof(zSeedAlignment), "zInitPairTrellis: trellis->blocks");
	trellis->mem_blocks = zMalloc(sizeof(zSeedAlignment), "zInitPairTrellis: trellis->mem_blocks");
	if (trellis->seed == NULL) {
		zHSP *hsp;
		zInitSeedAlignment(trellis->blocks);
		trellis->blocks->strand = '+';
		trellis->blocks->def = zMalloc(sizeof(char)*16, "zInitPairTrellis: trellis->blocks->def");
		strcpy(trellis->blocks->def, ">dummy blocks");
		trellis->blocks->hsps = 1;
		trellis->blocks->hsp = zMalloc(sizeof(zHSP), "zInitPairTrellis: trellis->blocks->hsp");
		trellis->blocks->gb_start = 1;
		trellis->blocks->gb_end   = genomic->length;

		hsp = &trellis->blocks->hsp[0];
		hsp->g_start = 1; 
		hsp->c_start = 1; 
		hsp->g_end = genomic->length;
		hsp->c_end = cdna->length;
		zTranslateSeedAlignment(trellis->blocks, trellis->padding - 1);
	} else {
		zTranslateSeedAlignment(trellis->seed, trellis->padding - 1);
		zPruneSeedAlignment(trellis->genomic, trellis->cdna, trellis->seed, 20);

		/* Alignment blocks for stepping stone */

		zInitSeedAlignment(trellis->blocks);
		zSeedAlignment2AlignmentBlocks(trellis->genomic, trellis->cdna, trellis->seed, trellis->blocks);
	}

	/* Nonoverlapping memory blocks for memory optimized trellis */

	zInitSeedAlignment(trellis->mem_blocks);
	zAlignmentBlocks2MemoryBlocks(trellis->blocks, trellis->mem_blocks);

	/* Set the blocks and mem_blocks end points to cover the whole sequence */
	/* gb_start and gb_end should only be read from trellis->blocks. trellis->seed could very well be NULL, and mem_blocks can be anything */

	trellis->mem_blocks->hsp[0].c_start = trellis->padding - 1;
	trellis->mem_blocks->hsp[0].g_start = MAX(1, trellis->blocks->gb_start)-1;
	trellis->mem_blocks->hsp[trellis->mem_blocks->hsps-1].g_end = MIN(trellis->genomic->length - 1, trellis->blocks->gb_end);
	trellis->mem_blocks->hsp[trellis->mem_blocks->hsps-1].c_end = trellis->cdna->length - 1;

	trellis->blocks->hsp[0].c_start = trellis->padding - 1;
	trellis->blocks->hsp[0].g_start = MAX(1, trellis->blocks->gb_start)-1;
	trellis->blocks->hsp[trellis->blocks->hsps-1].g_end = MIN(trellis->genomic->length - 1, trellis->blocks->gb_end);

	trellis->allocated_mem_blocks = -1;
	trellis->allocated_fwd_blocks = -1;
	trellis->allocated_bak_blocks = trellis->mem_blocks->hsps;

	trellis->hmm  = hmm;

	/* isocore group */
	trellis->tiso_group = -1;
	if(trellis->hmm->iso_transitions > 0) {
		for(i=0;i<trellis->hmm->iso_transitions;i++){
			if(zGetDNAGC(trellis->genomic) < trellis->hmm->iso_transition[i]){
				trellis->tiso_group = i;
				break;
			}
		}
	}
	trellis->iiso_group = -1;
	if(hmm->iso_states > 0){
		for(i = 0; i < hmm->iso_states; i++){
			if(zGetDNAGC(trellis->genomic) < hmm->iso_state[i]){
				trellis->iiso_group=i;
				break;
			}
		}
	}	

	/* create scanners */
	trellis->scanner  = zCalloc(hmm->feature_count, sizeof(zScanner*), 
								"zInitPairTrellis: scanner[]");
	
	for (i = 0; i < hmm->feature_count; i++) {
		zModel* model = hmm->mmap[i];
		if (model == NULL) continue;
		if (model->seq_type == DNA) {
			trellis->scanner[i] = zMalloc(sizeof(zScanner), "zInitPairTrellis scan");
			zInitScanner(trellis->scanner[i], trellis->cdna->seq, model);
		} else if (model->seq_type == GENOMIC) {
			trellis->scanner[i] = zMalloc(sizeof(zScanner), "zInitPairTrellis scan");
			zInitScanner(trellis->scanner[i], trellis->genomic->seq, model);
		} else if (model->seq_type == PAIR) {
			trellis->scanner[i] = zMalloc(sizeof(zScanner), "zInitPairTrellis scan");
			zInitScanner(trellis->scanner[i], trellis->genomic->seq, model);
		} else {
			zDie("non-DNA model in sequence models\n");
		}
	}

	/********************************************************************************
	 * Set the scanner boundaries. This is the correct way of doing it, considering *
	 * PADDING in the beginning and adjusting the boundaries. The standard way does *
	 * not do it that way. But I did not want to change the main branch since it    *
	 * a lot of work to convince people that it is the right way to do. The outputs *
	 * change in a few cases in the benchmark test sets if you did this for iscan   *
	 ********************************************************************************/
	if (hmm->mode == GPAIRHMM) {
		for (i = 0; i < hmm->states; i++) {
			zModel* model = hmm->mmap[hmm->state[i].model];
			zScanner *scanner = trellis->scanner[hmm->state[i].model];
			if (model->seq_type == GENOMIC) {
				scanner->min_pos = trellis->blocks->gb_start - 1 + model->length;
			} else if (model->seq_type == DNA) {
				scanner->min_pos = trellis->padding          - 1 + model->length;
			} else if (model->seq_type == PAIR) {
				scanner->min_gpos = trellis->blocks->gb_start - 1 + model->length/2;
				scanner->min_pos  = trellis->padding          - 1 + model->length/2;
			}

			/* Only states with init_prob > 0 can extend NULL alignment. Others can only extend something that those states started */
			if (zGetInitProb(hmm, i, trellis->tiso_group) == MIN_SCORE) {
				scanner->min_pos  += 1;
				scanner->min_gpos += 1;
			}
			if (model->seq_type == GENOMIC || model->seq_type == DNA) {
				/* Pre-compute DNA scanners for EXPLICIT    */
				/* states as a "running sum". Required by   */
				/* the implementation of EXPLICIT states    */
				/* (see zScoreScanners in zTransition.c for */
				/* details).                                */
				
				zPreComputeScanner(scanner);
			}
		}
	}

	if (hmm->mode != GPAIRHMM) zAllocFactories(trellis);
}

void zFreePairTrellis (zPairTrellis *trellis) {
	int i;
	
	for (i = 0; i < trellis->hmm->feature_count; i++) {
		zScanner* scanner = trellis->scanner[i];
		if (scanner  != NULL) {
			if (scanner->model->seq_type == GENOMIC || scanner->model->seq_type == DNA) {
				zDePreComputeScanner(scanner);
			}
			zFreeScanner(scanner);
			zFree(scanner);
		}
	}

	zFree(trellis->scanner);        trellis->scanner = NULL; 

	if (trellis->hmm->mode != GPAIRHMM) zFreeFactories(trellis);

	zFreeViterbiVars(trellis);
	zFreeForwardVars(trellis);
	zFreeBackwardVars(trellis);

	/* Stepping Stone stuff */
	zFreeSeedAlignment(trellis->mem_blocks);
	zFreeSeedAlignment(trellis->blocks);
	zFreeSeedAlignment(trellis->seed);
	zFree(trellis->mem_blocks);
	zFree(trellis->blocks);
	zFree(trellis->seed);
	
	zFreeDNA(trellis->fcdna);
	zFree(trellis->fcdna);
	trellis->fcdna = NULL;
	zFreeDNA(trellis->rcdna);
	zFree(trellis->rcdna);
	trellis->rcdna = NULL;
	zFreeDNA(trellis->genomic);
	zFree(trellis->genomic);
	trellis->genomic = NULL;
}

/*********************************************\
 Regular Viterbi Decoding
\*********************************************/

/* Start the alignment at position (gpos, cpos) by calculating the initial probabilities
   and making the nodes with that probabilities. Also update the trellis->cell for that position */

static void zStartAlignmentForward(zPairTrellis *trellis, coor_t gpos, coor_t cpos) {
	int state;
	score_t   score;
	coor_t    real_gpos, real_cpos;
	zHMM     *hmm = trellis->hmm;
	for (state = 0; state < trellis->hmm->states; state++) {
		/* Start state */
		zPairTrellisCell *cell = zGetCurrentCell(trellis, gpos, cpos, state);
		score  = zGetFixedInitProb(hmm, state, trellis->iiso_group, trellis->genomic->gc); /* Fixed Initial Probability */
		if (score > MIN_SCORE) {
			cell->score  = score;
			cell->length = 0;
			cell->trace  = state;
			cell->keep   = 1;
		}

		/* First instance of state should not have transition prob in it. *
                 * If you pass it to Viterbi, it will add transition prob to it.  */
		real_gpos = gpos + zGetGenomicIncrement(hmm, state);
		real_cpos = cpos + zGetCDnaIncrement(hmm, state);
		cell = zGetCurrentCell(trellis, real_gpos, real_cpos, state);
		score += zGetScannerScore(trellis, trellis->scanner[trellis->hmm->state[state].model], state, real_gpos, real_cpos); /* Score that pos */
		if (trellis->forward != NULL) trellis->forward[real_gpos][real_cpos][state] = cell->score;
		if (score > MIN_SCORE) {
			cell->score  = score;
			cell->length = 1;
			cell->trace  = state;
			cell->keep   = 1;
		}
	}
}

/* Finish the alignment at position (gpos, cpos) by calculating the initial probabilities
   and adding that probabilities to the trellis->cell for that position */

static void zFinishAlignmentForward(zPairTrellis *trellis, coor_t gpos, coor_t cpos) {
	int state;
	trellis->forward_score = MIN_SCORE;
	for (state = 0; state < trellis->hmm->states; state++) {
		zPairTrellisCell *cell = zGetCurrentCell(trellis, gpos, cpos, state);
		score_t score  = zGetInitProb(trellis->hmm, state, trellis->iiso_group);
		cell->score  += score;
		if (trellis->forward != NULL) trellis->forward_score = zFloatwiseScoreAdd(trellis->forward[gpos][cpos][state]+score, trellis->forward_score);
	}
}

static void zRunPartialPairViterbiAndForward(zPairTrellis *trellis, coor_t gstart, coor_t gend, coor_t cstart, coor_t cend) {
	zHMM         *hmm = trellis->hmm;
	coor_t        genomic, cdna; /* iterators for sequence */
	int           state;         /* iterator for internal states */
	int           prev;          /* iterator for previous states */
	zIVec*        jumps;

	zPairTrellisCell *cell;

	/* induction */
	zTrace2("beginning induction");

	/* Make sure that all the required blocks are loaded */
	zCheckViterbiVariables(trellis, gend, cend);
	zCheckForwardVariables(trellis, gend, cend);
	zTrace2("Calling (%u, %u) (%u, %u)", gstart, cstart, gend, cend);
	for (genomic = gstart; genomic <= gend; genomic++) {
		for (cdna = cstart; cdna <= cend; cdna++) {
			for (state = 0; state < hmm->states; state++) {
				if (NULL == (cell = zGetCurrentCell(trellis, genomic, cdna, state))) continue;
				if (cell->score != MIN_SCORE) continue; /* Has been done already */

#ifdef FORWARD
				if (trellis->forward != NULL) trellis->forward[genomic][cdna][state] = MIN_SCORE;
#endif
				jumps = hmm->jmap[state];
				for (prev = 0; prev < jumps->size; prev++) {
					zGetPairTransFunc(hmm->state[state].type)
						(trellis, jumps->elem[prev], state, genomic, cdna);
				}
			}
		}
	}
}

zAFVec* zRunPairViterbiAndForward (zPairTrellis *trellis, score_t* path_score) {
	zDNA         *cdna = trellis->cdna;
	zDNA         *genomic = trellis->genomic;

	zAFVec       *afv;
	int           i;
	coor_t        gmin, cmin, gmax, cmax;

	gmin = MAX(                   trellis->padding - 1, trellis->blocks->gb_start - 1);
	cmin =                        trellis->padding - 1;
	gmax = MIN( genomic->length - trellis->padding - 1, trellis->blocks->gb_end);
	cmax =         cdna->length - trellis->padding - 1;
	 	 
	/* 	Viterbi and Forward Alg initialization */

	zAllocViterbiVars(trellis);
	zCheckViterbiVariables(trellis, gmin, cmin);

	if (0) zAllocForwardVars(trellis);
	zCheckForwardVariables(trellis, gmin, cmin);

	/* induction - first pass */

	zTrace2("beginning first induction");
	zStartAlignmentForward(trellis, gmin, cmin);

	zTrace2("running viterbi");
	for (i = 0; i < trellis->blocks->hsps; i++) {
		zHSP *hsp = &trellis->blocks->hsp[i];
		zRunPartialPairViterbiAndForward(trellis, hsp->g_start, hsp->g_end, hsp->c_start, hsp->c_end);
		/* Garbage collect in regions that wont be required by the next hsp */
		if (i < trellis->blocks->hsps - 1) {
			coor_t k;
			zHSP *next_hsp = &trellis->blocks->hsp[i+1];
			int    garbage_collect = 0;
			for (k = 0; k <= 2*BLOCK_OVERLAP; k++) {
				if (next_hsp->c_start + 2*BLOCK_OVERLAP >  k) {
					zQuickTracePartialTrellis(trellis, gmin, next_hsp->g_start + k - 1, cmin, next_hsp->c_start + 2*BLOCK_OVERLAP - k);
					garbage_collect = 1;
				}
			}
			if (garbage_collect == 1) {
				/* This number 32 comes from max explicit duration length = 30 and buffer 2 */
				zGarbageCollectPartialTrellis(trellis, MAX(next_hsp->g_start, 32+trellis->padding) - 32, MAX(next_hsp->c_start, 32+trellis->padding) - 32);
			}
		}
	}

	zFinishAlignmentForward(trellis, gmax, cmax);
	
	/* Viterbi trace back */
	zTrace2("traceback");
	afv = zMalloc(sizeof(zAFVec),"zRunPairViterbiAndForward afv");
	zInitAFVec(afv,10);

/* FOR U12 SCANNING * 
cmax = cmin;
 * FOR U12 SCANNING */
	zTracePartialPairTrellis(trellis, gmin, gmax, cmin, cmax, -1, afv);
	*path_score = trellis->cell[gmax][cmax][afv->last->state].score;

	return afv;
}

static void zTracePartialPairTrellis(zPairTrellis *trellis, coor_t g_begin, coor_t g_end, coor_t c_begin, coor_t c_end, int state, zAFVec *afv) {
	zHMM *hmm = trellis->hmm;
	int trace, max_state, max_state_bak;
	coor_t g_current, c_current;
	score_t max_score, cell_score;
	zAlnFeature af;

	
	/* Get the optimal end state */
	if (state < 0) {
		max_score = MIN_SCORE;
		for (max_state = 0; max_state < hmm->states; max_state++) {
			if (NULL == zGetCellArray(trellis, g_end, c_end)) continue;   /* Uninitialized cell */
			cell_score = zGetCurrentCell(trellis, g_end, c_end, max_state)->score;
			if (cell_score > max_score) {
				state = max_state;
				max_score = cell_score;
			}
		}
	}

	max_state_bak = state;
	g_current = g_end;
	c_current = c_end;

	while (g_current >= g_begin || c_current >= c_begin) {
		coor_t length = zGetCurrentCell(trellis, g_current, c_current, state)->length;
		coor_t genomic_start = g_current - length*zGetGenomicIncrement(hmm, state);
		coor_t cdna_start    = c_current - length*zGetCDnaIncrement(hmm, state);
		score_t score;
       
		trace = zGetCurrentCell(trellis, g_current, c_current, state)->trace;

		/* Calculate this state's score = emission_score + duration_score */

		score = zGetCurrentCell(trellis, g_current, c_current, state)->score;

		/* Add back the exit probabilities since it had been fixed by zFixInternalTransitions() */
		if (hmm->state[state].type == INTERNAL || hmm->state[state].type == GINTERNAL) {
			score += zScoreDurationGroup(hmm->dmap[hmm->state[state].duration], 1, trellis->genomic->gc);
		}

		/* MANI now removing initial probabilities from output for the LAST state and transition probabilities for the others */
		if (g_current == g_end && c_current == c_end) {
			score -= zGetInitProb(hmm, state, trellis->iiso_group);
		} 

		/* Init Prob should also be removed from the FIRST state but it is included in the score for (g_begin,c_begin). 
		   So it will be removed later when we remove previous cell's score. Different from the way zTree2AFL() does this, 
		   since tree->root doesn't have any score. */

		/* Remove transition score from output for non-FIRST states */
		if (!(genomic_start == g_begin && cdna_start == c_begin)) {
			score -= zGetTransitionScore(trellis->hmm, trace, state, trellis->tiso_group);
		}

		/* Remove previous state's score */
		score -= zGetCurrentCell(trellis, genomic_start, cdna_start, trace)->score;

		zClearAlnFeature(&af);

		af.name          = hmm->state[state].name;
		af.state         = state;
		af.strand        = hmm->state[state].strand;
		af.genomic       = trellis->genomic;
		af.genomic_start = genomic_start;
		af.genomic_end   = g_current;
		af.cdna          = trellis->cdna;
		af.cdna_start    = cdna_start;
		af.cdna_end      = c_current;
		af.score         = score;
		af.padding       = trellis->padding;

		/********************************************************************
		 * In donor and acceptor states, length is the length of the state, *
		 * but the actual sequence length is 2 times that (that comes from  *
		 * the increment[state][0] part). So we have to adjust the length   *
		 * for that. This is useful later in printing the alignment and     *
		 * converting pairagon output to est_genome output. I used to       *
		 * hardcode multiplying by 2 in the zWriteAlignment() function, but *
		 * there is no need for that now. It also maked the statement       *
		 *  length = f.genomic_end - f.genomic_start                        *
		 * valid in all cases.	                                            *
		 ********************************************************************/

		if (zGetGenomicIncrement(hmm, state) > 0) {
			length *= zGetGenomicIncrement(hmm, state);
		} else if (zGetCDnaIncrement(hmm, state) > 0) {
			length *= zGetCDnaIncrement(hmm, state);
		}
		af.length        = length;

		/********************************************************************
		 * start is incremented to get a fully closed interval. e.g.,
		 * 50 -- 100 is converted to 51 -- 100, since the previous state
		 * has ended in 50, and it belongs to the previous state. To 
		 * differentiate a single base state from a gap/unaligned state,
		 * I decremented the end by 1, so that start > end. This is checked
		 * in the zWriteAlignment() function to decide whether or not the 
		 * begin and end coordinates are printed out for this state.
		 *********************************************************************/

		if (zIsMatch(hmm->state[state].name)) { 
			coor_t iterator, matches = 0;
			for (iterator = 0; iterator < af.length; iterator++) {
				/* add 1 since it is not a fully closed interval */
				if (zGetDNAS5(af.genomic, af.genomic_start + 1 + iterator) == zGetDNAS5(af.cdna, af.cdna_start + 1 + iterator)) {
					matches++;
				}
			}
			af.percent_identity = (100*matches/af.length);
		} else {
			af.percent_identity = 0.0;
		}

		zPushAFVec(afv, &af);

		g_current = genomic_start;
		c_current = cdna_start;
		state = trace; 

		/* Traceback termination */
		if (g_current == g_begin && c_current == c_begin) {
			break;
		}

		if (-1 == trace) {
			zDie("Impossible to perform traceback from this state.");
		}

	}
	qsort(afv->elem, afv->size, sizeof(zAlnFeature), zAFPtrCmp);
}

void zTracePairTrellis (zPairTrellis *trellis, int state, zAFVec *afv) {
	zTrace2("#Tracing: %u-%u, %u-%u", trellis->padding - 1, trellis->genomic->length - 1 - trellis->padding, trellis->padding - 1, trellis->cdna->length - 1 - trellis->padding);
	zTracePartialPairTrellis(trellis, trellis->padding - 1, trellis->genomic->length - 1 - trellis->padding, trellis->padding - 1, trellis->cdna->length - 1 - trellis->padding, state, afv);
}

static void zQuickTracePartialTrellis(zPairTrellis *trellis, coor_t g_begin, coor_t g_end, coor_t c_begin, coor_t c_end) {
	int end_state;
	zHMM *hmm = trellis->hmm;
	coor_t g_current, c_current;

	zTrace3("#QuickTrace: %u--%u, %u--%u", g_begin, g_end, c_begin, c_end);

	if (NULL == zGetCellArray(trellis, g_end, c_end)) return;   /* Uninitialized cell array*/
	for (end_state = 0; end_state < hmm->states; end_state++) {
		int state;
		zPairTrellisCell *cell;

		/* Mark the end position */
		cell      = zGetCurrentCell(trellis, g_end, c_end, end_state);
		if (cell->score == MIN_SCORE) continue;
		zGetCurrentCell(trellis, g_end, c_end, 0)->keep = 1;

		/* Start tracing back */
		g_current = g_end;
		c_current = c_end;
		state     = end_state;
		while (g_current >= g_begin || c_current >= c_begin) {
			int length = cell->length;
			zTrace3("%s: %u, %u, %d from %s", zStrIdx2Char(trellis->hmm->state[state].name), g_current, c_current, length, zStrIdx2Char(hmm->state[cell->trace].name));
			g_current -= length*zGetGenomicIncrement(hmm, state);
			c_current -= length*zGetCDnaIncrement(hmm, state);
			state = cell->trace;
			cell      = zGetCurrentCell(trellis, g_current, c_current, state);
			zGetCurrentCell(trellis, g_current, c_current, 0)->keep = 1;

			/* Traceback termination */
			if (g_current == g_begin && c_current == c_begin) {
				break;
			}
			if (-1 == state) {
				zDie("Impossible to perform traceback from this state.");
			}
		}
	}
}

/*********************************************\
  Backward Probabilities
\*********************************************/

void zStartAlignmentBackward(zPairTrellis *trellis, coor_t gpos, coor_t cpos) {
	int     state;
	for (state = 0; state < trellis->hmm->states; state++) {
		if (trellis->backward != NULL) trellis->backward[gpos][cpos][state] = zGetInitProb(trellis->hmm, state, trellis->iiso_group);
	}
}

void zFinishAlignmentBackward(zPairTrellis *trellis, coor_t gpos, coor_t cpos) {
	int state;
	trellis->backward_score = MIN_SCORE;
	for (state = 0; state < trellis->hmm->states; state++) {
		score_t score  = zGetFixedInitProb(trellis->hmm, state, trellis->iiso_group, trellis->genomic->gc);
		trellis->backward_score = zFloatwiseScoreAdd(trellis->backward[gpos][cpos][state]+score, trellis->forward_score);
	}
}

static void zRunPartialPairBackward(zPairTrellis *trellis, coor_t gstart, coor_t gend, coor_t cstart, coor_t cend) {
	zHMM         *hmm = trellis->hmm;
	coor_t        genomic, cdna; /* iterators for sequence */
	int           state;         /* iterator for internal states */
	int           next;          /* iterator for next states */
	zIVec*        fjumps;

	/* induction */
	zTrace2("beginning induction");

	/* Make sure that all the required blocks are loaded */
	zCheckBackwardVariables(trellis, gstart, cstart);
	zTrace2("Calling (%u, %u) (%u, %u)", gstart, cstart, gend, cend);
	for (genomic = gend; genomic >= gstart; genomic--) {
		for (cdna = cend; cdna >= cstart; cdna--) {
			for (state = 0; state < hmm->states; state++) {
				fjumps = hmm->fmap[state];
				for (next = 0; next < fjumps->size; next++) {
					zGetPairBackTransFunc(hmm->state[state].type)
						(trellis, state, fjumps->elem[next], genomic, cdna);
				}
			}
		}
	}
}

void zRunPairBackward (zPairTrellis *trellis) {
	zDNA         *cdna = trellis->cdna;
	zDNA         *genomic = trellis->genomic;

	int           i;
	coor_t        gmin, cmin, gmax, cmax;

	gmin = trellis->padding - 1;
	cmin = trellis->padding - 1;
	gmax = genomic->length - 1 - trellis->padding;
	cmax = cdna->length - 1 - trellis->padding;
	 	 
	/* 	Viterbi and Forward Alg initialization */

	zAllocBackwardVars(trellis);
	zCheckBackwardVariables(trellis, gmax, cmax);

	/* induction - first pass */

	zTrace2("beginning first induction");
	zStartAlignmentBackward(trellis, gmax + 1, cmax + 1);

	zTrace2("running viterbi");
	for (i = trellis->blocks->hsps - 1; i >= 0; i--) {
		zHSP *hsp = &trellis->blocks->hsp[i];
		zRunPartialPairBackward(trellis, hsp->g_start, hsp->g_end, hsp->c_start, hsp->c_end);
	}
	zFinishAlignmentBackward(trellis, gmin, cmin);
fprintf(stdout, "backward: %f\n", (trellis->backward_score));
}

/*********************************************\
  Posterior Probability
\*********************************************/

void zComputePosteriorProbability(zPairTrellis* trellis, zAFVec* path) {
	int i, from;
	zAlnFeature* f = &path->elem[0];
	f->score += (trellis->backward[f->genomic_end][f->cdna_end][f->state]);
	f->score -= trellis->forward_score;
	f->score = 100*zScore2Float(f->score);
zWriteAlnFeature(stdout, f, 0, 0);
	for (i = 1; i < path->size - 1; i++) {
		f = &path->elem[i];
		from = path->elem[i-1].state;
		f->score += (trellis->forward[f->genomic_start][f->cdna_start][from] + zGetTransitionScore(trellis->hmm, from, f->state, trellis->tiso_group) + trellis->backward[f->genomic_end][f->cdna_end][f->state]);
		f->score -= trellis->forward_score;
		f->score = 100*zScore2Float(f->score);
zWriteAlnFeature(stdout, f, 0, 0);
	}
	f = &path->elem[i];
	from = path->elem[i-1].state;
	f->score += (trellis->forward[f->genomic_end][f->cdna_start][from] + zGetTransitionScore(trellis->hmm, from, f->state, trellis->tiso_group));
	f->score -= trellis->forward_score;
	f->score = 100*zScore2Float(f->score);
}

/*********************************************\
 Small Helper Functions 
\*********************************************/

void zShowPairTrellisCell(zPairTrellis *trellis, coor_t i, coor_t j, int k) {
	zPairTrellisCell *cell = zGetCurrentCell(trellis, i, j, k);
	fprintf(stdout, "trellis->cell[%u][%u][%s]={%f,%d,%d}\n", i, j, zStrIdx2Char(trellis->hmm->state[k].name), cell->score, cell->length, cell->trace);
}

void ShowPairTrellis(zPairTrellis *trellis){
	zHMM         *hmm = trellis->hmm;
	zDNA         *cdna = trellis->cdna;
	zDNA         *genomic = trellis->genomic;

	coor_t        i, j;        /* iterator for sequence */
	int           k;        /* iterator for internal states */

	printf("===================== Show PairTrellis ===================================\n");
	for(i = 0; i < genomic->length; i ++) { 
		for (j = 0; j < cdna->length; j++) {
			for (k = 0; k < hmm->states; k++) {
				zShowPairTrellisCell(trellis, i, j, k);
			}
		}
	}
	printf("======================= End PairTrellis =================================\n");

}

void zClearPairColumn(zPairTrellis* trellis, int gpos, int cpos) {
	int state;
	for (state=0; state < trellis->hmm->states; ++state) {
		trellis->cell[gpos][cpos][state].length = 0;
		trellis->cell[gpos][cpos][state].score  = MIN_SCORE;
		trellis->cell[gpos][cpos][state].trace = -1;
	}
}

#ifdef DEBUG

/* If not in DEBUG mode, these are #defines in zPairTrellis.h */

zPairTrellisCell* zGetCurrentCell(zPairTrellis *trellis, coor_t genomic, coor_t cdna, int state) {
	return &trellis->cell[genomic][cdna][state];
}

zPairTrellisCell* zGetNextCell(zPairTrellis *trellis, coor_t genomic, coor_t cdna, int next_state) {
	return &trellis->cell[genomic+zGetGenomicIncrement(trellis->hmm, next_state)][cdna+zGetCDnaIncrement(trellis->hmm, next_state)][next_state];
}

zPairTrellisCell* zGetPreviousCell(zPairTrellis *trellis, coor_t genomic, coor_t cdna, int current_state, int previous_state) {
	return &trellis->cell[genomic-zGetGenomicIncrement(trellis->hmm, current_state)][cdna-zGetCDnaIncrement(trellis->hmm, current_state)][previous_state];
}

zPairTrellisCell* zGetCellArray(zPairTrellis *trellis, coor_t genomic, coor_t cdna) {
	return trellis->cell[genomic][cdna];
}

#endif /* DEBUG */

#endif /* ZOE_PAIR_TRELLIS_C */
