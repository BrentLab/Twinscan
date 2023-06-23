/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
	zStopSeq.c - part of the ZOE library for genomic analysis
 
	Copyright (C) 2001-2002

\******************************************************************************/

#ifndef ZOE_StopSeq_C
#define ZOE_StopSeq_C

#include "zStopSeq.h"

static const coor_t BLOCK_SIZE = 100000;
static const int BLOCK_COUNT = 4;

static short id = 1;

struct zStopBlock{
	coor_t*  prev;
	coor_t*  next;
	short*   var;
	coor_t   pos;
}; 
typedef struct zStopBlock zStopBlock;

void* zCreateStopBlock(){
	zStopBlock* block = zMalloc(sizeof(zStopBlock),"zCreateStopBlock block");
	block->pos = -1;
	block->prev = zMalloc(BLOCK_SIZE*sizeof(coor_t),"zCreateStopBlock prev");
	block->next = zMalloc(BLOCK_SIZE*sizeof(coor_t),"zCreateStopBlock next");
	block->var = zMalloc(BLOCK_SIZE*sizeof(short),"zCreateStopBlock var");
	return (void*)block;
}

void zFreeStopBlock(void* v){
	zStopBlock* block = (zStopBlock*)v;
	block->pos = -1;
	zFree(block->prev);
	zFree(block->next);
	zFree(block->var);
	zFree(block);
}

void zResetStopBlock(void* v){
	zStopBlock* block = (zStopBlock*)v;
	block->pos = -1;
}

static bool zCheckDNAForCurrentStopCodon(zDNA* dna, coor_t pos){
	if((zGetDNAS5(dna,pos) == 3) && /* pos 0 is T */
	   (((zGetDNAS5(dna,pos+1) == 0) &&/* pos 1 is (T)A */
		 ((zGetDNAS5(dna,pos+2) == 0) || /* pos 2 is (TA)A */
		  (zGetDNAS5(dna,pos+2) == 2))) || /* pos 2 is (TA)G */
		((zGetDNAS5(dna,pos+1) == 2) && /* pos 1 is (T)G */
		 (zGetDNAS5(dna,pos+2) == 0)))){ /* pos 2 is (TG)A */
		return true;
	}
	return false;
}

static bool zCheckDNAForPossibleStopCodon(zDNA* dna, coor_t pos){
	bool t = zDNACheckT(dna,pos);
	bool a1 = zDNACheckA(dna,pos+1);
	bool a2 = zDNACheckA(dna,pos+2);
	bool g1 = zDNACheckG(dna,pos+1);
	bool g2 = zDNACheckG(dna,pos+2);
	if(t && ((a1 && (a2 || g2)) || (g1 && a2))){
		return true;
	}
	return false;
}

static coor_t zCalculatePrevStopSeqValue(zStopSeq* ss, coor_t pos){
    while(pos > 2 && 
		  zCheckDNAForPossibleStopCodon(ss->dna,pos) == false){
		pos -= 3;
    }
    return pos;	
}

static coor_t zCalculateNextStopSeqValue(zStopSeq* ss, coor_t pos){
    while(pos < ss->dna->length - 3 && 
		  zCheckDNAForPossibleStopCodon(ss->dna,pos) == false){
		pos += 3;
    }
    return pos;	
}

static zStopBlock* zLoadStopSeq(zStopSeq* ss, coor_t pos){
	coor_t i,j,max_pos;
	zStopBlock* block;
	if(zListGetSize(&ss->value) < BLOCK_COUNT){
		block = zListAddLast(&ss->value);
	}
	else{
		block = zListMoveLast(&ss->value);
		zListMoveCurrentToFirst(&ss->value);
	}
	max_pos = BLOCK_SIZE-1;
	if(ss->dna->length < pos+max_pos){
		max_pos = ss->dna->length - pos;
	}
	block->pos = pos;

	/* clear var array */
	for(i = 0; i < BLOCK_SIZE; i++){
		block->var[i] = -1;
	}
	for(i = 0; (int)i < ss->dna->seq->var_count; i++){
		j = ss->dna->seq->variants[i]->pos;
		if(j < block->pos) continue;
		if(j >= block->pos + BLOCK_SIZE) break;
		block->var[j-block->pos] = i;		
	}
	
	/* fill in prev array */
	block->prev[0] = zCalculatePrevStopSeqValue(ss,pos);
	block->prev[1] = zCalculatePrevStopSeqValue(ss,pos+1);
	block->prev[2] = zCalculatePrevStopSeqValue(ss,pos+2);
	j = 2;
	for(i = pos+3; i <= pos+max_pos; i++){
		j++;
		if(zCheckDNAForPossibleStopCodon(ss->dna,i) == true){
			block->prev[j] = i;
		}
		else{
			block->prev[j] = block->prev[j-3];
		}
	}
	/* fill in next array */
	if(max_pos < BLOCK_SIZE-1){
		block->next[max_pos] = ss->dna->length;
		block->next[max_pos-1] = ss->dna->length - 1;
		block->next[max_pos-2] = ss->dna->length - 2;
	}
	else{
		block->next[max_pos] = zCalculateNextStopSeqValue(ss,pos+max_pos);
		block->next[max_pos-1] = zCalculateNextStopSeqValue(ss,pos+max_pos-1);
		block->next[max_pos-2] = zCalculateNextStopSeqValue(ss,pos+max_pos-2);
	}
	j = max_pos-2;
	for(i = pos+max_pos-3; i != pos; i--){
		j--;
		if(block->prev[j] == i){
			block->next[j] = i;
		}
		else{
			block->next[j] = block->next[j+3];
		}
	}
	/* pull out base case to avoid infinite loop when pos=0 since i is unsigned */
	if(block->prev[0] == pos){
		block->next[0] = pos;
	}
	else{
		block->next[0] = block->next[3];
	}
	return block;
}

static zStopBlock* zGetStopBlock(zStopSeq* ss, coor_t pos){
	zStopBlock* block = zListMoveFirst(&ss->value);
	while(block != NULL){
		if((block->pos <= pos) && (block->pos + BLOCK_SIZE > pos)){
			zListMoveCurrentToFirst(&ss->value);
			return block;
		}
		block = zListMoveNext(&ss->value);
	}
	pos = (int)(pos/BLOCK_SIZE) * BLOCK_SIZE;
	return zLoadStopSeq(ss,pos);
}

void zInitStopSeq (zStopSeq* ss,zDNA* dna){
	ss->id = id++;
	ss->dna = dna;
	zInitList(&ss->value,zCreateStopBlock,zFreeStopBlock,zResetStopBlock);
	zLoadStopSeq(ss,0);	
}

void zFreeStopSeq (zStopSeq* ss){
	zFreeList(&ss->value);
}

coor_t zGetStopSeqPrevPos(zStopSeq* ss, coor_t pos){
	zStopBlock* block;
	coor_t sc_pos;
	block = zGetStopBlock(ss,pos);
	sc_pos = block->prev[pos - block->pos];
	if(sc_pos < block->pos){
		return zGetStopSeqPrevPos(ss, sc_pos);
	}
	else{
		if((block->var[sc_pos - block->pos] != -1) ||
		   (block->var[sc_pos - block->pos+1] != -1) ||
		   (block->var[sc_pos - block->pos+2] != -1)){
			if(zCheckDNAForCurrentStopCodon(ss->dna,sc_pos)){
				return sc_pos;
			}
			else{
				if(sc_pos < 3){
					return sc_pos;
				}
				return zGetStopSeqPrevPos(ss, sc_pos-3);
			}
		}
		else{
			return sc_pos;
		}
	}
}

coor_t zGetStopSeqNextPos(zStopSeq* ss, coor_t pos){
	zStopBlock* block;
	coor_t sc_pos;
	block = zGetStopBlock(ss,pos);
	sc_pos = block->next[pos - block->pos];
	if(sc_pos >= block->pos+BLOCK_SIZE){
		return zGetStopSeqNextPos(ss, sc_pos);
	}
	else{
		if((block->var[sc_pos - block->pos] != -1) ||
		   (block->var[sc_pos - block->pos+1] != -1) ||
		   (block->var[sc_pos - block->pos+2] != -1)){
			if(zCheckDNAForCurrentStopCodon(ss->dna,sc_pos)){
				return sc_pos;
			}
			else{
				if(sc_pos > ss->dna->length-3){
					return sc_pos;
				}
				return zGetStopSeqNextPos(ss, sc_pos+3);
			}
		}
		else{
			return sc_pos;
		}
	}
}

coor_t zGetStopSeqMinPrevPos(zStopSeq* ss, coor_t pos){
	zStopBlock* block;
	coor_t sc_pos;
	block = zGetStopBlock(ss,pos);
	sc_pos = block->prev[pos - block->pos];
	if(sc_pos < block->pos){
		return zGetStopSeqMinPrevPos(ss, sc_pos);
	}
	else{
		if((block->var[sc_pos - block->pos] != -1) ||
		   (block->var[sc_pos - block->pos+1] != -1) ||
		   (block->var[sc_pos - block->pos+2] != -1)){
			if(sc_pos < 3){
				return sc_pos;
			}
			return zGetStopSeqMinPrevPos(ss, sc_pos-3);
		}
		else{
			return sc_pos;
		}
	}
}

coor_t zGetStopSeqMaxNextPos(zStopSeq* ss, coor_t pos){
	zStopBlock* block;
	coor_t sc_pos;
	block = zGetStopBlock(ss,pos);
	sc_pos = block->next[pos - block->pos];
	if(sc_pos >= block->pos+BLOCK_SIZE){
		return zGetStopSeqMaxNextPos(ss, sc_pos);
	}
	else{
		if((block->var[sc_pos - block->pos] != -1) ||
		   (block->var[sc_pos - block->pos+1] != -1) ||
		   (block->var[sc_pos - block->pos+2] != -1)){
			if(sc_pos > ss->dna->length-3){
				return sc_pos;
			}
			return zGetStopSeqMaxNextPos(ss, sc_pos+3);
		}
		else{
			return sc_pos;
		}
	}
}

#endif
