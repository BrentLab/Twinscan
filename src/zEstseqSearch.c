/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
	zEstseqSearch.c - part of the ZOE library for genomic analysis
 
	Copyright (C) 2001-2002

\******************************************************************************/

#ifndef ZOE_EstseqSearch_C
#define ZOE_EstseqSearch_C

#include "zEstseqSearch.h"


static const coor_t BLOCK_SIZE = 100000;
static const int BLOCK_COUNT = 5;

static short MAX_EST_SYMBOLS = 3;

static short id = 1;

static short char_map[256];

struct zEstseqSearchBlock{
	coor_t**  prev;
	coor_t**  next;
	coor_t   pos;
}; 
typedef struct zEstseqSearchBlock zEstseqSearchBlock;

static void set_char_map(){
	int i;
	for(i = 0; i < 256; i++){
		char_map[i] = -1;
	}
	for(i = 0; i < MAX_EST_SYMBOLS; i++){
		char_map[(int)'0'+i] = i;
	}
}

void* zCreateEstseqSearchBlock(){
	coor_t i;
	zEstseqSearchBlock* block = 
		zMalloc(sizeof(zEstseqSearchBlock),"zCreateEstseqSearchBlock block");
	block->pos = -1;
	block->prev = zMalloc(BLOCK_SIZE*sizeof(coor_t*),"zCreateEstseqSearchBlock prev");
	block->next = zMalloc(BLOCK_SIZE*sizeof(coor_t*),"zCreateEstseqSearchBlock next");
	for(i = 0; i < BLOCK_SIZE; i++){
		block->prev[i] = zMalloc(MAX_EST_SYMBOLS*sizeof(coor_t),"zCreateEstseqSearchBlock prev[i]");
		block->next[i] = zMalloc(MAX_EST_SYMBOLS*sizeof(coor_t),"zCreateEstseqSearchBlock next[i]");
	}
	return (void*)block;
}

void zFreeEstseqSearchBlock(void* v){
	coor_t i;
	zEstseqSearchBlock* block = (zEstseqSearchBlock*)v;
	block->pos = -1;
	for(i = 0; i < BLOCK_SIZE; i++){
		zFree(block->prev[i]);
		zFree(block->next[i]);
	}
	zFree(block->prev);
	zFree(block->next);
	zFree(block);
}

void zResetEstseqSearchBlock(void* v){
	zEstseqSearchBlock* block = (zEstseqSearchBlock*)v;
	block->pos = -1;
}

static zEstseqSearchBlock* zLoadEstseqSearch(zEstseqSearch* es, coor_t pos){
	coor_t i,max_pos;
	int j;
	
	zEstseqSearchBlock* block;
	if(zListGetSize(&es->value) < BLOCK_COUNT){
		block = zListAddLast(&es->value);
	}
	else{
		block = zListMoveLast(&es->value);
		zListMoveCurrentToFirst(&es->value);
	}
	max_pos = BLOCK_SIZE-1;
	if(es->estseq->length < pos+max_pos){
		max_pos = es->estseq->length - pos;
	}
	block->pos = pos;

	if(block->pos == 0){
		/*first position has value 0 */
		for(j = 0; j < es->symbols; j++){
			es->temp[j] = 0;
		}
	}
	else{
		for(j = 0; j < es->symbols; j++){
			es->temp[j] = -1;
		}
	}
	/* fill in prev array */
	for(i = 0; i < BLOCK_SIZE; i++){
		es->temp[char_map[(int)zGetEstseqSeq(es->estseq,i+block->pos)]] = i+block->pos;
		for(j = 0; j < es->symbols; j++){
			block->prev[i][j] = es->temp[j];
			/*printf("PREV %d %d = %d\n",i,j,es->temp[j]);*/
		}
	}


	if(block->pos+BLOCK_SIZE > es->estseq->length){
		/*set values beyond end of estseq */
		for(j = 0; j < es->symbols; j++){
			es->temp[j] = es->estseq->length;
		}
		max_pos = es->estseq->length - block->pos - 1;
		for(i = max_pos+1; i < BLOCK_SIZE; i++){
			for(j = 0; j < es->symbols; j++){
				block->next[i][j] = es->estseq->length;
			}
		}
	}
	else{
		max_pos = BLOCK_SIZE-1;
		for(j = 0; j < es->symbols; j++){
			es->temp[j] = -1;
		}
	}
	/*fill in next array*/
	for(i = max_pos;i > 0; i--){
		es->temp[char_map[(int)zGetEstseqSeq(es->estseq,i+block->pos)]] = i+block->pos;
		for(j = 0; j < es->symbols; j++){
			block->next[i][j] = es->temp[j];
		}
	}
	es->temp[char_map[(int)zGetEstseqSeq(es->estseq,block->pos)]] = block->pos;
	for(j = 0; j < es->symbols; j++){
		block->next[0][j] = es->temp[j];
	}

	return block;
}

static zEstseqSearchBlock* zGetEstseqSearchBlock(zEstseqSearch* es, coor_t pos){
	zEstseqSearchBlock* block = zListMoveFirst(&es->value);
	while(block != NULL){
		if((block->pos <= pos) && (block->pos + BLOCK_SIZE > pos)){
			zListMoveCurrentToFirst(&es->value);
			return block;
		}
		block = zListMoveNext(&es->value);
	}
	pos = (int)(pos/BLOCK_SIZE) * BLOCK_SIZE;
	return zLoadEstseqSearch(es,pos);
}

void zInitEstseqSearch (zEstseqSearch* es,zEstseq* estseq){
	es->id = id++;
	es->estseq = estseq;
	es->symbols = 3; /*EVAN this shouldn't be hard coded*/
	MAX_EST_SYMBOLS = es->symbols;
	set_char_map();
	es->temp = zMalloc(sizeof(coor_t)*es->symbols,"zLoadEstseqSearch temp");
	zInitList(&es->value,zCreateEstseqSearchBlock,zFreeEstseqSearchBlock,zResetEstseqSearchBlock);
	zLoadEstseqSearch(es,0);	
}

void zFreeEstseqSearch (zEstseqSearch* es){
	zFree(es->temp);
	zFreeList(&es->value);
}

coor_t zGetEstseqPrevPos(zEstseqSearch* es, char estchar, coor_t pos){
	return zGetEstseqPrevPosBounded(es,estchar,pos,0);
}

coor_t zGetEstseqNextPos(zEstseqSearch* es, char estchar, coor_t pos){
	return zGetEstseqNextPosBounded(es,estchar,pos,es->estseq->length-1);
}

coor_t zGetEstseqPrevPosBounded(zEstseqSearch* es, char estchar, coor_t pos, coor_t min){
	zEstseqSearchBlock* block;
	coor_t old_pos;
	coor_t char_pos;
	block = zGetEstseqSearchBlock(es,pos);
	char_pos = block->prev[pos - block->pos][char_map[(int)estchar]];
	if(char_pos == (coor_t)-1){
		if(min >= block->pos){
			return min;
		}
		old_pos = block->pos;
		/* need top check previous block */
		char_pos = zGetEstseqPrevPosBounded(es, estchar, block->pos-1,min);
		/* set value  so we don't need to do this again */
		if(old_pos == block->pos && char_pos != min){
			/* in case block gets reloaded by recursive zGetEstseqPrevPos */
			block->prev[pos-block->pos][char_map[(int)estchar]] = char_pos;
		}
		return char_pos;
	}
	else{
		return zCoorMax(char_pos,min);
	}
}

coor_t zGetEstseqNextPosBounded(zEstseqSearch* es, char estchar, coor_t pos, coor_t max){
	zEstseqSearchBlock* block;
	coor_t char_pos;
	coor_t old_pos;
	block = zGetEstseqSearchBlock(es,pos);
	char_pos = block->next[pos - block->pos][char_map[(int)estchar]];
	if(char_pos == (coor_t)-1){
		if(block->pos+BLOCK_SIZE > max){
			return max;
		}
		old_pos = block->pos;
		/* need top check previous block */
		char_pos = zGetEstseqNextPosBounded(es, estchar, block->pos+BLOCK_SIZE,max);
		/* set value  so we don't need to do this again */
		if(old_pos == block->pos && char_pos != max){
			/* in case block gets reloaded by recursive zGetEstseqNextPos */
			block->next[pos-block->pos][char_map[(int)estchar]] = char_pos;
		}
		return char_pos;
	}
	else{
		return zCoorMin(char_pos,max);
	}
}

#endif
