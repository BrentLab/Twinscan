/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zEstseq.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_Estseq_C
#define ZOE_Estseq_C

#include "zEstseq.h"

void zInitEstseq (zEstseq *estseq, int bits){
	estseq->seq = NULL;
	estseq->bits = bits;
	estseq->digits = ceil(log10(estseq->bits));
	estseq->length = 0;
	estseq->def = NULL;
}

long zInitEstseqSequence(zSequence* seq, void* parent){
	coor_t length;
	int def_size;
	char c;
	char defline[4096];     /* should be big enough */
	int map_index;
	zEstseq* estseq = (zEstseq*)parent;

	seq->parent = parent;
	estseq->seq = seq;
	
	/* initial check for correct format */
	c = fgetc(seq->fp);
	if (c == EOF) return 0;
	if (c != '>') {
		zWarn("zInitEstseqSequence \">\" not found");
		return 0;
	}
	(void)ungetc(c, seq->fp);
	
    /* read the def line */
	(void)fgets(defline, sizeof(defline), seq->fp);
    def_size = strlen(defline);
	estseq->def = zMalloc(def_size +1,"zInitEstseqSequence def");
	(void)strcpy(estseq->def, defline);
	estseq->def[strlen(estseq->def) -1] = '\0'; /* remove newline */
	
	length = 0;
	map_index = 0;
	while ((c = fgetc(seq->fp)) != EOF) {
		if (c == '>') {
			(void)ungetc(c, seq->fp);
			break; /* next record found */
		}
		if (isspace((int)c)) continue; /* skip spaces */
		
		if(length == seq->file_map[map_index].seq_pos){
			ungetc(c,seq->fp);
			seq->file_map[map_index].file_pos = ftell(seq->fp);
			c = fgetc(seq->fp);
			map_index++;
		}
		length++;
	}
	
	seq->real_length = length;
	estseq->length = seq->real_length;

	return 1;
}

int zReadEstseqSequence(zSequence* seq, zSeqBlock* block){

	char c;
	coor_t index;
	/*zEstseq* estseq = (zEstseq*)seq->parent;*/

	if (ferror(seq->fp)) {
		zDie("zReadEstseqSequence file error 0");
	}

	index = 0;
	while ((c = fgetc(seq->fp)) != EOF) {
		if (c == '>') {
			(void)ungetc(c, seq->fp);
			break; /* next record found */
		}
		if (isspace((int)c)) continue; /* skip spaces */

		block->data[index] = c;
		
		index++;

		if (ferror(seq->fp)) {
			zDie("zReadEstseqSequence file error 2");
		}

		if(index == seq->block_size){
			break;
		}
	}

	while(index != seq->block_size){
		block->data[index] = '2';
		index++;
	}
	
	/* last check for errors */
	if (ferror(seq->fp)) {
		zDie("zReadEstseqSequence file error");
	}

	return 1;
}

void zLoadEstseqFromFasta(zEstseq* estseq, char* filename){
	zInitSequence(filename,estseq,zInitEstseqSequence,zReadEstseqSequence);
	estseq->seq->default_char = '0';
}

void zFreeEstseq(zEstseq *estseq) {
	if(estseq->seq != NULL){
		zFreeSequence(estseq->seq);
		zFree(estseq->seq);
		estseq->seq = NULL;
	}
	if(estseq->def != NULL){
		zFree(estseq->def);
		estseq->def = NULL;
	}
}

void zCopyEstseq(const zEstseq *estseq, zEstseq *copy) {
	copy->length = estseq->length;
	copy->bits = estseq->bits;
	copy->digits = estseq->digits;
	copy->def = zMalloc(strlen(estseq->def)+1,"zCopyEstseq def");
	(void)strcpy(copy->def, estseq->def);

	copy->seq = zMalloc(sizeof(zSequence), "zCopyEstseq: copy->seq");
	zCopySequence(estseq->seq,copy->seq);
	copy->seq->parent = copy;
}

void zReverseEstseq (zEstseq *estseq) {
	zReverseSequence(estseq->seq);
}

void zSetEstseqPadding(zEstseq* estseq,coor_t padding){
	zSetSequencePadding(estseq->seq,padding);
	estseq->length = estseq->seq->length;
}

char zGetEstseqSeq(zEstseq* estseq,coor_t pos){
	return zGetSequencePos(estseq->seq,pos);
}

char zGetEstseqS10(zEstseq* estseq,coor_t pos){
	char c = zGetSequencePos(estseq->seq,pos);
	return c - '0';
}


#endif
