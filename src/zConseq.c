/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zConseq.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_Conseq_C
#define ZOE_Conseq_C

#include "zConseq.h"

void zInitConseq (zConseq *conseq, int bits){
	conseq->seq = NULL;
	conseq->bits = bits;
	conseq->digits = ceil(log10(conseq->bits));
	conseq->length = 0;
	conseq->def = NULL;
}

long zInitConseqSequence(zSequence* seq, void *parent){
	coor_t length;
	int def_size;
	char c;
	char defline[4096];     /* should be big enough */
	int map_index;
	zConseq* conseq = (zConseq*)parent;

	seq->parent = parent;
	conseq->seq = seq;
	
	/* initial check for correct format */
	c = fgetc(seq->fp);
	if (c == EOF) return 0;
	if (c != '>') {
		zWarn("zInitConseqSequence \">\" not found");
		return 0;
	}
	(void)ungetc(c, seq->fp);
	
    /* read the def line */
	(void)fgets(defline, sizeof(defline), seq->fp);
    def_size = strlen(defline);
	conseq->def = zMalloc(def_size +1,"zInitConseqSequence def");
	(void)strcpy(conseq->def, defline);
	conseq->def[strlen(conseq->def) -1] = '\0'; /* remove newline */
	
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
	conseq->length = seq->real_length;

	return 1;
}

int zReadConseqSequence(zSequence* seq, zSeqBlock* block){

	char c;
	coor_t index;
	/*zConseq* conseq = (zConseq*)seq->parent;*/

	if (ferror(seq->fp)) {
		zDie("zReadConseqSequence file error 0");
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
			zDie("zReadConseqSequence file error 2");
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
		zDie("zReadConseqSequence file error");
	}

	return 1;
}

void zLoadConseqFromFasta(zConseq* conseq, char* filename){
	zInitSequence(filename,conseq,zInitConseqSequence,zReadConseqSequence);
	conseq->seq->default_char = '2';
}

void zFreeConseq(zConseq *conseq) {
	if(conseq->seq != NULL){
		zFreeSequence(conseq->seq);
		zFree(conseq->seq);
		conseq->seq = NULL;
	}
	if(conseq->def != NULL){
		zFree(conseq->def);
		conseq->def = NULL;
	}
}

void zCopyConseq(const zConseq *conseq, zConseq *copy) {
	copy->length = conseq->length;
	copy->bits = conseq->bits;
	copy->digits = conseq->digits;
	copy->def = zMalloc(strlen(conseq->def)+1,"zCopyConseq def");
	(void)strcpy(copy->def, conseq->def);

	copy->seq = zMalloc(sizeof(zSequence), "zCopyConseq: copy->seq");
	zCopySequence(conseq->seq,copy->seq);
	copy->seq->parent = copy;
}

void zReverseConseq (zConseq *conseq) {
	zReverseSequence(conseq->seq);
}

void zSetConseqPadding(zConseq* conseq,coor_t padding){
	zSetSequencePadding(conseq->seq,padding);
	conseq->length = conseq->seq->length;
}

char zGetConseqSeq(zConseq* conseq,coor_t pos){
	return zGetSequencePos(conseq->seq,pos);
}

char zGetConseqS10(zConseq* conseq,coor_t pos){
	char c = zGetSequencePos(conseq->seq,pos);
	return c - '0';
}


#endif
