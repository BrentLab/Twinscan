/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zSfeature.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_SFEATURE_C
#define ZOE_SFEATURE_C

#include "zSfeature.h"

void zClearSfeature (zSfeature *f) {
	f->name   = -1; 
	f->state  = -1; 
	f->start  = UNDEFINED_COOR; 
	f->end    = UNDEFINED_COOR;
	f->strand = UNDEFINED_STRAND; 
	f->score  = MIN_SCORE; 
	f->lfrag  = UNDEFINED_FRAME;
	f->rfrag  = UNDEFINED_FRAME; 
	f->frame  = UNDEFINED_FRAME; 
	f->group  = NULL;
	f->cdscons_score = MIN_SCORE;
	f->cdsest_score = MIN_SCORE;
	f->align_score = MIN_SCORE;
#ifdef DEBUG
	f->begin_score   = MIN_SCORE; 
	f->end_score     = MIN_SCORE; 
	f->content_score = MIN_SCORE; 
#endif /*DEBUG */
}

void zFreeSfeature (zSfeature *f) {
	zFree(f->group);
	zClearSfeature(f);
}

void zWriteSfeature (FILE *stream, const zSfeature *f) {
	char start[16], end[16], strand[8], score[8],
		left[8], right[8], frame[8];

	zCoor2Text(f->start +1, start);  /* note: must add 1 to start and end because */
	zCoor2Text(f->end +1, end);      /* most people use 1-based coordinates */
	zScore2Text(f->score, score);
	zFrame2Text(f->lfrag, left);
	zFrame2Text(f->rfrag, right);
	zFrame2Text(f->frame, frame);
	zStrand2Text(f->strand, strand);

	(void)fprintf(stream, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
				  zStrIdx2Char(f->name), start, end, strand, score, left, right, frame);

	if (f->group == NULL) (void)fprintf(stream, "\n");
	else                  (void)fprintf(stream, "\t%s\n", f->group);
}

int zReadSfeature (FILE *stream, zSfeature *f) {
	char line[2048], name[17], start[17], end[17], strand[9], score[9],
		left[9], right[9], frame[9], group[65];
	int group_defined, short_form;
	coor_t tmp;

	zClearSfeature(f);
	
	do {
		if (fgets(line, sizeof(line)-1, stream) == NULL) return 0;
	} while ('#' == line[0]);
	
	if (sscanf(line, "%16s %16s %16s %8s %8s %8s %8s %8s %64s", /* group */
			   name, start, end, strand, score, left, right, frame, group) == 9) {
		group_defined = 1;
		short_form    = 0;
	} else if (sscanf(line, "%16s %16s %16s %8s %8s %8s %8s %8s", /* no group */
					  name, start, end, strand, score, left, right, frame) == 8) {
		group_defined = 0;
		short_form    = 0;
	} else if (sscanf(line, "%16s %16s %16s %64s", /* short form w/ group */
					  name, start, end, group) == 4) {
		group_defined = 1;
		short_form    = 1;
	} else if (sscanf(line, "%16s %16s %16s", /* short form w/o group */
					  name, start, end) == 3) {
		group_defined = 0;
		short_form    = 1;
	} else {
		zWarn("zReadSfeature format not valid");
		return 0;
	}
	/* consider adding gff someday */
	
	/* set name to enum */
	f->name = zChar2StrIdx(name);
	
	if (short_form) {
		f->start  = zText2Coor(start);
		f->end    = zText2Coor(end);
		if (f->start < f->end) {
			f->strand = '+';
		} else {
			f->strand = '-';
			tmp = f->start;
			f->start = f->end;
			f->end = tmp;
		}
		f->score  = MIN_SCORE;
		f->lfrag  = 0;
		f->rfrag  = 0;
		f->frame  = 0;
	} else {
	
		f->start  = zText2Coor(start);
		f->end    = zText2Coor(end);
		f->strand = zText2Strand(strand);
		f->score  = zText2Score(score);
		f->lfrag  = zText2Frame(left);
		f->rfrag  = zText2Frame(right);
		f->frame  = zText2Frame(frame);
	}
	
	/* subtract 1 from start and end, zoe uses 0-based coordinates internally */
	f->start--;
	f->end--;
	
	/* set group */
	if (group_defined) {
		f->group = zMalloc(strlen(group) +1, "zReadSfeature group");
		(void)strcpy(f->group, group);
	} else {
		f->group = NULL;
	}
	
	if (zVerifySfeature(f)) {
		return 1;
	} else {
		zWarn("zReadSfeature parse error (%s)", line);
		zFreeSfeature(f);
		return 0;
	}
}

int zVerifySfeature (const zSfeature *f) {
	if (f->start > f->end)          return 0;
	if (f->start == UNDEFINED_COOR) return 0;
	if (f->end   == UNDEFINED_COOR) return 0;
	if (f->strand != '-' && f->strand != '+' && f->strand != '=')      return 0;
	if (f->lfrag != UNDEFINED_FRAME && (f->lfrag < 0 || f->lfrag > 2)) return 0;
	if (f->lfrag != UNDEFINED_FRAME && (f->rfrag < 0 || f->rfrag > 2)) return 0;
	/*if (f->lfrag != UNDEFINED_FRAME && (f->frame < 0 || f->frame > 2)) return 0;*/
	
	/* passed the tests */
	return 1;
}

void zCopySfeature (const zSfeature *orig, zSfeature *copy) {
	copy->name   = orig->name;
	copy->start  = orig->start;
	copy->end    = orig->end;
	copy->strand = orig->strand;
	copy->score  = orig->score;
	copy->lfrag  = orig->lfrag;
	copy->rfrag  = orig->rfrag;
	copy->frame  = orig->frame;
	copy->state  = orig->state;

#ifdef DEBUG
	copy->begin_score = orig->begin_score;
	copy->end_score = orig->end_score;
	copy->content_score = orig->content_score;
#endif /* DEBUG */
	copy->cdscons_score = orig->cdscons_score;
	copy->cdsest_score = orig->cdsest_score;
	copy->align_score = orig->align_score;
	/* 	copy->begincons_score = orig->begincons_score; */
	/* 	copy->endcons_score = orig->endcons_score; */


	copy->group  = NULL;
	/*if (orig->group != NULL) {
	  copy->group = zMalloc(strlen(orig->group) +1, "zCopySfeature group");	
	  (void)strcpy(copy->group, orig->group);
	  }*/
}

void zAntiSfeature(zSfeature *f, coor_t length) {
	coor_t tmp_start;

	tmp_start = f->start;
	f->start = length - f->end - 1;
	f->end = length - tmp_start -1;
}

int zSfeatureCmp (const zSfeature *f1, const zSfeature *f2) {

	if      (f1->start < f2->start) return -1;
	else if (f1->start > f2->start) return 1;
	else {
		if      (f1->end < f2->end) return -1;
		else if (f1->end > f2->end) return 1;
		else {
			if      (f1->strand < f2->strand) return -1;
			else if (f1->strand > f2->strand) return 1;
			else return f1->name - f2->name;
		}
	}
}

int zSfPtrCmp (const void *v1, const void *v2) {
	return zSfeatureCmp( *(zSfeature**)v1, *(zSfeature**)v2 );
}

/* zSFVec stuff */

void zInitSFVec (zSFVec *vec, int limit) {
	vec->size = 0;
	vec->limit = limit;
	if (vec->limit > 0) vec->elem = zMalloc(limit * sizeof(zSfeature), "zInitSFVec");
	else                vec->elem = NULL;
}

void zPushSFVec (zSFVec *vec, const zSfeature *feature) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zRealloc(vec->elem, vec->limit * sizeof(zSfeature), "zPushSFVec");
	}
	zCopySfeature(feature, &vec->elem[vec->size]);
	vec->last = &vec->elem[vec->size];
	vec->size++;
}

void zResetSFVec (zSFVec *vec) {
	vec->size = 0;
}

void zFreeSFVec (zSFVec *vec) {
	zFree(vec->elem);
}

void zFlipSFVec (zSFVec* sfv) {
	int i,other;
	zSfeature tmp;

	for (i = 0; i < sfv->size/2; ++i) {
		other = sfv->size - 1 - i;
		zCopySfeature(&sfv->elem[i],     &tmp);
		zCopySfeature(&sfv->elem[other], &sfv->elem[i]);
		zCopySfeature(&tmp,              &sfv->elem[other]);
	}
}

void zWriteSFVec (FILE* stream, const zSFVec* sfv) {
	int i;
	for (i = 0; i < sfv->size; ++i) {
		zWriteSfeature(stream, &sfv->elem[i]);
	}
}

zSFVec* zReadSFVec (FILE* stream) {
	zSFVec    *sfv;
	zSfeature  f;

	sfv = zMalloc(sizeof(zSFVec), "zReadSFVec");
	zInitSFVec(sfv, 1);
	while(zReadSfeature(stream, &f)) {
		zPushSFVec(sfv, &f);
	}
	return sfv;
}

void zTranslateSFVec (zSFVec* sfv, int dist) {
	int i;
	for (i=0; i < sfv->size; ++i) {
		sfv->elem[i].start += dist;
		sfv->elem[i].end   += dist;
	}
}

void zTranslateSFList (zSFList* sfl, int dist) {
	zSfeature* f = zSFListMoveFirst(sfl);
	while(f != NULL){
		f->start += dist;
		f->end   += dist;
		f = zSFListMoveNext(sfl);
	}
}

/*
 * Returns a phase type that directly corresponds to the right fragment
 */
zPhase_t zSFEndPhase (zDNA *dna, zSfeature *exon) {


	/* 	return exon->rfrag; */


	int terminal_base;
	coor_t end;
	/*  	end = ('-' == exon->strand) ? dna->length - 1 - exon->start : exon->end; */
 	end = ('-' == exon->strand) ? dna->length - 1 - exon->start : exon->end;
		
	switch (exon->rfrag) {
	case 1:  terminal_base = zGetDNAS5(dna,end) +1; break;
	case 2:  terminal_base = (zGetDNAS5(dna,end -1) +1) * 5
				 +  zGetDNAS5(dna,end) +1; break;
	default: terminal_base = 0;
	}

	switch (terminal_base) {
	case 0:
		return Phase0;
	case 1: case 2: case 3: case 5: /* terminal base, not T */
		return Phase1;
	case 4:  /* T  */
		return Phase1T;
	case 21: /* TA */
		return Phase2TA;
	case 23: /* TG */
		return Phase2TG;
	default:
		return Phase2;
	}
}

/*
 * Assigns the following phases to the beginning of a feature. These assignments
 * are somewhat arbitrary, and should be investigated further. (!!!)
 *
 *            lfrag     bases
 *            -----     -----
 * Phase0:      0
 * Phase1:      2        NOT AA, AG, or GA
 * Phase1T:     2        AA, AG, or GA
 * Phase2:      1        C or T
 * Phase2TA:    1        A
 * Phase2TG:    1        G
 */
zPhase_t zSFBeginPhase (zDNA *dna, zSfeature *exon) {

	
	
	/* 	return exon->lfrag; */

	int initial_base;
	coor_t start;



 	start = ('-' == exon->strand) ? dna->length - 1 - exon->end : exon->start;
	/*  	start = ('-' == exon->strand) ? dna->length - 1 - exon->end : exon->start; */
	switch (exon->lfrag) {
	case 1:  initial_base =  zGetDNAS5(dna,start)   + 1; break;
	case 2:  initial_base = (zGetDNAS5(dna,start)   + 1) *5
				 +  zGetDNAS5(dna,start+1) + 1; break;
	default: initial_base = 0;
	}
	
	switch (initial_base) {
	case 0:
		return Phase0;
	case 1: /* single A */
		return Phase2TA;
		/* 	case 2: case 4: /\* single C or single T *\/ */
	case 2: case 4: case 5:/* single C or single T or single N */
		return Phase2;
	case 3: /* single G*/
		return Phase2TG;
	case 6: case 8: case 16: /* AA, AG, or GA */
		return Phase1T;
	default:
		return Phase1;
	}
}


int zCompatPhases(int positive_strand, zPhase_t from, zPhase_t to) {


	/* 	return ((from+to) % 3 == 0); */

	zPhase_t temp;
	if (0 == positive_strand) {
		temp = from;
		from = to;
		to = temp;
	}


	/*   switch (from) { */
	/*   case Phase0: */
	/* 	return (Phase0 == to); */
	/*   case Phase1: */
	/* 	return (Phase1 == to) | (Phase1T == to); */
	/*   case Phase1T: */
	/*  	  return (Phase1 == to) | (Phase1 == to); */
	/*   case Phase2: */
	/* 	return (Phase2 == to) | (Phase2TA == to) | (Phase2TG == to); */
	/*   case Phase2TA: */
	/* 	return (Phase2 == to) | (Phase2TA == to) | (Phase2TG == to); */
	/*   case Phase2TG: default: */
	/* 	  return (Phase2 == to) | (Phase2TA == to) | (Phase2TG == to); */
	/*   } */

	switch (from) {
	case Phase0:
		return (Phase0 == to);
	case Phase1:
		return (Phase1 == to) | (Phase1T == to);
	case Phase1T:
		return (Phase1 == to);
	case Phase2:
		return (Phase2 == to) | (Phase2TA == to) | (Phase2TG == to);
	case Phase2TA:
		return (Phase2 == to);
	case Phase2TG: default:
		return (Phase2 == to) | (Phase2TG == to);
	}

}


/**********************************************************************\
  zSFList - linked list of zSfeature objects
\**********************************************************************/

static unsigned int zSFListCount = 0; 
static zSFListNode* zSFListDeadNodes;

void zIncrementSFListPool(){
	if(zSFListCount == 0){
		zSFListDeadNodes = NULL;	
	}
	zSFListCount++;
}

void zDecrementSFListPool(){
	zSFListNode* n;
    zSFListCount--;
	if(zSFListCount == 0){
		while(zSFListDeadNodes != NULL){
			n = zSFListDeadNodes;
			zSFListDeadNodes = zSFListDeadNodes->next;    
			zFree(n);
		}
	}
}

zSFListNode* zSFListAcquireNode(){
	zSFListNode* temp;
	if(zSFListDeadNodes == NULL){
		temp = zMalloc(sizeof(zSFListNode),"zSFListAcquireNode head");
	}
	else{
		temp = zSFListDeadNodes;
		zSFListDeadNodes = zSFListDeadNodes->next;
	}		
	zClearSfeature(&temp->data);
	temp->next = NULL;
	temp->prev = NULL;
	return temp;
}

void zSFListReleaseNode(zSFListNode* n){
	n->next = zSFListDeadNodes;
	n->prev = NULL;
	zSFListDeadNodes = n;
}

void zInitSFList(zSFList* l){
	zIncrementSFListPool();
	l->head = zSFListAcquireNode();
	l->tail = zSFListAcquireNode();
	l->current = l->head;
	l->head->prev = NULL;
	l->head->next = l->tail;
	l->tail->prev = l->head;
	l->tail->next = NULL;
	l->size = 0;
}

void zFreeSFList(zSFList* l){
	l->current = l->head->next;
	while(l->current != l->tail){
		l->current = l->current->next;
		zSFListReleaseNode(l->current->prev);
	}
	zSFListReleaseNode(l->head);
	zSFListReleaseNode(l->tail);
	l->size = 0;
    zDecrementSFListPool();
}

void zResetSFList(zSFList* l){
	l->current = l->head->next;
	while(l->current != l->tail){
		l->current = l->current->next;
		zSFListReleaseNode(l->current->prev);
	}
	l->head->next = l->tail;
	l->tail->prev = l->head;
    l->current = l->head;
	l->size = 0;
}

zSfeature* zSFListMoveFirst(zSFList* l){
	l->current = l->head->next;
	if(l->current == l->tail){
		return NULL;
	}
	return &l->current->data;
}

zSfeature* zSFListMoveLast(zSFList* l){
	l->current = l->tail->prev;
	if(l->current == l->head){
		return NULL;
	}
	return &l->current->data;
}

zSfeature* zSFListMoveNext(zSFList* l){
    if(l->current != l->tail){
		l->current = l->current->next;
		if(l->current != l->tail){
			return &l->current->data;
		}
	}
	return NULL;
}

zSfeature* zSFListMovePrev(zSFList* l){
    if(l->current != l->head){
		l->current = l->current->prev;
		if(l->current != l->head){
			return &l->current->data;
		}
	}
	return NULL;
}

zSfeature* zSFListGetCurrent(zSFList* l){
	if((l->current == l->head) ||
	   (l->current == l->tail)){
		return NULL;
	}
	else{
		return &l->current->data;
	}
}

zSfeature* zSFListAppend(zSFList* l){
	zSFListNode* n;
	n = zSFListAcquireNode();
	n->next = l->tail;
	n->prev = l->tail->prev;
	l->tail->prev->next = n;
	l->tail->prev = n;
	l->size++;
	return &n->data;
}

zSfeature* zSFListPrepend(zSFList* l){
	zSFListNode* n;
	n = zSFListAcquireNode();
	n->next = l->head->next;
	n->prev = l->head;
	l->head->next->prev = n;
	l->head->next = n;
	l->size++;
	return &n->data;
}

zSfeature* zSFListInsert(zSFList* l,zSfeature* f){
	zSFListNode* n;
	n = zSFListAcquireNode();
	zCopySfeature(f,&n->data);
	n->next = l->tail;
	n->prev = l->tail->prev;
	l->tail->prev->next = n;
	l->tail->prev = n;
	l->size++;
	return &n->data;
}

zSfeature* zSFListInsertNext(zSFList* l,zSfeature* f){
	zSFListNode* n;
	n = zSFListAcquireNode();
	zCopySfeature(f,&n->data);
	if(l->current == l->tail){
		return zSFListInsertPrev(l,f);
	}
	n->next = l->current->next;
	n->prev = l->current;
	l->current->next->prev = n;
	l->current->next = n;
	l->size++;
	return &n->data;
}

zSfeature* zSFListInsertPrev(zSFList* l,zSfeature* f){
	zSFListNode* n;
	n = zSFListAcquireNode();
	zCopySfeature(f,&n->data);
	if(l->current == l->head){
		return zSFListInsertNext(l,f);
	}
	n->next = l->current;
	n->prev = l->current->prev;
	l->current->prev->next = n;
	l->current->prev = n;
	l->size++;
	return &n->data;
}

void zSFListRemoveLast(zSFList* l){
    zSFListNode* n = l->tail->prev;
	if(n != l->head){
		n->prev->next = l->tail;
		l->tail->prev = n->prev;
		zSFListReleaseNode(n);
	}
	l->size--;
}

void zSFListRemoveFirst(zSFList* l){
    zSFListNode* n = l->head->next;
	if(n != l->tail){
		n->next->prev = l->head;
		l->head->next = n->next;
		zSFListReleaseNode(n);
	}
	l->size--;
}

void zSFListRemoveCurrent(zSFList* l){
	zSFListNode* n = l->current;
	l->current = n->next;
	if((n == l->head) || (n == l->tail)){ 
		return;
	}
	n->prev->next = n->next;
	n->next->prev = n->prev;
	zSFListReleaseNode(n);
	l->size--;
}

void zSFList2SFVec(zSFList* l,zSFVec* v){
	zSFListNode* n = l->head->next;
	while(n != l->tail){
		zPushSFVec(v,&n->data);
		n = n->next;
	}
}

#endif
