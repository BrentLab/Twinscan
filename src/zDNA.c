/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zDNA.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_DNA_C
#define ZOE_DNA_C

#include "zDNA.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define Z_DEFAULT_GC_LEVEL 0.50f

static coor_t BLOCK_SIZE = 1000000;
static int BLOCK_COUNT = 2;
static bool BLOCK_FIXED = false;

static bool MAP_READY = 0;
static char S5MAP[256];
static char S16MAP[256];

int zSetDNABlockSize(coor_t size){
	if(BLOCK_FIXED == true){
		return 0;
	}
	BLOCK_SIZE = size;
	return 1;
}

int zSetDNABlockCount(int count){
	if(BLOCK_FIXED == true){
		return 0;
	}
	BLOCK_COUNT = count;
	return 1;
}

static void load_maps(){
	int i;
	for(i = 0;i < 256;i++){
		S5MAP[i] = -1;
		S16MAP[i] = -1;		
	}
	S5MAP['A'] = 0;
	S5MAP['a'] = 0;
	S5MAP['C'] = 1;
	S5MAP['c'] = 1;
	S5MAP['G'] = 2;
	S5MAP['g'] = 2;
	S5MAP['T'] = 3;
	S5MAP['t'] = 3;
	S5MAP['R'] = 4;
	S5MAP['r'] = 4;
	S5MAP['Y'] = 4;
	S5MAP['y'] = 4;
	S5MAP['M'] = 4;
	S5MAP['m'] = 4;
	S5MAP['K'] = 4;
	S5MAP['k'] = 4;
	S5MAP['W'] = 4;
	S5MAP['w'] = 4;
	S5MAP['S'] = 4;
	S5MAP['s'] = 4;
	S5MAP['B'] = 4;
	S5MAP['b'] = 4;
	S5MAP['D'] = 4;
	S5MAP['d'] = 4;
	S5MAP['H'] = 4;
	S5MAP['h'] = 4;
	S5MAP['V'] = 4;
	S5MAP['v'] = 4;
	S5MAP['N'] = 4;
	S5MAP['n'] = 4;

	S16MAP['A'] = 8;
	S16MAP['a'] = 8;
	S16MAP['C'] = 4;
	S16MAP['c'] = 4;
	S16MAP['G'] = 2;
	S16MAP['g'] = 2;
	S16MAP['T'] = 1;
	S16MAP['t'] = 1;
	S16MAP['R'] = 10;
	S16MAP['r'] = 10;
	S16MAP['Y'] = 5;
	S16MAP['y'] = 5;
	S16MAP['M'] = 12;
	S16MAP['m'] = 12;
	S16MAP['K'] = 3;
	S16MAP['k'] = 3;
	S16MAP['W'] = 9;
	S16MAP['w'] = 9;
	S16MAP['S'] = 6;
	S16MAP['s'] = 6;
	S16MAP['B'] = 7;
	S16MAP['b'] = 7;
	S16MAP['D'] = 11;
	S16MAP['d'] = 11;
	S16MAP['H'] = 13;
	S16MAP['h'] = 13;
	S16MAP['V'] = 14;
	S16MAP['v'] = 14;
	S16MAP['N'] = 15;
	S16MAP['n'] = 15;
	
	MAP_READY = 1;		
}

void* zCreateDNASeqBlock(){
	return zCreateSeqBlockWithSize(BLOCK_SIZE);
}

static int zIsGC(char c) {
	return (c == 'c' || c == 'C' || c == 'G' || c == 'g');
}

static int zIsAT(char c) {
	return (c == 'a' || c == 'A' || c == 'T' || c == 't');
}

static int zIsNotMasked(char c) {
	return (c == 'a' || c == 'A' || c == 'c' || c == 'C' || 
			c == 'g' || c == 'G' || c == 't' || c == 'T');
}

static char zComplementDNAChar(char c){
	
	switch (c) {
	case 'A': return 'T';
	case 'C': return 'G';
	case 'G': return 'C';
	case 'T': return 'A';
	case 'R': return 'Y';
	case 'Y': return 'R';
	case 'M': return 'K';
	case 'K': return 'M';
	case 'W': return 'S';
	case 'S': return 'W';
	case 'B': return 'V';
	case 'D': return 'H';
	case 'H': return 'D';
	case 'V': return 'B';
	case 'a': return 't';
	case 'c': return 'g';
	case 'g': return 'c';
	case 't': return 'a';
	case 'r': return 'y';
	case 'y': return 'r';
	case 'm': return 'k';
	case 'k': return 'm';
	case 'w': return 's';
	case 's': return 'w';
	case 'b': return 'v';
	case 'd': return 'h';
	case 'h': return 'd';
	case 'v': return 'b';
	}
	return c;
}	

static void zComplementDNABlock(zSeqBlock* block){
	coor_t i;
	
	for (i = 0; i < block->size; i++) {
		block->data[i] = zComplementDNAChar(block->data[i]);
	}	
}

void zComplementDNA (zDNA *dna) {
	int i,j;
	zSeqBlock *block;
	if(dna->complement == true){
		dna->complement = false;
	}
	else{
		dna->complement = true;
	}
	block = zListMoveFirst(&dna->seq->seq);
	while(block != NULL){
		zComplementDNABlock(block);
		block = zListMoveNext(&dna->seq->seq);
	}
	for(i = 0;i < dna->seq->var_count;i++){
		for(j = 0;j < dna->seq->variants[i]->variants;j++){
			dna->seq->variants[i]->values[j] = 
				zComplementDNAChar(dna->seq->variants[i]->values[j]);
		}
	}
}

void zReverseDNA (zDNA *dna) {
	zReverseSequence(dna->seq);
}

void zInitDNA(zDNA* dna){
	if(MAP_READY == 0){
		load_maps();
	}
	BLOCK_FIXED = true;
	dna->complement = false;
	dna->gc = -1;
	dna->gcs = NULL;
	dna->length = 0;
	dna->seq = NULL;
	dna->def = NULL;
}

void zInitGenericDNA(void *dna) {
	zInitDNA((zDNA*) dna);
}


long zInitDNASequence(zSequence* seq, void *parent){
	coor_t length,gc,unmasked;   
	int def_size;
	char c;
	char defline[4096];     /* should be big enough */
	int map_index;
	int more_records = 0;
	zDNA* dna = (zDNA*)parent;
	seq->parent = parent;
	dna->seq = seq;
	dna->gcs = zMalloc(sizeof(float)*seq->block_count,"zInitDNASequence gcs");
	
	/* initial check for fasta format */
	c = fgetc(seq->fp);
	if (c == EOF) return 0;
	if (c != '>') {
		zWarn("zInitDNASequence \">\" not found");
		return 0;
	}
	(void)ungetc(c, seq->fp);
	
	/* read the def line */
	(void)fgets(defline, sizeof(defline), seq->fp);
	def_size = strlen(defline);
	dna->def = zMalloc(def_size +1,"zInitDNASequence def");
	(void)strcpy(dna->def, defline);
	dna->def[strlen(dna->def) -1] = '\0'; /* remove newline */
	
	length = 0;
	gc = 0;
	unmasked = 0;
	map_index = 0;
	while ((c = fgetc(seq->fp)) != EOF) {
		if (c == '>') {
			(void)ungetc(c, seq->fp);
			more_records = 1;
			break; /* next record found */
		}
		if (isspace((int)c)) continue; /* skip spaces */
		
		if (zIsNotMasked(c)){
			unmasked++;
			if (zIsGC(c)){
				gc++;
			}
		}
		
		if(length == seq->file_map[map_index].seq_pos){
			ungetc(c,seq->fp);
			seq->file_map[map_index].file_pos = ftell(seq->fp);
			c = fgetc(seq->fp);
			map_index++;
		}
		length++;
	}
	
	dna->length = length;	
	seq->real_length = length;
	seq->default_char = 'N';
	
	if(unmasked == 0){
		dna->gc = 0;
	}
	else{
		dna->gc = (float)gc/(float)unmasked;
	}
	dna->complement = false;

	return (more_records == 1)?ftell(seq->fp):more_records;
}

int zReadDNASequence(zSequence* seq, zSeqBlock* block){

	char c;
	coor_t index;
	coor_t gc = 0;
	coor_t unmasked = 0;
	zDNA* dna = (zDNA*)seq->parent;

	if (ferror(seq->fp)) {
		zDie("zReadDNASequence file error 0");
	}

	index = 0;
	while ((c = fgetc(seq->fp)) != EOF) {
		if (c == '>') {
			(void)ungetc(c, seq->fp);
			break; /* next record found */
		}
		if (isspace((int)c)) continue; /* skip spaces */
		
		if (zIsNotMasked(c)){
			unmasked++;
			if (zIsGC(c)){
				gc++;
			}
		}
		
		block->data[index] = toupper(c);
		
		index++;

		if (ferror(seq->fp)) {
			zDie("zReadDNASequence file error 2");
		}

		if(index == seq->block_size){
			break;
		}
	}

	dna->gcs[block->id] = (float)gc/unmasked;
	
	/*fprintf(stderr,"DNA(%d) Load pos %d(%d): %d/%d = %.4f (%.4f total)\n",dna->id,block->pos,block->id,gc,unmasked,dna->gcs[block->id],dna->gc);*/

	while(index != seq->block_size){
		block->data[index] = 'N';
		index++;
	}
	
	/* last check for errors */
	if (ferror(seq->fp)) {
		zDie("zReadDNASequence file error");
	}

	if(dna->complement == true){
		zComplementDNABlock(block);
	}
	
	return 1;
}

void zLoadDNAFromFasta(zDNA* dna, char* filename, char* snp_filename){
	zInitSequenceSpecialized(filename,dna,zInitDNASequence,zReadDNASequence,
							 zCreateDNASeqBlock,BLOCK_SIZE,BLOCK_COUNT,0,-1);
	dna->seq->default_char = 'N';
	if(snp_filename != NULL){
		zLoadSeqVariants(dna->seq,snp_filename);
	}
}

int zLoadMultiDNAFromMultiFasta(zVec *multi_dna, char* filename, char* snp_filename){
	if(snp_filename != NULL){
		zDie("Sequence variation cannot be handled for multifasta files yet");
	}
	return zInitMultiSequenceSpecialized(filename, multi_dna, zInitGenericDNA, sizeof(zDNA), 
					zInitDNASequence, zReadDNASequence, zCreateDNASeqBlock,
					BLOCK_SIZE, BLOCK_COUNT);
}

void zFreeDNA(zDNA *dna) {
	if(dna->seq != NULL){
		zFreeSequence(dna->seq);
		zFree(dna->seq);
		dna->seq = NULL;
	}
	if(dna->gcs != NULL) {
		zFree(dna->gcs);
		dna->gcs = NULL;
	}
	if(dna->def != NULL){
		zFree(dna->def);
		dna->def = NULL;
	}
}

void zCopyDNA(zDNA *dna, zDNA *copy) {
	copy->length = dna->length;
	copy->gc = dna->gc;
	copy->complement = dna->complement;
	copy->gcs = zMalloc(sizeof(float)*dna->seq->block_count,"zCopyDNA gcs");
	memcpy(copy->gcs,dna->gcs,sizeof(float)*dna->seq->block_count);
	copy->def = zMalloc(strlen(dna->def)+1,"zCopyDNA def");
	(void)strcpy(copy->def, dna->def);

	copy->seq = zMalloc(sizeof(zSequence), "zCopyDNA: copy->seq");
	zCopySequence(dna->seq,copy->seq);
	copy->seq->parent = copy;
}

void zAntiDNA (zDNA *dna) {
	zReverseDNA(dna);
	zComplementDNA(dna);
}

void zSetDNAPadding(zDNA* dna,coor_t padding){
	zSetSequencePadding(dna->seq,padding);
	dna->length = dna->seq->length;
}

char zGetDNASeq(zDNA* dna, coor_t pos){
	return zGetSequencePos(dna->seq,pos);
}

char* zGetDNASeqRange(zDNA* dna, coor_t from, coor_t to){
	char* buffer;
	coor_t size = to - from + 1;
	coor_t i;
	buffer = zMalloc(sizeof(char)*(size + 1), "zGetDNASeqRange: buffer");
	for (i = 0; i < size; i++) {
		buffer[i] = zGetDNASeq(dna, from + i);
	}
	buffer[size] = '\0';
	return buffer;
}

char zGetDNAUCSeq(zDNA* dna, coor_t pos){
	char c = zGetSequencePos(dna->seq,pos);
	return toupper(c);
}

char zGetDNAS5(zDNA* dna, coor_t pos){
	char c = zGetSequencePos(dna->seq,pos);
	return S5MAP[(int)c];
}

char zGetDNAS16(zDNA* dna, coor_t pos){
	char c = zGetSequencePos(dna->seq,pos);
	return S16MAP[(int)c];
}

float zGetDNAGC(zDNA* dna){
	return dna->gc;
}

float zGetDNAWindowedGC(zDNA* dna, coor_t pos){
	int i = zGetSeqBlockID(dna->seq,pos);
	if(i == -1){
		zGetSequencePos(dna->seq,pos);
		i = zGetSeqBlockID(dna->seq,pos);
	}
	return dna->gcs[i];
}

bool zDNACheckA(zDNA* dna, coor_t pos){
	return zSeqCheckChar(dna->seq,pos,'A');
}

bool zDNACheckC(zDNA* dna, coor_t pos){
	return zSeqCheckChar(dna->seq,pos,'C');
}

bool zDNACheckG(zDNA* dna, coor_t pos){
	return zSeqCheckChar(dna->seq,pos,'G');
}

bool zDNACheckT(zDNA* dna, coor_t pos){
	return zSeqCheckChar(dna->seq,pos,'T');
}

void zDNASetSNPs(zDNA* dna){ 
	zSetSequenceVariants(dna->seq); 
}

void zAddDNACompanion(zDNA* dna, const char* key, void* ptr) {
	void* old = NULL;
	if ((old = zGetHash(dna->companion, key)) != NULL) {
		zFree(old);
	}
	zSetHash(dna->companion, key, ptr);
}

void *zGetDNACompanion(zDNA* dna, const char* key) {
	return zGetHash(dna->companion, key);
}

float* zCalcSlidingGCLevel(zDNA* dna, size_t win) {
	coor_t left_bound = 0;
	coor_t i = 0;
	coor_t right_bound = 0;

	size_t gc_count = 0;
	size_t at_count = 0;
	size_t gc_total = 0;
	size_t at_total = 0;

	float level;
	float *gcs;

	size_t hwin = win / 2; 
	win = 2 * hwin + 1;
	
	/* Allocate and clear */
	gcs = zMalloc(dna->length * sizeof(float), "zCalcWindowedGCLevel");
	/* I use MIN_SCORE instead of 0.0 so that it's clear where the bugs are */
	for (i = 0; i < dna->length; i++) gcs[i] = -42.00f;
	i = 0;

	if (dna->length < win) {
		win = dna->length; /* don't read past the end of the sequence */
	}

	/* 5' section */
	for ( ; right_bound < win; right_bound++) {
		if (zIsGC(zGetDNAS5(dna,right_bound))) { ++gc_count; ++gc_total; }
		if (zIsAT(zGetDNAS5(dna,right_bound))) { ++at_count; ++at_total; }
	}
	level = (at_count+gc_count == 0) 
		? Z_DEFAULT_GC_LEVEL 
		: (float)gc_count / (at_count+gc_count);
	for ( ; i <= hwin && i < dna->length; i++) {
		gcs[i] = level;
	}

	/* Middle section (this is skipped if win > dna->length) */
	for ( ; right_bound < dna->length; ++right_bound, ++left_bound, ++i) {
		if (zIsGC(zGetDNAS5(dna,right_bound))) { ++gc_count; ++gc_total; }
		if (zIsAT(zGetDNAS5(dna,right_bound))) { ++at_count; ++at_total; }

		if (zIsGC(zGetDNAS5(dna,left_bound)))  --gc_count;
		if (zIsAT(zGetDNAS5(dna,left_bound)))  --at_count;

		level = (at_count+gc_count == 0) 
			? Z_DEFAULT_GC_LEVEL 
			: (float)gc_count / (at_count+gc_count);

		gcs[i] = level;
	}

	/* 3' tail */
	for ( ; i < dna->length; ++i) {
		gcs[i] = level; /* level is the same as the last chunk */
	}

	/* set the DNA's total gc level */
	dna->gc = (at_total + gc_total)
		? (float)gc_total / (float)(at_total + gc_total)
		: Z_DEFAULT_GC_LEVEL;

	/* return the result */
	return gcs;
}

#endif

