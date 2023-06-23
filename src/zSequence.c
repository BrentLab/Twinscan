/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zSequence.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002

\******************************************************************************/

#ifndef ZOE_SEQUENCE_C
#define ZOE_SEQUENCE_C

#include "zSequence.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

static const int BLOCK_SIZE = 100000;
static const int BLOCK_COUNT = 5;

void* zCreateSeqBlock(){
	return zCreateSeqBlockWithSize(BLOCK_SIZE);
}

void* zCreateSeqBlockWithSize(coor_t block_size){
	zSeqBlock* block = zMalloc(sizeof(zSeqBlock),"zCreateSeqBlock block");
	block->size = block_size;
	block->pos = -1;
	block->id = -1;
	block->map_idx = -1;
	block->data = zMalloc(block->size*sizeof(char),"zCreateSeqBlock block");
	return (void*)block;
}

void zFreeSeqBlock(void* v){
	zSeqBlock* block = (zSeqBlock*)v;
	block->pos = -1;
	zFree(block->data);
	zFree(block);
}

void zResetSeqBlock(void* v){
	zSeqBlock* block = (zSeqBlock*)v;
	block->map_idx = -1;
	block->pos = -1;
}

void zCopySeqBlock(zSeqBlock* orig, zSeqBlock* copy){
	copy->id = orig->id;
	copy->pos = orig->pos;
	copy->size = orig->size;
	copy->map_idx = orig->map_idx;
	memcpy(copy->data, orig->data, orig->size*sizeof(char));
}

static zSeqBlock* zLoadSeq(zSequence* seq, coor_t pos){
	int i;
	zSeqBlock* block;
	
	if (ferror(seq->fp)) {
		zDie("zLoadSeq file error 0");
	}
	
	for(i = 0;i < seq->map_size;i++){
		if(pos < seq->file_map[i].seq_pos + seq->block_size){
			break;
		}
	}
	fseek(seq->fp,seq->file_map[i].file_pos,SEEK_SET);
	if (ferror(seq->fp)) {
		zDie("zLoadSeq file error 1");
	}
	block = zListMoveLast(&seq->seq);
	zListMoveCurrentToFirst(&seq->seq);
	if(block->map_idx != -1){
		seq->file_map[block->map_idx].block = NULL;
	}
	block->pos = seq->file_map[i].seq_pos;
	block->map_idx = i;
	seq->file_map[i].block = block;

	/* reads from the current position until EOF or block is full */
	seq->read_func(seq,block);

	/* load the current seq variants */
	zSetSequenceBlockVariants(seq,block);

	/* last check for errors */
	if (ferror(seq->fp)) {
		zDie("zLoadSeq file error");
	}

	return block;
}

int zInitSequence (char* filename, void* parent,
		   zSequenceInitFunc init_func, zSequenceReadFunc read_func){
	return zInitSequenceSpecialized(filename,parent,init_func,read_func,NULL,0,0,0,-1);
}

long zInitSequenceSpecialized (char* filename, void* parent,
				  zSequenceInitFunc init_func, zSequenceReadFunc read_func,
				  zListInitFunc init_block_func, coor_t block_size,
				  coor_t block_count, long fpos_start, long fpos_end){
	zSequence *seq = zMalloc(sizeof(zSequence), "zInitSequenceSpecialized: seq");
	int i, more_records;
	struct stat stbuf;
	zSeqBlock *block;

	if(init_block_func == 0){
		block_size = BLOCK_SIZE;
		block_count = BLOCK_COUNT;
		init_block_func = zCreateSeqBlock;
	}

	seq->var_count   = 0;
	seq->variants    = NULL;
	seq->hap_count   = 0;
	seq->haps        = NULL;
	seq->block_size  = block_size; 
	seq->block_count = block_count;
	seq->parent      = parent;
	seq->init_func   = init_func;
	seq->read_func   = read_func;
	seq->init_block_func = init_block_func;

	seq->reverse     = false;
	seq->padding     = 0;

	seq->filename = zMalloc(strlen(filename) +1,"zInitSeq filename");
	strcpy(seq->filename, filename);
	seq->filename[strlen(seq->filename)] = '\0'; /* remove newline */

	stat(seq->filename,&stbuf);
	/* This size may be too large as it counts the header and the newline chars 
	   but it also saves an extra pass through the file */
	if (fpos_end != -1) {
		seq->map_size = (fpos_end - fpos_start + 1)/seq->block_size + 1;
	} else {
		seq->map_size = stbuf.st_size/seq->block_size + 1;
	}
	if(seq->map_size < seq->block_count){
		seq->map_size = seq->block_count;
	}
	seq->file_map = zMalloc(sizeof(zSeqFileMap)*seq->map_size, "zInitSeq seq_map");

	for(i = 0;i < seq->map_size;i++){
		seq->file_map[i].seq_pos = i*seq->block_size;
		seq->file_map[i].file_pos = -1;
		seq->file_map[i].vars = NULL;
		seq->file_map[i].var_count = 0;
		seq->file_map[i].block = NULL;
	}
	
	if((seq->fp = fopen(seq->filename,"r")) == 0){
		zDie("zInitSeq: Couldn't open sequence file, %s\n",seq->filename);
	}

	/* this function should read the sequence, strip the header, and  and build a map */
	fseek(seq->fp,fpos_start,SEEK_SET);
	more_records = init_func(seq, parent);
	

	zInitList(&seq->seq,seq->init_block_func,zFreeSeqBlock,zResetSeqBlock);
	for(i = 0;i < seq->block_count;i++){
		if(seq->file_map[i].file_pos == -1){
			break;
		}
		block = zListAddLast(&seq->seq);
		block->pos = seq->file_map[i].seq_pos;
		block->id = i;
		block->map_idx = -1;
		seq->file_map[i].block = block;
		zLoadSeq(seq,i*seq->block_size);
	}
	
	if (ferror(seq->fp)) {
		zWarn("zInitSeq ferror");
		return 0;
	}

	if (seq->real_length == 0) {
		zWarn("zInitSeq no sequence");
		return 0;
	}
	seq->length = seq->real_length;

	return more_records;
}

/*********************************************************************** 

Remember to free the vector multi_parent and its contents too

	for (i = 0; i < multi_parent->size; i++) {
		zFreeObject((zObject*)multi_parent->elem[i]);
		zFree(multi_parent->elem[i]);
	}
************************************************************************/
int zInitMultiSequenceSpecialized (char* filename, zVec* multi_parent,
					zGenericInitFunc init_parent_func,
					size_t parent_size,
					zSequenceInitFunc init_sequence_func, 
					zSequenceReadFunc read_func,
					zListInitFunc init_block_func, coor_t block_size,
					coor_t block_count){
	void *parent;
	FILE *fp;
	long fpos_start;
	long more_records;

	if ((fp = fopen(filename, "r")) == NULL) {
		zDie("Cannot open Sequence file %s", filename);
	}


	fseek(fp,0,SEEK_SET);
	fpos_start = 0;
	more_records = 1;
	while (more_records != 0) {
		parent = zMalloc(parent_size, "zInitMultiSequenceSpecialized: parent");
		init_parent_func(parent);
		more_records = zInitSequenceSpecialized(filename,parent,init_sequence_func,read_func,init_block_func,block_size,block_count, fpos_start, -1);
		zPushVec(multi_parent, parent);
		fpos_start = more_records; 
	}
	fclose(fp);
	return multi_parent->size;
}

void zFreeSequence (zSequence* seq){
	int i;
	if(seq->fp != NULL && ferror(seq->fp) == 0){
		fclose(seq->fp);
	}
	for(i = 0;i < seq->map_size;i++){
		zFree(seq->file_map[i].vars);
	}
	zFree(seq->file_map);
	seq->file_map = NULL;
	if(seq->var_count > 0){
		for(i = 0;i < seq->var_count;i++){
			zFree(seq->variants[i]->values);
			zFree(seq->variants[i]->freqs);
			zFree(seq->variants[i]);
		}
		zFree(seq->variants);
	}
	seq->variants = NULL;
	if(seq->hap_count > 0){
		for(i = 0;i < seq->hap_count;i++){
			zFree(seq->haps[i].snp_vals);
		}
		zFree(seq->haps);
	}
	seq->haps = NULL;
	zFreeList(&seq->seq);
	zFree(seq->filename);
	seq->filename = NULL;
}

void zCopySeqVariant(zSeqVariant* orig, zSeqVariant* copy){
	copy->pos = orig->pos;
	copy->real_pos = orig->real_pos;
	copy->variants = orig->variants;
	copy->cur = orig->cur;
	copy->map_idx = orig->map_idx;
	copy->values = zMalloc(sizeof(char)*copy->variants,
						   "zCopySeqVariants copy->values");
	copy->freqs = zMalloc(sizeof(float)*copy->variants,
						   "zCopySeqVariants copy->values");
	memcpy(copy->values, orig->values, orig->variants*sizeof(char));
	memcpy(copy->freqs, orig->freqs, orig->variants*sizeof(float));
}

void zCopySeqHaplotype(zSeqHaplotype* orig,zSeqHaplotype* copy){
	int i;
	copy->first_snp = orig->first_snp;
	copy->last_snp = orig->last_snp;
	copy->snp_count = orig->snp_count;
	copy->snp_vals = zMalloc(sizeof(short)*copy->snp_count,
							 "zCopySeqHaplotype copy->snp_vals");
	for(i = 0; i < copy->snp_count;i++){
		copy->snp_vals[i] = orig->snp_vals[i];
	}
}

void zCopySequence (zSequence* orig, zSequence* copy){

	zSeqBlock *block, *new_block;
	int i,j,k;

	copy->length      = orig->length;
	copy->real_length = orig->real_length;
	copy->map_size    = orig->map_size;
	copy->reverse     = orig->reverse;
	copy->padding     = orig->padding;

	copy->init_func   = orig->init_func;
	copy->read_func   = orig->read_func;
	copy->init_block_func = orig->init_block_func;

	copy->block_size  = orig->block_size;
	copy->block_count  = orig->block_count;
	copy->default_char = orig->default_char;

	zInitList(&copy->seq,copy->init_block_func,zFreeSeqBlock,zResetSeqBlock);
	copy->file_map = zMalloc(sizeof(zSeqFileMap)*copy->map_size,
							 "zCopySeq seq_map");

	copy->var_count = orig->var_count;
	if(orig->variants != NULL){
		copy->variants = zMalloc(sizeof(zSeqVariant*)*orig->var_count,
							"zCopySequence copy->variants");
		for(i = 0;i < copy->var_count;i++){
			copy->variants[i] = zMalloc(sizeof(zSeqVariant),
							"zCopySequence copy->variants[i]");
			zCopySeqVariant(orig->variants[i],copy->variants[i]);
		}
	} else {
		copy->variants = NULL;
	}
	copy->hap_count = orig->hap_count;
	if(orig->haps != NULL){
		copy->haps = zMalloc(sizeof(zSeqHaplotype)*orig->hap_count,
							"zCopySequence copy->haps");
		for(i = 0;i < copy->hap_count;i++){
			zCopySeqHaplotype(&orig->haps[i],&copy->haps[i]);
		}
	} else {
		copy->haps = NULL;
	}
	
	k = 0;
	for(i = 0;i < copy->map_size;i++){
		copy->file_map[i].seq_pos = orig->file_map[i].seq_pos;
		copy->file_map[i].file_pos = orig->file_map[i].file_pos;
		copy->file_map[i].block = NULL;
		copy->file_map[i].var_count = orig->file_map[i].var_count;
		copy->file_map[i].vars = zMalloc(sizeof(zSeqVariant*)*copy->file_map[i].var_count,
							"zCopySequence copy->file_map[i]");
		for(j = 0;j < copy->file_map[i].var_count; j++){
			copy->file_map[i].vars[j] = copy->variants[k++];
		}
	}
	copy->filename = zMalloc(strlen(orig->filename)+1,"zCopy filename");
	strcpy(copy->filename, orig->filename);
	copy->fp = fopen(copy->filename,"r");

	block = zListMoveFirst(&orig->seq);
	while(block != NULL){
		new_block = zListAddLast(&copy->seq);
		zCopySeqBlock(block,new_block);
		copy->file_map[new_block->map_idx].block = new_block;
		block = zListMoveNext(&orig->seq);
	}
}

static zSeqBlock* zGetSeqBlock(zSequence* seq, coor_t pos){
	zSeqBlock* block;
	block = zListMoveFirst(&seq->seq);
	while(block != NULL){
		if((block->pos <= pos) && (block->pos + block->size > pos)){
			zListMoveCurrentToFirst(&seq->seq);
			return block;
		}
		block = zListMoveNext(&seq->seq);
	}
	return zLoadSeq(seq,pos);
}

int zGetSeqBlockID(zSequence* seq, coor_t pos){
	zSeqBlock* block;

	if(seq->reverse == true){
		pos = seq->length-pos-1;
	} 
	if(pos < seq->padding){
		pos = seq->padding;
	}
	else{
		pos -= seq->padding;
	}
	if(pos >= seq->real_length){
		pos = seq->real_length-1;
	}

	block = zListMoveFirst(&seq->seq);
	while(block != NULL){
		if((block->pos <= pos) && (block->pos + block->size > pos)){
			zListMoveCurrentToFirst(&seq->seq);
			return block->id;
		}
		block = zListMoveNext(&seq->seq);
	}
	return -1;
}

char zGetSequencePos(zSequence* seq,coor_t pos){
	zSeqBlock *block;
	coor_t block_pos;
	
	if(seq->reverse == true){
		pos = seq->length-pos-1;
	} 
	if(pos < seq->padding){
		return seq->default_char;
	}
	pos -= seq->padding;
	if(pos >= seq->real_length){
		return seq->default_char;
	}
	
	block = zGetSeqBlock(seq,pos);
	block_pos = pos - block->pos;
	return block->data[block_pos];
}

void zSetSequencePadding(zSequence* seq,coor_t padding){
	int i;
	seq->padding = padding;
	seq->length = seq->real_length + 2*padding;
	for(i = 0;i < seq->var_count;i++){
		seq->variants[i]->pos = seq->variants[i]->real_pos + padding;
	}
}

void zReverseSequence (zSequence *seq) {
	int i;
	if(seq->reverse == true){
		seq->reverse = false;
	}
	else{
		seq->reverse = true;
	}
	for(i = 0;i < seq->var_count;i++){
		seq->variants[i]->pos = seq->length - seq->variants[i]->pos -1;
	}
}

void zLoadSeqVariants(zSequence *seq,char *filename){
	FILE* fp;
	int i,j,map,var_count,next_var,hap_count;
	coor_t last_pos = 0;
	coor_t next_map;
	char c;
	int good;
	int* var_map;
    /* open file */
	if((fp = fopen(filename,"r")) == 0){
		zDie("zLoadSeqVariants: Couldn't open sequence variant file, %s\n",filename);
	}
	/* read total number of varying positions */
	if(!(fscanf(fp,"%d\n",&var_count)))
		zDie("Bad sequence variant file, %s.  Can't get var_count\n",filename); 
	seq->var_count = var_count;
	/* allocate temp storage for zSeqVariant objects */
	seq->variants = zMalloc(sizeof(zSeqVariant*)*var_count,"zLoadSeqVariants seq->variants");
	var_map = zMalloc(sizeof(int)*(var_count+1),"zLoadSeqVariants var_map");
	/* count of varying position for each sequence block */	
	seq->file_map[0].var_count = 0;
	map = 0;
	next_map = seq->block_size;
	/*get haplotype count */
	if(!(fscanf(fp,"%d\n",&hap_count)))
		zDie("Bad sequence variant file, %s.  Can't get hap_count\n",filename);
	seq->hap_count = hap_count;
	/* read each sequence variant */
	i = 0;
	for(next_var = 0;next_var < var_count; next_var++){
		seq->variants[i] = zMalloc(sizeof(zSeqVariant),"zLoadSeqVariant seq->variants[i]");
		seq->variants[i]->freqs = NULL;
		seq->variants[i]->values = NULL;
		seq->variants[i]->cur = 0;
		if(!(fscanf(fp,"%u",&seq->variants[i]->real_pos)))
			zDie("Bad sequence variant file format %s\n",filename);
		seq->variants[i]->real_pos--; /* adjust for 0 based indexing */
		seq->variants[i]->pos = seq->variants[i]->real_pos;
		/* verify this position is on the sequence */
		if(seq->variants[i]->real_pos >= seq->real_length){
			zDie("Bad sequence variant position %ul: beyond end of sequence",
				 seq->variants[i]->real_pos); 
		}
		/* verify this position is in order */
		if(seq->variants[i]->real_pos < last_pos){
			zDie("Bad sequence variant file: out of order");
		}
		last_pos = seq->variants[i]->real_pos;
		/* read the number and value of the values for this variant */
		if(!(fscanf(fp,"%hd",&seq->variants[i]->variants)))
			zDie("Bad sequence variant file format %s\n",filename);
		seq->variants[i]->values = zMalloc(sizeof(char)*(seq->variants[i]->variants+1),
										   "zLoadSeqVariant seq->variants[i].values");
		if(!(fscanf(fp,"%s",seq->variants[i]->values)))
			zDie("Bad sequence variant file format %s\n",filename);
		seq->variants[i]->freqs = zMalloc(sizeof(float)*(seq->variants[i]->variants),
										  "zLoadSeqVariant seq->variants[i].freqs");
		c = zGetSequencePos(seq,seq->variants[i]->pos);
		good = -1;
		for(j = 0;j < seq->variants[i]->variants; j++){
			if(seq->variants[i]->values[j] == c){
				good = j;
			}
			if(!(fscanf(fp,"%f",&seq->variants[i]->freqs[j])))
				zDie("Bad sequence variant file format %s\n",filename);
		}
		if(good == -1){
			/*EVAN this is bad.  zSequence knows nothing about N masking */
			if(c == 'N'){
				zWarn("zSequence: masked vairant at pos %d\n\tseq val: %c\tvar val: %s\n",
					  seq->variants[i]->real_pos+1,c,seq->variants[i]->values);
			}
			else{
				zWarn("zSequence: bad vairant at pos %d\n\tseq val: %c\tvar val: %s\n",
					  seq->variants[i]->real_pos+1,c,seq->variants[i]->values);
			}
			zFree(seq->variants[i]->values);
			zFree(seq->variants[i]->freqs);
			seq->variants[i]->freqs = NULL;
			seq->variants[i]->values = NULL;
			zFree(seq->variants[i]);
			/* +1 changes from 0 to 1 indexing */
			var_map[next_var+1] = -1;
		}
		else{ 
			if(good != 0){
				zWarn("zSequence: variant at pos %d\n\thas default value (%c) which is not the same as the reference (%c)\n",
					  seq->variants[i]->real_pos+1,seq->variants[i]->values[0],c);
			}
			/* if the variant is on the next block update the vars */
			while(seq->variants[i]->real_pos >= next_map){
				map++;
				seq->file_map[map].var_count = 0;
				next_map += seq->block_size;
			}
			/* update the variant count for this block */
			seq->file_map[map].var_count++;
			seq->variants[i]->map_idx = map;
			/* +1 changes from 0 to 1 indexing */
			var_map[next_var+1] = i;
			i++;
		}
	}
	seq->var_count = i; /* only good snps were kept */
	/* allocate and fill the variant arrays */
	j = 0;
	for(map = 0;map < seq->map_size;map++){
		seq->file_map[map].vars = zMalloc(sizeof(zSeqVariant*)*seq->file_map[map].var_count,
										  "zLoadSeqVariant seq->file_map->vars[map]");
		for(i = j;i < j+seq->file_map[map].var_count;i++){
			seq->file_map[map].vars[i-j] = seq->variants[i];
		}
		j += seq->file_map[map].var_count;
	}
	zSetSequenceVariants(seq);
	/* if haplotype data is inlcuded read it*/
	if(seq->hap_count > 0){
		zSeqHaplotype* hap;
		int next_hap, next_hap_idx;
		char* val = zMalloc(sizeof(char)*(var_count+10),"zLoadSeqVariant val");
		seq->haps = zMalloc(sizeof(zSeqHaplotype)*seq->hap_count,"zLoadSeqVariant seq->haps");
		next_hap_idx = 0;
		for(next_hap = 0; next_hap < seq->hap_count; next_hap++){
			int first_snp,last_snp,snp_count,next_snp,next_snp_idx;
			first_snp = -1; /* insure compiler shush */
			last_snp = -1; /* insure compiler shush */
			if(!(fscanf(fp,"%d",&first_snp)))
				zDie("Bad sequence variant file format %s\n",filename);
			if(!(fscanf(fp,"%d",&last_snp)))
				zDie("Bad sequence variant file format %s\n",filename);
			snp_count = 0;
			hap = &seq->haps[next_hap_idx];
			/* count the snps which passed the above filters */
			for(next_snp = first_snp; next_snp <= last_snp; next_snp++){
				if(var_map[next_snp] != -1){
					snp_count++;
				}
			}
			/* if none of the snps passed the above filters then just 
			   skip this haplotype */
			if(snp_count == 0){
				continue;
			}
			/* allocate and fill the vals array */
			hap->snp_count = snp_count;
			hap->first_snp = -1;
			hap->last_snp = -1;
			hap->snp_vals = zMalloc(sizeof(int)*hap->snp_count,
									"zLoadSeqVariant hap->snp_vals");
			next_snp_idx = 0;
			if(!(fscanf(fp,"%s",&(val[2]))))
				zDie("Bad sequence variant file format %s\n",filename);
			for(next_snp = first_snp; next_snp <= last_snp; next_snp++){
				val[1] = '\0';
				/* if the snp passed the filter (-1 to get 0 based index)*/
				if(var_map[next_snp] != -1){
					val[0] = val[next_snp-first_snp+2];
					hap->snp_vals[next_snp_idx] = (short)atoi(val);
					next_snp_idx++;
					hap->last_snp = var_map[next_snp];
				}
				if(hap->first_snp == -1){
					hap->first_snp = var_map[next_snp];
				}
			}
			next_hap_idx++;
		}
	}
	zFree(var_map);
}

void zSetSequenceVariant(zSequence *seq, int var, int val){
	zSeqBlock* block;
	zSeqFileMap* map;
	zSeqVariant* variant;
	/* get the variant */ 
	variant = seq->variants[var];
	/* set the cur value */
	variant->cur = val;
	/* get the file map */
	map = &seq->file_map[variant->map_idx];
	/* get the block */
	block = map->block;
	/* if the sequence is in memory set the value */
	if(block != NULL){
		block->data[variant->real_pos - block->pos] = 
			variant->values[variant->cur];
	}
} 

void zSetSequenceBlockVariants(zSequence *seq,zSeqBlock *block){
	int i;
	zSeqVariant* var;
	zSeqFileMap* map;
	map = &seq->file_map[block->map_idx];
	/* set all variants in map to cur value*/
	for(i = 0;i < map->var_count;i++){
		var = map->vars[i];
		block->data[var->real_pos-block->pos] = var->values[var->cur];
	}	
}

void zSetSequenceVariants(zSequence *seq){
	zSeqBlock* block;
	block = zListMoveFirst(&seq->seq);
	while(block != NULL){
		zSetSequenceBlockVariants(seq,block);
		block = zListMoveNext(&seq->seq);
	}
}

bool zSeqCheckChar(zSequence* seq,coor_t pos,char check){
	zSeqBlock *block;
	zSeqFileMap *map;
	coor_t block_pos;
	int i,j;
	
	if(seq->reverse == true){
		pos = seq->length-pos-1;
	} 
	if(pos < seq->padding){
		return (check == seq->default_char);
	}
	pos -= seq->padding;
	if(pos >= seq->real_length){
		return (check == seq->default_char);
	}
	
	block = zGetSeqBlock(seq,pos);
	block_pos = pos - block->pos;
	if(block->data[block_pos] == check){
		return true;
	}
	i = 0;
	map = &seq->file_map[block->map_idx];
	if(map->var_count == 0){
		return false;
	}
	while(i < map->var_count-1 && map->vars[i]->real_pos < pos){
		i++;
	}
	if(map->vars[i]->real_pos == pos){
		for(j = 0; j < map->vars[i]->variants; j++){
			if(map->vars[i]->values[j] == check){
				return true;
			}
		}
	}
	return false;	
}

#endif



