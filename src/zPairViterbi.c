#include <time.h>
#include "zTBTree.h"
#include "zViterbi.h"
#include "zHardCoding.h"

/**********************************************************************\
  zViterbi

  This hold the data needed for the memory optimized viterbi run

\**********************************************************************/

/* zViterbi structure holds the trellis, traceback tree and tbnodes */ 

static const coor_t  VITERBI_CACHE_LENGTH             = 3;
static const int     MAX_CONCURRENT_SEQUENCE_VARIANTS = 16384; /*2^14*/
static const int     MAX_CONCURRENT_SNPS              = 6;
static const int     MAX_ACTIVE_SNPS                  = 6;
static const float   HAPLOTYPE_FREQUENCY_THRESHOLD    = 0.;

extern coor_t any_new_cpoint;

/* this moves live nodes from links ending at the current pos from the live node list
	 into the current cells.  */
int zViterbiIncorporateLiveNodes(zViterbi* v, int allele){
	int removed_node_count = 0;
	zPtrList* l = v->live_nodes[allele];
	zTBTreeNode** cells = v->cache[allele][0][v->cdna_pos];
	zTBTreeNode* n = NULL;

	n = zPtrListMoveFirst(l);
	/* until the first live node comes after pos or we reach the end of the list */
	while(n != NULL && n->pos == v->pos){
		/* incorporate the live node */
		cells[n->state] = n;
		/* move to the next node */
		n = zPtrListMoveNext(l);
		/* remove the incorporated live node */
		zPtrListRemovePrev(l);
		removed_node_count++;
	}

	return removed_node_count;
}


/* copy allele a1 to a2 including live_nodes, tbtree, and cells */
void zInitAlleleDecode(zAlleleDecode* ad, size_t header_size, size_t val_size){
	ad->sfv = zMalloc(sizeof(zSFVec),"zInitAlleleDecode ad->sfv");
	zInitSFVec(ad->sfv,1);
	ad->sfl = zMalloc(sizeof(zSFList),"zInitAlleleDecode ad->sfl");
	zInitSFList(ad->sfl);

	ad->afv = zMalloc(sizeof(zAFVec),"zInitAlleleDecode ad->afv");
	zInitAFVec(ad->afv,1);
	ad->afl = zMalloc(sizeof(zAFList),"zInitAlleleDecode ad->afl");
	zInitAFList(ad->afl);

	ad->header = zMalloc(sizeof(char)*header_size,"zInitAlleleDecode ad->header");
	ad->vals = zMalloc(sizeof(short)*val_size,"zInitAlleleDecode ad->vals");
}

void zFreeAlleleDecode(zAlleleDecode* ad){
	if(ad->sfv != NULL){
		zFreeSFVec(ad->sfv);
		zFree(ad->sfv);
	}
	if(ad->sfl != NULL){
		zFreeSFList(ad->sfl);
		zFree(ad->sfl);
	}
	if(ad->afv != NULL){
		zFreeAFVec(ad->afv);
		zFree(ad->afv);
	}
	if(ad->afl != NULL){
		zFreeAFList(ad->afl);
		zFree(ad->afl);
	}
	if(ad->header != NULL){
		zFree(ad->header);
		ad->header = NULL;
	}
	if(ad->vals != NULL){
		zFree(ad->vals);
		ad->vals = NULL;
	}
}

void zSNPCleanUpDecodes(zViterbi *v){
	zAlleleDecode *ad;
	zAlleleDecode *main_ad = zPtrListMoveFirst(v->traceback);
	zSfeature* f1;
	zSfeature* f2;	

	f1 = zSFListMoveFirst(main_ad->sfl);
	while(f1 != NULL){
		f1 = zSFListMoveNext(main_ad->sfl);
	}

	
	ad = zPtrListMoveNext(v->traceback);
	/*ad = NULL; EVAN this blocks clean up */
	while(ad != NULL){
		
		f2 = zSFListMoveFirst(ad->sfl);
		while(f2 != NULL){
			f2 = zSFListMoveNext(ad->sfl);
		}

		/*EVAN this stupidly traveres the entire main_ad->sfl for each ad it checks*/
		/* remove same features at start of ad->sfl */
		f1 = zSFListMoveFirst(main_ad->sfl);
		f2 = zSFListMoveFirst(ad->sfl);
		if(f2 == NULL){
			ad = zPtrListMoveNext(v->traceback);
			continue;
		}
		
		while(f1 != NULL && f1->start > f2->end){
			/*printf("check %d < %d\n",f1->start,f2->end);*/
			f1 = zSFListMoveNext(main_ad->sfl);
		}
		/* if last state is internal it only needs to match at the start */
		if((v->hmm->state[f2->state].type == INTERNAL || 
			v->hmm->state[f2->state].type == GINTERNAL) && 
		   (f1 != NULL && f1->start == f2->start && f1->state == f2->state)){
			/*printf("REMOVE LAST %d\n",f2->start);*/
			zSFListRemoveFirst(ad->sfl);
			f2 = zSFListMoveFirst(ad->sfl);
			f1 = zSFListMoveNext(main_ad->sfl);		
		}
		while(f1 != NULL && f2 != NULL && zSfeatureCmp(f1,f2) == 0){
			/*printf("REMOVE LAST %d\n",f2->start);*/
			zSFListRemoveFirst(ad->sfl);
			f2 = zSFListMoveFirst(ad->sfl);
			f1 = zSFListMoveNext(main_ad->sfl);
		}
		/* remove same features at end of ad->sfl */		
		f1 = zSFListMoveLast(main_ad->sfl);
		f2 = zSFListMoveLast(ad->sfl);
		if(f2 == NULL){
			ad = zPtrListMoveNext(v->traceback);
			continue;
		}

		while(f1 != NULL && f1->end < f2->start){
			/*printf("check %d < %d\n",f1->end,f2->start);*/
			f1 = zSFListMovePrev(main_ad->sfl);
		}
		/* if first state is internal it only needs to match at the start */
		if((v->hmm->state[f2->state].type == INTERNAL ||
			v->hmm->state[f2->state].type == GINTERNAL) && 
			 (f1 != NULL && f1->end == f2->end && f1->state == f2->state)){
			/*printf("REMOVE FIRST %d\n",f2->start);*/
			zSFListRemoveLast(ad->sfl);
			f2 = zSFListMoveLast(ad->sfl);
			f1 = zSFListMovePrev(main_ad->sfl);		
		}
		while(f1 != NULL && f2 != NULL && zSfeatureCmp(f1,f2) == 0){
			/*printf("REMOVE FIRST %d\n",f2->start);*/
			zSFListRemoveLast(ad->sfl);
			f2 = zSFListMoveLast(ad->sfl);
			f1 = zSFListMovePrev(main_ad->sfl);
		}
		/* get next ad */
		ad = zPtrListMoveNext(v->traceback);
	}
	
	ad = zPtrListMoveFirst(v->traceback);
	while(ad != NULL){
		
		if(ad != main_ad){
			f1 = zSFListMoveFirst(ad->sfl);
			if(f1 == NULL){
				ad->diff_start = 0;
				ad->diff_stop = 0;
			}
			else{
				ad->diff_stop = f1->end+1;
				f1 = zSFListMoveLast(ad->sfl);
				ad->diff_start = f1->start+1;
			}
			zSNPFillAlleleHeader(v,ad);
		}
		zSFList2SFVec(ad->sfl,ad->sfv);
		ad = zPtrListMoveNext(v->traceback);
	}
}

/**********************************************************************\
  zViterbi / zPairTrellis functions
\**********************************************************************/

static void zInitPairViterbi(zViterbi* v,zPairTrellis* trellis){
	int i,j,k;
	coor_t m;
	zHMM         *hmm       = trellis->hmm;
	zSeqVariant **vars      = trellis->genomic->seq->variants;
	int           var_count = trellis->genomic->seq->var_count;

	v->pair_trellis         = trellis;
	v->hmm                  = hmm;

	v->sfl                  = NULL;

	v->max_trellis_count    = MAX_CONCURRENT_SEQUENCE_VARIANTS;
	v->max_snp_count        = MAX_CONCURRENT_SNPS;
	v->max_active_snp_count = MAX_ACTIVE_SNPS;
	v->hap_threshold        = HAPLOTYPE_FREQUENCY_THRESHOLD;
	v->cache_length         = VITERBI_CACHE_LENGTH;

	v->trellis_count        = 0;
	v->pos                  = 0;
	v->cpos                 = 0;
	v->cdna_pos             = 0;
	v->cdna_min             = 0;
	v->cdna_max             = trellis->cdna->length;

	if (var_count == 0) v->max_trellis_count = 1;

	if((unsigned int)v->max_snp_count > sizeof(int)*8 + 1){
		zDie("zViterbi max_snp_count %d > %d which is max for this architecture\n",
			 v->max_snp_count,sizeof(int)*8 + 1);
	}


	for(i = 0; i < hmm->states; i++){
		if(hmm->state[i].type == EXPLICIT){
			zDurationGroup *dg = hmm->dmap[hmm->state[i].duration];
			int group_index    = zGetDurationIsochoreGroup(dg, trellis->genomic->gc);
			v->cache_length    = MAX(v->cache_length, (int)dg->duration[group_index].distribution[0].end + 2);
		}
	}

	v->snp_count = zMalloc(sizeof(short)*v->max_trellis_count,
						   "zInitViterbi snp_count");
	v->active_snp_count = zMalloc(sizeof(short)*v->max_trellis_count,
								  "zInitViterbi active_snp_count");
	v->first_snp = zMalloc(sizeof(int)*v->max_trellis_count,
						   "zInitViterbi first_snp");
	v->temp = zMalloc(sizeof(int)*v->max_trellis_count,
					  "zInitViterbi temp");
	v->snp_vals = zMalloc(sizeof(short*)*v->max_trellis_count,
						  "zInitViterbi snp_vals");
	v->snp_int_vals = zMalloc(sizeof(int)*v->max_trellis_count,
							  "zInitViterbi snp_int_vals");
	v->tbtrees = zMalloc(sizeof(zTBTree*)*v->max_trellis_count,
						 "zInitViterbi tbtree");
	v->live_nodes = zMalloc(sizeof(zPtrList*)*v->max_trellis_count,
							"zInitViterbi live_nodes");
	v->cache = zMalloc(sizeof(zTBTreeNode****)*v->max_trellis_count,
					   "zInitViterbi cache");  
	v->next_tbtn = zMalloc(sizeof(zTBTreeNode*)*v->max_trellis_count,
						   "zInitViterbi next_tbtn");
	v->next_tbti = zMalloc(sizeof(int)*v->max_trellis_count,
						   "zInitViterbi next_tbti");	
	v->tbt_use_count = zMalloc(sizeof(short)*v->max_trellis_count,
							   "zInitViterbi tbt_use_count");
	v->start_pos = zMalloc(sizeof(coor_t)*v->max_trellis_count,
						   "zInitViterbi start_pos");
	v->active = zMalloc(sizeof(short)*v->max_trellis_count,
						"zInitViterbi active");

	for(i = 0;i < v->max_trellis_count; i++){
		v->snp_count[i] = 0;
		v->active_snp_count[i] = 0;
		v->first_snp[i] = -1;
		v->snp_int_vals[i] = 0;
		v->snp_vals[i] = zMalloc(sizeof(short)*v->max_snp_count,
								 "zInitViterbi snp_vals[i]");
		for(j = 0; j < v->max_snp_count; j++){
			v->snp_vals[i][j] = 0;
		}
		v->tbtrees[i] = zMalloc(sizeof(zTBTree),
								"zInitViterbi tbtrees[i]");
		zInitTBTree(v->tbtrees[i]);
		v->tbtrees[i]->root->pos = trellis->padding - 1;
		v->tbtrees[i]->root->cdna_pos = trellis->padding - 1;
		v->live_nodes[i] = zMalloc(sizeof(zPtrList),
								   "zInitViterbi lide_nodes[i]");
		zInitPtrList(v->live_nodes[i]);
		v->cache[i] = zMalloc(sizeof(zTBTreeNode***)*v->cache_length, "zInitViterbi cache[i]");
		for (k = 0; k < v->cache_length; k++) {
			v->cache[i][k] = zMalloc(sizeof(zTBTreeNode**)*(v->cdna_max - v->cdna_min + 1), "zInitViterbi cache[i][k]"); 
			for (m = v->cdna_min; m <= v->cdna_max; m++) {
				v->cache[i][k][m] = zMalloc(sizeof(zTBTreeNode*)*hmm->states, "zInitViterbi cache[i][k][m]"); 
				for(j = 0; j < hmm->states; j++){
					v->cache[i][k][m][j] = NULL;
				}
			}
		}

		v->next_tbtn[i] = NULL;
		v->next_tbti[i] = 0;
		v->tbt_use_count[i] = 0;
		v->start_pos[i] = (coor_t)-1;
		v->active[i] = 0;
	}
	
	/*EVAN what should we do with the first snp in haplotypes? If only a single value is ever observed 
	  this info is not represented here.  this doesn't matter now since the hapmap.org haplotypes cover
	  the full chromosomes but could potentially matter in the future */
	/* allocate and fill hapmap and hapmax arrays */
	v->hapmap = zMalloc(sizeof(float**)*var_count,
						"zInitViterbi hapmap");
	v->hapmax = zMalloc(sizeof(short*)*var_count,
						"zInitViterbi hapmax");
	for(i = 0; i < var_count; i++){
		int count = 0;
		v->hapmap[i] = zMalloc(sizeof(float*)*vars[i]->variants,
							   "zInitViterbi hapmap[i]");
		v->hapmax[i] = zMalloc(sizeof(short)*vars[i]->variants,
							   "zInitViterbi hapmax[i]");
		if(i == var_count -1){
			for(j = 0;j < vars[i]->variants; j++){
				v->hapmap[i][j] = zMalloc(sizeof(float)*1,
										  "zInitViterbi hapmap[i][j]");
				v->hapmap[i][j][0] = 0;
			}
			continue;
		}
		for(j = 0;j < vars[i]->variants; j++){
			v->hapmap[i][j] = zMalloc(sizeof(float)*vars[i+1]->variants,
									  "zInitViterbi hapmap[i][j]");
			for(k = 0;k < vars[i+1]->variants; k++){
				v->hapmap[i][j][k] = 0;
			}
		}
		for(k = 0;k < trellis->genomic->seq->hap_count; k++){
			zSeqHaplotype* hap;			
			hap = &trellis->genomic->seq->haps[k];
			if(i >= hap->first_snp && i+1 <= hap->last_snp){
				int v1 = hap->snp_vals[i-hap->first_snp];
				int v2 = hap->snp_vals[i+1-hap->first_snp];
				v->hapmap[i][v1][v2]++;
				count++;
			}
		}
		if(count == 0){
			/* if no haplotypes cover this region */
			for(j = 0;j < vars[i]->variants; j++){
				for(k = 0;k < vars[i+1]->variants; k++){
					/*fprintf(stderr,"var %d val %d nextval %d\t%d/%d =~ %.3f\n",
					  i,j,k,(int)v->hapmap[i][j][k],count,(float)1/vars[i+1]->variants);
					EVAN*/
					v->hapmap[i][j][k] = (float)1/vars[i+1]->variants;
				}
				v->hapmax[i][j] = 0;
				/*fprintf(stderr,"hapmax[%d][%d] = %d\n",i,j,0);
				 EVAN*/
			}
		}
		else{
			for(j = 0;j < vars[i]->variants; j++){
				int max_idx = 0;
				int max_val = 0;
				int val_count = 0;
				for(k = 0;k < vars[i+1]->variants; k++){
					if(v->hapmap[i][j][k] > max_val){
						max_val = v->hapmap[i][j][k];
						max_idx = k;
					}
					val_count += v->hapmap[i][j][k];
				}
				v->hapmax[i][j] = max_idx;
				/*fprintf(stderr,"hapmax[%d][%d] = %d\n",i,j,v->hapmax[i][j]);EVAN*/
				if(val_count == 0){
					/*for(k = 0;k < vars[i+1]->variants; k++){
						fprintf(stderr,"var %d val %d nextval %d\t%d/0 = 0\n",
							   i,j,k,(int)v->hapmap[i][j][k]);
							   }EVAN*/
				}
				else{
					for(k = 0;k < vars[i+1]->variants; k++){
						/*fprintf(stderr,"var %d val %d nextval %d\t%d/%d = %.3f\n",
							   i,j,k,(int)v->hapmap[i][j][k],val_count,v->hapmap[i][j][k]/val_count);
							   EVAN*/
						v->hapmap[i][j][k] /= val_count;
					}
				}
			}
		}
	}

	v->afl = zMalloc(sizeof(zAFList),"zInitViterbi afl");
	zInitAFList(v->afl);	
	v->traceback = zMalloc(sizeof(zPtrList),
						   "zInitViterbi traceback");
	zInitPtrList(v->traceback);

	/* don't use PIN stuff */
	v->next = NULL;
	v->dead_count = -1;
	v->active_count = -1;
}

static void zFreePairViterbi(zViterbi* v){
	int i, j, k;
	coor_t m;
	/* free trellis stuff */
	for(i = 0;i < v->max_trellis_count; i++){
		for (k = 0; k < v->cache_length; k++) {
			for (m = v->cdna_min; m <= v->cdna_max; m++) {
				for(j = 0; j < v->pair_trellis->hmm->states; j++){
					v->cache[i][k][m][j] = NULL;
				}
				zFree(v->cache[i][k][m]);
			}
			zFree(v->cache[i][k]);
		}
		zFree(v->cache[i]);
		zFreePtrList(v->live_nodes[i]);
		zFree(v->live_nodes[i]);
		zFreeTBTree(v->tbtrees[i]);
		zFree(v->tbtrees[i]);
	}
	zFree(v->live_nodes);
	zFree(v->cache);
	zFree(v->tbtrees);
	zFree(v->active);
	zFree(v->start_pos);

	zFreeAFList(v->afl);
	zFree(v->afl);

	/* free SNP specific stuff */
	if(v->snp_vals != NULL){
		for(i = 0;i < v->max_trellis_count; i++){
			zFree(v->snp_vals[i]);
		}	
		zFree(v->snp_vals);
		zFree(v->snp_int_vals);
		zFree(v->first_snp);
		zFree(v->snp_count);
		zFree(v->active_snp_count);
		zFree(v->next_tbtn);
		zFree(v->next_tbti);
		zFree(v->tbt_use_count);
		zFree(v->temp);
		for(i = 0;i < v->pair_trellis->genomic->seq->var_count; i++){
			for(j = 0; j < v->pair_trellis->genomic->seq->variants[i]->variants; j++){
				zFree(v->hapmap[i][j]);
			}
			zFree(v->hapmap[i]);
			zFree(v->hapmax[i]);
		}
		zFree(v->hapmap);
		zFree(v->hapmax);
	}
	
	/* free PIN specific stuff */
	if(v->active != NULL){
		zFree(v->next);
	}

	if(v->traceback != NULL){
		zFreePtrList(v->traceback);
		zFree(v->traceback);
	}
}

/* sets all snps in the DNA seq for this allele */
static void zSNPSetSeqForAllele(zViterbi* v,int allele){
	int i,j,f;
	short* vals = v->snp_vals[allele];
	if(v->snp_count[allele] == 0){
		return;
	}
	f = v->first_snp[allele];
	j = 0;
	for(i = 0; i < v->snp_count[allele]; i++){
		v->pair_trellis->genomic->seq->variants[i+f]->cur = vals[j];
		j++;
	}
	zDNASetSNPs(v->pair_trellis->genomic);
}

/* unsets all snps in the DNA seq for this allele */
static void zSNPUnSetSeqForAllele(zViterbi* v,int allele){
	int i,f;
	if(v->snp_count[allele] == 0){
		return;
	}
	f = v->first_snp[allele];
	for(i = 0; i < v->snp_count[allele]; i++){
		v->pair_trellis->genomic->seq->variants[i+f]->cur = 0;
	}
	zDNASetSNPs(v->pair_trellis->genomic);
}

/* copies a zTBTree object into a zAFList tracing back from a zTBTreeNode */
static void zTBTree2AFL(zAFList* afl, zTBTree* tree, zTBTreeNode* in, zHMM* hmm, zPairTrellis* trellis){
	zAlnFeature* f = NULL;
	zTBTreeNode* n = in;

	while(n != tree->root){
		f = zAFListAppend(afl);
		f->name   = hmm->state[n->state].name;
		f->genomic = trellis->genomic;
		f->genomic_end    = n->pos;
		f->genomic_start  = n->parent->pos;
		f->cdna = trellis->cdna;
		f->cdna_end    = n->cdna_pos;
		f->cdna_start  = n->parent->cdna_pos;
		f->strand = hmm->state[n->state].strand;
		f->state  = n->state;
		f->score  = n->score;
		f->padding = trellis->padding;

		if(n == in) {
			/* MANI now removing initial probabilities from output for the LAST state */
			f->score -= zGetInitProb(trellis->hmm, n->state, trellis->iiso_group);
		}

		while(n != tree->root && n->parent->state == n->state){
			n = n->parent;
		}
		f->genomic_start = n->parent->pos;
		f->cdna_start = n->parent->cdna_pos;

		f->length = zCoorMax((f->genomic_end - f->genomic_start), (f->cdna_end - f->cdna_start));
                if (zIsMatch(hmm->state[n->state].name)) {
                        coor_t iterator, matches = 0;
                        for (iterator = 0; iterator < f->length; iterator++) {
				/* add 1 since it is not a fully closed interval */
                                if (zGetDNAS5(f->genomic, f->genomic_start + 1 + iterator) == zGetDNAS5(f->cdna, f->cdna_start + 1 + iterator)) {
                                        matches++;
                                }
                        }
                        f->percent_identity = (100*matches/f->length);
                } else {
                        f->percent_identity = 0.0;
                }

		if(n->parent != tree->root){
			/*EVAN now removing transition scores from output */
			f->score -= (n->parent->score + zGetTransitionScore(trellis->hmm, n->parent->state, n->state,trellis->tiso_group));
			/* Add back the exit probabilities since it had been fixed in tmap*/
			if (hmm->state[n->state].type == INTERNAL || hmm->state[n->state].type == GINTERNAL) {
				f->score += zScoreDurationGroup(hmm->dmap[hmm->state[n->state].duration], 1, trellis->genomic->gc);
			}
		} else {
			/* MANI now removing initial probabilities from output of the FIRST state*/
			f->score -= zGetInitProb(trellis->hmm, n->state, trellis->iiso_group);

		}
		n = n->parent;
	}
}

/* create all features starting at pos and moving towards the end of the sequence */
static int zSNPViterbiCreateLiveNodes(zViterbi* v,int allele){
	zPairTrellis* trellis = v->pair_trellis;
	coor_t l;
	int j;
	zSFList* sfl = v->sfl;
	zSfeature *f;

	int new_node_count = 0;
	
	for (j = 0; j < trellis->hmm->feature_count; j++) {
		l = v->pos;
		if(trellis->factory[j] != NULL){
			trellis->factory[j]->create5(trellis->factory[j], l, sfl);				
			f = zSFListMoveFirst(sfl);
			while(f != NULL){
				f->strand = '+';
				f->state  = j;
				/*new_node_count += zViterbiAddFeature(v,allele,j,f);*/
				f = zSFListMoveNext(sfl);
			}
			zResetSFList(sfl);
		}
	}		
	
	return new_node_count;
			zSNPViterbiCreateLiveNodes(v,allele); 
                /* Compiler hush: unused function */
}

void zSNPFillAlleleHeader(zViterbi* v,zAlleleDecode* ad){
	int next_pos,i,val;
	zSeqVariant* var;
	
	/* fill in snp values for this allele */
	sprintf(ad->header,"# %d->%d GOOD (%d->%d) ",ad->start,ad->stop,ad->diff_start,ad->diff_stop);
	next_pos = strlen(ad->header);
	/* for each snp in this range */
	for(i = 0; i < ad->snp_count; i++){
		var = v->trellis->dna->seq->variants[i+ad->first_snp];
		/* if this snp is not set in this allele */
		val = ad->vals[i];
		if(val == 0){ 
			continue;
		}
		sprintf(&(ad->header[next_pos]),"SNP %d = %c(%d), ",
						/* +1 adjusts from 0 to 1 based indexing */
						var->real_pos+1,var->values[val],val);
		next_pos += strlen(&(ad->header[next_pos]));
	}
	sprintf(&(ad->header[next_pos]),"score diff = %.2f\n",ad->score_diff);
}

zAFVec* zRunPairViterbi(zPairTrellis* trellis, score_t *path_score){
	zPtrList *list;
	zAlleleDecode *ad;
	zAFVec* afv;
	zAlnFeature* f;
	if(trellis->genomic->seq->var_count > 0){
		zDie("Use zRunSNPPairViterbi not zRunPairViterbi on sequences containing snps.");
	}
	list = zRunSNPPairViterbi(trellis,path_score);
	ad = zPtrListMoveFirst(list);
	afv = zMalloc(sizeof(zAFVec),"zRunViterbi afv");
	zInitAFVec(afv,ad->afl->size);
	f = zAFListMoveFirst(ad->afl);
	while(f != NULL){
		zPushAFVec(afv,f);
		f = zAFListMoveNext(ad->afl);
	}
	qsort(afv->elem, afv->size, sizeof(zAlnFeature), zAFPtrCmp);
	zFreeAlleleDecode(ad);
	zFree(ad);
	zFreePtrList(list);
	zFree(list);
	return afv;
}

/* Functions for initializing and finishing alignments */

static void zExtendState(zViterbi *viterbi, int allele, zPairTrellis *trellis, coor_t gpos, coor_t cpos, int from_state, int to_state) {
	score_t         score;
	zHMM           *hmm           = trellis->hmm;
	zTBTreeNode ****cache         = viterbi->cache[allele];
	zTBTree        *tree          = viterbi->tbtrees[allele];
	int             cache_index   = gpos%viterbi->cache_length;
	zTBTreeNode    *previous_cell = cache[(gpos-zGetGenomicIncrement(hmm,to_state))%viterbi->cache_length][cpos-zGetCDnaIncrement(hmm,to_state)][from_state];
	zTBTreeNode    *tbtn          = zGetTBTreeNode(tree);
	if (previous_cell == NULL) zDie("Cannot extend to state=%d for (%u, %u)", to_state, gpos, cpos);

	/* Get the score */
	score =  zGetScannerScore(trellis,trellis->scanner[trellis->hmm->state[to_state].model],to_state,gpos,cpos);
	score += zGetTransitionScore(hmm,from_state,to_state,trellis->tiso_group);
	score += previous_cell->score;

	/* Set it in the node and cache it */
	cache[cache_index][cpos][to_state] = tbtn;
	tbtn->pos = gpos;
	tbtn->cdna_pos = cpos;
	tbtn->state = to_state;
	tbtn->score = score;
	tbtn->frame_data = 0;
	zTBTreeSetChild(previous_cell,tbtn);
}

/* Start the alignment at position (gpos, cpos) by calculating the initial probabilities
   and making the nodes with that probabilities. Also update the cache for that position 
   Extend the unaligned cDNA till the end of cDNA sequence
*/

static void zStartAlignmentForward(zViterbi* viterbi, int allele, zPairTrellis *trellis, coor_t gpos, coor_t cpos) {
	zHMM           *hmm         = trellis->hmm;
	zTBTreeNode ****cache       = viterbi->cache[allele];
	zTBTree        *tree        = viterbi->tbtrees[allele];
	zTBTreeNode    *tbtn;
	int             state;
	score_t         score;
	coor_t          real_gpos, real_cpos;
	for (state = 0; state < hmm->states; state++) {
		/* Start state */
		score  = zGetFixedInitProb(hmm, state, trellis->iiso_group, trellis->genomic->gc); /* Fixed Initial Probability */
		if(score > MIN_SCORE){
			tbtn           = zGetTBTreeNode(tree);
			tbtn->pos      = gpos;
			tbtn->cdna_pos = cpos;
			tbtn->state    = state;
			tbtn->score    = score;
			zTBTreeSetChild(tree->root,tbtn);
			cache[gpos%viterbi->cache_length][cpos][state] = tbtn;
		} else {
			cache[gpos%viterbi->cache_length][cpos][state] = NULL;
			continue;
		}

		/* First instance of state should not have transition prob in it. *
                 * If you pass it to Viterbi, it will add transition prob to it.  */
		real_gpos = gpos + zGetGenomicIncrement(hmm, state);
		real_cpos = cpos + zGetCDnaIncrement(hmm, state);
		score += zGetScannerScore(trellis, trellis->scanner[trellis->hmm->state[state].model], state, real_gpos, real_cpos); /* Score that pos */
		if(score > MIN_SCORE){
			tbtn           = zGetTBTreeNode(tree);
			tbtn->pos      = real_gpos;
			tbtn->cdna_pos = real_cpos;
			tbtn->state    = state;
			tbtn->score    = score;
			zTBTreeSetChild(tree->root,tbtn);
			cache[real_gpos%viterbi->cache_length][real_cpos][state] = tbtn;
		} else {
			cache[real_gpos%viterbi->cache_length][real_cpos][state] = NULL;
		}
	}

	/* Get the 5' unaligned cDNA state */

	state = zGetFivePrimeCDna(trellis->hmm);
	if (state == -1) {
		zDie("Cannot find 5' unaligned cDNA");
	}
	if(hmm->state[state].type != INTERNAL){
		zDie("GPAIRHMM doesn't support non-INTERNAL states in unaligned regions");							
	}

	/* Extend the unaligned 5' cDNA */

	/* cpos + 1 is done */
	for (real_cpos = cpos + 2; real_cpos < trellis->cdna->length; real_cpos++) {
		zTBTreeNode* previous_cell = cache[gpos%viterbi->cache_length][real_cpos-zGetCDnaIncrement(hmm,state)][state];
		zTBTreeNode* tbtn;

		/* Get the score */
		score =  zGetScannerScore(trellis,trellis->scanner[trellis->hmm->state[state].model],state,gpos,real_cpos);
		if (score == MIN_SCORE) continue;

		score += zGetTransitionScore(hmm,state,state,trellis->tiso_group);

		if (previous_cell == NULL) zDie("Cannot traceback from (%u, %u), state=%d", gpos, real_cpos, state);
		score += previous_cell->score;

		/* Set it in the node and cache it */
		tbtn             = zGetTBTreeNode(tree);
		tbtn->pos        = gpos;
		tbtn->cdna_pos   = real_cpos;
		tbtn->state      = state;
		tbtn->score      = score;
		tbtn->frame_data = 0;
		zTBTreeSetChild(previous_cell,tbtn);
		cache[gpos%viterbi->cache_length][real_cpos][state] = tbtn;
	}
	
}

/* Finish the alignment at position (gpos, cpos) by calculating the initial probabilities
   and adding that probabilities to the cache for that position */

static void zFinishAlignmentForward(zTBTreeNode ***cache, zPairTrellis *trellis, coor_t cpos) {
	int state;
	zHMM *hmm = trellis->hmm;
	for (state = 0; state < hmm->states; state++) {
		score_t score = zGetInitProb(hmm, state, trellis->iiso_group);
		if (cache[cpos][state] != NULL)
			cache[cpos][state]->score += score;
	}
}

/* End Functions for initializing and finishing alignments */

void zRunSNPPairViterbiOnBlock(zViterbi* viterbi, zPairTrellis* trellis, zHSP* block){
	zHMM*         hmm = trellis->hmm;
	zDNA*         genomic = trellis->genomic;
	zDNA*         cdna = trellis->cdna;
	score_t       score,best_score,this_score;
	int           state,best_state,prestate,i;
	zIVec         *jumps;
	zTBTreeNode  *tbtn;
	coor_t        next_snp_pos;
	int           snp_idx;
	int           allele;
	zTBTree*      tree;
	zPtrList*     nodes;
	coor_t        cpos;
	int           cstate;
	coor_t        cdna_pos;
	coor_t        cache_index;
	zTBTreeNode ****cache;

	/* init snp tracking */
	snp_idx = 0;
	if(genomic->seq->var_count == 0){
		next_snp_pos = genomic->length+1;
	}
	else{
		next_snp_pos = genomic->seq->variants[0]->pos;
	}
	
	cpos = -1;
	cstate = -1;

	/* walk forward through sequence */
	for(viterbi->pos = MAX(block->g_start, trellis->blocks->gb_start); viterbi->pos <= MIN(block->g_end, genomic->length - trellis->padding - 1); viterbi->pos++){
		
		cache_index = viterbi->pos%viterbi->cache_length;
		for(allele = 0; allele < viterbi->trellis_count; allele++){
			if(viterbi->active[allele] == 0) continue;
			cache = viterbi->cache[allele];
			zSNPSetSeqForAllele(viterbi,allele); /* No effect */
			/* STEP 1 - create new live nodes */
/*
			zSNPViterbiCreateLiveNodes(viterbi,allele);
*/
			/* STEP 2 - step one pos forward */			
			tree = viterbi->tbtrees[allele];
			nodes = viterbi->live_nodes[allele];

			/* for EVAN's comparison *
			if (viterbi->pos % 100 == 0) {
				fprintf(stderr, "GPos: %u, Tree size: %d\n", viterbi->pos, tree->size);
			}
			*/

			state = zGetFivePrimeGenomic(trellis->hmm);
			zExtendState(viterbi, allele, trellis, viterbi->pos, trellis->padding - 1, state, state);

			for (viterbi->cdna_pos = MAX(block->c_start, trellis->padding); viterbi->cdna_pos <= MIN(block->c_end, cdna->length - trellis->padding - 1); viterbi->cdna_pos++) {
				for(state = 0;state < hmm->states;state++){
					/* handle INTERNAL state transitions */ 
					if (cache[cache_index][viterbi->cdna_pos][state] != NULL) continue;
					if(hmm->state[state].type == INTERNAL){
						zTBTreeNode** previous_array;
						int           previous_cache_index;
						best_state = -1;
						best_score = MIN_SCORE;
						jumps = hmm->jmap[state];
						this_score = zGetScannerScore(trellis,trellis->scanner[trellis->hmm->state[state].model],state,viterbi->pos,viterbi->cdna_pos);
						if (this_score == MIN_SCORE) continue;
					/* I added viterbi->cache_length here because the modulus of a negative number is returned as
					   negative. For example, -1%3 = -1, -4%3 = -1, 1%3 = 1, 4%3 = 1.
					   So I do -1%3 => 2%3 = 2, etc.  Note that viterbi->cache_length is ALWAYS greater than the
					   MAXIMUM genomicincrement. If that changes, this will break */
						previous_cache_index = viterbi->cache_length + cache_index - zGetGenomicIncrement(hmm,state);
						previous_cache_index %= viterbi->cache_length;
						previous_array = cache[previous_cache_index][viterbi->cdna_pos-zGetCDnaIncrement(hmm,state)];
						for(i = 0; i < jumps->size; i++) {
							zTBTreeNode* previous_cell;
							prestate = jumps->elem[i];
							previous_cell = previous_array[prestate];
							if(previous_cell == NULL) continue;
							
							score = this_score + zGetTransitionScore(hmm,prestate,state,trellis->tiso_group);

							if(score != MIN_SCORE){
								score += previous_cell->score;
								if(score > best_score){
									best_score = score;
									best_state = prestate;
								}
							}
						}
						if(best_score != MIN_SCORE){
							/* fill in traceback stuff here */
							zTBTreeNode* previous_cell = previous_array[best_state];
							if (previous_cell == NULL) zDie("Cannot find previous state for (%u, %u), state=%d", viterbi->pos, viterbi->cdna_pos, state);
							tbtn = zGetTBTreeNode(tree);
							cache[cache_index][viterbi->cdna_pos][state] = tbtn;
							tbtn->pos = viterbi->pos;
							tbtn->cdna_pos = viterbi->cdna_pos;
							tbtn->state = state;
							tbtn->score = best_score;
							tbtn->frame_data = 0;
							zTBTreeSetChild(previous_cell,tbtn);
						}
					} else if(hmm->state[state].type == EXPLICIT){
						zTBTreeNode** previous_array;
						int           previous_cache_index;
						int best_length = 0;
						zStrIdx    name;
						zScanner *scanner;
						int gincrement, cincrement;
						zDurationGroup *group;
						zDistribution *d;
						coor_t     gpos    = viterbi->pos, cpos = viterbi->cdna_pos;
						best_state = -1;
						best_score = MIN_SCORE;
					/* I added viterbi->cache_length here because the modulus of a negative number is returned as
					   negative. For example, -1%3 = -1, -4%3 = -1, 1%3 = 1, 4%3 = 1.
					   So I do -1%3 => 2%3 = 2, etc.  Note that viterbi->cache_length is ALWAYS greater than the
					   MAXIMUM genomicincrement. If that changes, this will break */
						previous_cache_index = viterbi->cache_length + cache_index - zGetGenomicIncrement(hmm,state);
						previous_cache_index %= viterbi->cache_length;
						previous_array = cache[previous_cache_index][viterbi->cdna_pos-zGetCDnaIncrement(hmm,state)];

						gincrement = zGetGenomicIncrement(trellis->hmm, state);
						cincrement = zGetCDnaIncrement(trellis->hmm, state);
						group      = trellis->hmm->dmap[trellis->hmm->state[state].duration];
						d          = &group->duration[0].distribution[0];
						name       = trellis->hmm->state[state].model;
						scanner    = trellis->scanner[name];
									

						jumps = hmm->jmap[state];
						for(i = 0; i < jumps->size; i++) {
							int        from_state = jumps->elem[i];
							score_t    tscore, scan_score, total_score;
							coor_t     length;

							coor_t gmin, gmax, cmin, cmax;
							coor_t steps;

							/* Set min and max start coordinates for the range */
							steps = d->end - d->start + 1;
							if ((int)(gpos - steps*gincrement) < (int)trellis->blocks->gb_start - 1) {
								steps = (gpos - trellis->blocks->gb_start + 1)/gincrement;
							}
							if ((int)(cpos - steps*cincrement) < PADDING - 1) {
								steps = (cpos - PADDING + 1)/cincrement;
							}
							
							scan_score = 0.;
							best_score  = MIN_SCORE;
							tscore  = zGetTransitionScore(trellis->hmm, from_state, state, trellis->tiso_group);

							for (length = 1, gmax = gpos, cmax = cpos, gmin = gmax - gincrement, cmin = cmax - cincrement;
							     length <= steps;
							     length++, gmax = gmin, cmax = cmin, gmin = gmax - gincrement, cmin = cmax - cincrement) {

								coor_t i, j;
								zTBTreeNode* tbtn = cache[gmin%viterbi->cache_length][cmin][from_state];
								if (tbtn == NULL) continue;

								total_score = tscore + zScoreDistribution(d, length) + tbtn->score;
									
								for (i = gmin+1, j = cmin+1; i <= gmax; i+=gincrement, j+=cincrement) {
									scan_score += zGetScannerScore(trellis, scanner, state, i, j);
								}
								total_score += scan_score;

								if (total_score >= best_score) {
									best_score  = total_score;
									best_length = length;
									best_state  = from_state;
								}
							  
							}
						}
						if(best_score != MIN_SCORE){
							/* fill in traceback stuff here */
							zTBTreeNode* previous_cell = cache[(gpos-best_length*gincrement)%viterbi->cache_length][cpos-best_length*cincrement][best_state];
							if (previous_cell == NULL) zDie("Cannot find previous state for (%u, %u), state=%d", viterbi->pos, viterbi->cdna_pos, state);
							tbtn = zGetTBTreeNode(tree);
							cache[cache_index][viterbi->cdna_pos][state] = tbtn;
							tbtn->pos = viterbi->pos;
							tbtn->cdna_pos = viterbi->cdna_pos;
							tbtn->state = state;
							tbtn->score = best_score;
							tbtn->frame_data = 0;
							zTBTreeSetChild(previous_cell,tbtn);
						}
					} else if(hmm->state[state].type == EXTERNAL){
					/* EXTERNAL states dealt with in when incorporating live nodes below */
					}
					else{
						zDie("zRunSNPPairViterbi only supports INTERNAL and EXTERNAL states");
					}
				}
			}

			/* cleanup unused/unneeded old_cells, and reinit next cache_index */
			for (cdna_pos = 0; cdna_pos < cdna->length; cdna_pos++) {
				for (state = 0;state < hmm->states;state++){
					tbtn = cache[(cache_index+1)%viterbi->cache_length][cdna_pos][state];
					if(tbtn == NULL) continue;
					if(tbtn->children == 0){
						zReleaseDeadTBTreeNode(tree,tbtn);
					}
					else if((tbtn->children == 1) &&
									(tbtn->child->state == tbtn->state)){
						zReleaseRedundantTBTreeNode(tree,tbtn);
					}
					cache[(cache_index+1)%viterbi->cache_length][cdna_pos][state] = NULL;
				}
			}
			
			/* STEP 5 - remove/incorporate old live nodes */
			zSNPUnSetSeqForAllele(viterbi,allele);
		}
	}
}

zPtrList* zRunSNPPairViterbi(zPairTrellis* trellis, score_t *path_score){
	zViterbi*     viterbi;
	zHMM*         hmm = trellis->hmm;
	zDNA*         genomic = trellis->genomic;
	zAlleleDecode* ad;
	score_t       best_score;
	int           state,best_state,i;
	coor_t        next_snp_pos;
	int           snp_idx;
	int           allele;
	zTBTree*      tree;
	zPtrList*     nodes;
	coor_t        cpos;
	int           cstate;
	zPtrList*     ret_list;
	zAlnFeature*    f;
	zTBTreeNode **cells;

	/* Prepare Viterbi Vars */
	viterbi = zMalloc(sizeof(zViterbi),"zRunViterbi viterbi");
	zInitPairViterbi(viterbi,trellis);

	/* Init allele 0 */
	allele = 0;
	viterbi->first_snp[0] = -1;
	viterbi->snp_count[0] = 0;
	viterbi->snp_int_vals[0] = 0;
	tree = viterbi->tbtrees[0];
	nodes = viterbi->live_nodes[0];
	viterbi->trellis_count = 1;
	viterbi->tbt_use_count[0] = 1;
	viterbi->start_pos[0] = trellis->padding-1;
	viterbi->active[0] = 1;

	for (state = 0; state < hmm->states; state++) {
		if(hmm->state[state].type == GINTERNAL){
			zDie("PairHMM doesn't support GINTERNAL");							
		}
	}

	/* Init allele 0 cells */
	/* Sending cache[0] for gpos=trellis->gb_start-1 and viterbi starts at gpos=trellis->padding->gb_start. Remember this when you handle cache inside zStartAlignmentForward() */
	zStartAlignmentForward(viterbi, allele, trellis, trellis->blocks->gb_start-1, trellis->padding-1);

	/* init snp tracking */
	snp_idx = 0;
	if(genomic->seq->var_count == 0){
		next_snp_pos = genomic->length+1;
	}
	else{
		next_snp_pos = genomic->seq->variants[0]->pos;
	}
	
	cpos = -1;
	cstate = -1;

	/* for EVAN's comparison *
	{
		long area = 0;
		for (i = 0; i < trellis->mem_blocks->hsps; i++) {
			zHSP* block = &trellis->mem_blocks->hsp[i];
			area += (MIN(block->g_end, genomic->length - trellis->padding - 1) - MAX(block->g_start, trellis->blocks->gb_start) + 1)*(MIN(block->c_end , trellis->cdna->length - trellis->padding - 1) - MAX(block->c_start, trellis->padding) + 1);
		}
		fprintf(stderr, "Viterbi Non-Stepping-stone cells: %16ld\n", (long)(trellis->blocks->gb_end - trellis->blocks->gb_start + 1)*(trellis->cdna->length - 2*trellis->padding)*trellis->hmm->states);
		fprintf(stderr, "Viterbi     Stepping-stone cells: %16ld\n", area*trellis->hmm->states);
	}
	*/

	/* walk forward through sequence */
	for (i = 0; i < trellis->mem_blocks->hsps; i++) {
		zHSP* block = &trellis->mem_blocks->hsp[i];
		zRunSNPPairViterbiOnBlock(viterbi, trellis, block);
	}
	
	/*zPrintNodeGraph();EVAN*/
	zFinishAlignmentForward(viterbi->cache[0][(viterbi->pos-1)%viterbi->cache_length], trellis, viterbi->cdna_pos - 1); /* viterbi->cdna_pos is 1 more than max */

	/* Viterbi trace back */
	zTrace2("traceback\n");
	*path_score = 0;

	cells = viterbi->cache[0][(viterbi->pos-1)%viterbi->cache_length][viterbi->cdna_pos - 1]; /* viterbi->cdna_pos is 1 more than max */
	best_state = -1;
	best_score = MIN_SCORE;
	for(state = 0;state < hmm->states;state++){
		if(cells[state] != NULL && 
			cells[state]->score > best_score){
			best_state = state;
			best_score = cells[state]->score;
		}
	}
	if(best_state == -1){
		zDie("No live states at end of trellis.  This should never happen");
	}
	*path_score = best_score;
	
	viterbi->trellis_count = 1;
	
	ad = zMalloc(sizeof(zAlleleDecode),"zRunSNPViterbi ad");
	zInitAlleleDecode(ad,25,1);
	ad->vals[0] = 0;
	
	sprintf(ad->header,"## Main SNP Sequence\n"); 
	
	zTBTree2AFL(ad->afl,viterbi->tbtrees[0],cells[best_state],hmm,trellis);
	f = zAFListMoveLast(ad->afl);
	zPtrListAddFirst(viterbi->traceback,ad);
	zSNPCleanUpDecodes(viterbi);
	ret_list = viterbi->traceback;
	viterbi->traceback = NULL;
	zFreePairViterbi(viterbi);
	zFree(viterbi);
	return ret_list;
}
