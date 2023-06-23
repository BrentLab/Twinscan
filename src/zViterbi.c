#include <time.h>
#include "zTBTree.h"
#include "zTrellis.h"
#include "zViterbi.h"

/**********************************************************************\
  zViterbi

  This holds the data needed for the memory optimized viterbi run

\**********************************************************************/

/* zViterbi structure holds the trellis, traceback tree and tbnodes */ 

static const coor_t  VITERBI_CACHE_LENGTH         = 2;
static const int     MAX_CONCURRENT_SEQUENCE_VARIANTS = 16384; /*2^14*/
static const int     MAX_CONCURRENT_SNPS              = 6;
static const int     MAX_ACTIVE_SNPS                  = 6;
static const float   HAPLOTYPE_FREQUENCY_THRESHOLD    = 0.;

extern coor_t any_new_cpoint;

static void zSNPTraceAllele(zViterbi* v, int a, bool cptrace, score_t base_score, zSFList* base_list,char* text);

/*EVAN
void zCheckSNPViterbi(zViterbi *v){	
	int i;
	for(i = 0; i < v->max_trellis_count; i++){
		if(i < v->trellis_count){
			if(v->active[i] == 0 &&
			   v->tbt_use_count[i] == 0){
				v->temp[i] = 0;
			}
			else{
				v->temp[i] = 1;
			}
		}
		else{
			v->temp[i] = 0;
			if(v->tbt_use_count[i] != 0){
				zDie("zCSV 1 %d",i);
			}
		}
	}
	for(i = 0; i < v->max_trellis_count; i++){
		v->temp[v->next_tbti[i]]++;
		if(v->next_tbtn[i] != NULL && 
		   v->next_tbtn[i]->tree != v->tbtrees[v->next_tbti[i]]){
			zDie("zCSV 2 %d",i);
		}
	}
	for(i = 1; i < v->max_trellis_count; i++){
		if(v->temp[i] != v->tbt_use_count[i]){
			zDie("zCSV 3 %d",i);
		}
	}
	fprintf(stderr,"zCheckSNPViterbi PASSED\n");
}
*/

/* EVAN live nodes would be better also stored in a hash table, then we woudln't 
	 have to search for them */
/* get/create a live node for a given state/pos */
zTBTreeNode* zFindViterbiLiveNode(zViterbi *v,int allele, coor_t pos,int state){
	zPtrList* l = v->live_nodes[allele];
	zTBTree* t = v->tbtrees[allele];
	zTBTreeNode* n;

	/* search through until we find the node or find where it should go */
	n = zPtrListMoveFirst(l);
	while(n != NULL && n->pos < pos){
		n = zPtrListMoveNext(l);
	}
	if(n != NULL && n->pos == pos){
		while(n != NULL && n->pos == pos && n->state < state){
			n = zPtrListMoveNext(l);
		}
		if(n != NULL && n->pos == pos && n->state == state){
			return n;
		}
	}
	/* didn't find it so insert a new one */
	n = zGetTBTreeNode(t);
	n->pos = pos;
	n->state = state;
	zPtrListAddPrev(l,n);
	return n;
}

void zRemoveViterbiLiveNode(zViterbi *v,int t, zTBTreeNode* rem){
	zPtrList* l = v->live_nodes[t];
	zTBTreeNode* n;

	/* search through until we find the node or find where it should go */
	n = zPtrListMoveFirst(l);
	while(n != NULL && n->pos <= rem->pos){
		if(n == rem){
			zPtrListRemoveCurrent(l);
			return;
		}
		n = zPtrListMoveNext(l);
	}
}

/* inserts node n into the sorted live_nodes list */
static void zViterbiLiveNodeInsert(zPtrList* live_nodes,zTBTreeNode* n){
	zTBTreeNode* ln = zPtrListMoveFirst(live_nodes);
	while(ln != NULL){
		if((ln->pos > n->pos) ||
			 (ln->pos == n->pos && ln->state > n->state)){
			zPtrListAddPrev(live_nodes,n);
			return;
		}
		ln = zPtrListMoveNext(live_nodes);
	}
	zPtrListAddLast(live_nodes,n);
}

/* walks the tbtree copying nodes into the cell array and live_node list */
static void zViterbiFillNodesHelper(coor_t pos,zTBTreeNode** cells,
									zPtrList* live_nodes, zTBTreeNode* n){
	zTBTreeNode* c = n->child;
	if(n->pos == pos-1){
		cells[n->state] = n;
	}
	else if(n->pos >= pos){
		zViterbiLiveNodeInsert(live_nodes,n);
	}
	while(c != NULL){
		zViterbiFillNodesHelper(pos,cells,live_nodes,c);
		c = c->rsib;
	}
}

/* fills the cell array and live_node list with the nodes from 
	 the tbtree */
static void zViterbiFillNodes(zViterbi* v,int allele){
	zTBTree* tree = v->tbtrees[allele];
	zTBTreeNode** cells = v->cache[allele][0][v->cdna_pos];
	zPtrList* live_nodes = v->live_nodes[allele];
	zViterbiFillNodesHelper(v->pos,cells,live_nodes,tree->cpoint);
}

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

int zBackwardViterbiIncorporateLiveNodes(zViterbi* v, int allele){
	int removed_node_count = 0;
	zPtrList* l = v->live_nodes[allele];
	zTBTreeNode** cells = v->cache[allele][0][v->cdna_pos];
	zTBTreeNode* n = NULL;

	n = zPtrListMoveLast(l);
	if(n != NULL && n->pos > v->pos){
		zDie("WTF");
	}
	/* until the first live node comes after pos or we reach the end of the list */
	while(n != NULL && n->pos == v->pos){
		/* incorporate the live node */
		cells[n->state] = n;
		/* move to the next node */
		n = zPtrListMovePrev(l);
		/* remove the incorporated live node */
		zPtrListRemoveNext(l);
		removed_node_count++;
	}
	return removed_node_count;
}

/* test for same live_node lists after copy
void zCheckSameLists(zPtrList* l1, zPtrList* l2){
	zTBTreeNode* n1;
	zTBTreeNode* n2;
	int pos;
	n1 = zPtrListMoveFirst(l1);
	n2 = zPtrListMoveFirst(l2);
	pos = 0;
	while(n1 != NULL && n2 != NULL){
		if(n1->pos != n2->pos || n1->state != n2->state){
			printf("pos %d (%d %d != %d %d)\n",pos,n1->pos,n1->state,n2->pos,n2->state);
			pos++;
		}
		n1 = zPtrListMoveNext(l1);
		n2 = zPtrListMoveNext(l2);
	}
	if((n1 != NULL)||(n2 != NULL)){
		printf("last pos %p %p\n",(void*)n1,(void*)n2);
	}
}
*/

/* copy allele a1 to a2 including live_nodes, tbtree, and cells */
void zSNPCopyAllele(zViterbi* v, int a1, int a2){
	int i;
	if(a2 >= v->max_trellis_count){
		zDie("Cannot split allele %d to %d since max allele is %d\n",
				 a1,a2,v->max_trellis_count-1); 
	}

	/*fprintf(stderr," SPLIT %d->%d\n",a1,a2);*/

	v->first_snp[a2] = v->first_snp[a1];
	v->snp_count[a2] = v->snp_count[a1];
	v->active_snp_count[a2] = v->active_snp_count[a1];
	v->start_pos[a2] = v->start_pos[a1];
	v->next_tbtn[a2] = NULL;
	v->next_tbti[a2] = 0;
	v->tbt_use_count[a2] = 1;

	for(i = 0;i < v->snp_count[a2]; i++){
		v->snp_vals[a2][i] = v->snp_vals[a1][i];
	}
	v->snp_int_vals[a2] = v->snp_int_vals[a1];
	/* need to copy using next_tbtn info */
	zCopyTBTreeFromCPoint(v->tbtrees[a1],v->tbtrees[a2],v->start_pos[a1]);
	if(v->next_tbtn[a1] != NULL){	
		zTBTreeNode* node = v->next_tbtn[a1];
		int next = v->next_tbti[a1];
		while(node != NULL){
			zAppendTBTreeFromNode(v->tbtrees[next],node,v->tbtrees[a2],v->tbtrees[a2]->cpoint);
			node = v->next_tbtn[next];
			next = v->next_tbti[next];
		}
	}
	zViterbiFillNodes(v,a2);
}

void zSNPSwapAlleles(zViterbi* v, int a1,int a2){
	zTBTreeNode** cells = v->cache[a1][0][v->cdna_pos];
	zTBTree* tree = v->tbtrees[a1];
	zPtrList* list = v->live_nodes[a1];
	short* snp_vals = v->snp_vals[a1];
	int snp_int_val = v->snp_int_vals[a1];
	int first_snp = v->first_snp[a1];
	int snp_count = v->snp_count[a1];
	zTBTreeNode* next_tbtn = v->next_tbtn[a1];
	int next_tbti = v->next_tbti[a1];
	int active = v->active[a1];
	int use_count = v->tbt_use_count[a1];
	int active_snp_count = v->active_snp_count[a1];

	v->cache[a1][0][v->cdna_pos] = v->cache[a2][0][v->cdna_pos];
	v->tbtrees[a1] = v->tbtrees[a2];
	v->live_nodes[a1] = v->live_nodes[a2];
	v->snp_vals[a1] = v->snp_vals[a2];
	v->snp_int_vals[a1] = v->snp_int_vals[a2];
	v->first_snp[a1] = v->first_snp[a2];
	v->snp_count[a1] = v->snp_count[a2];
	v->next_tbtn[a1] = v->next_tbtn[a2];
	v->next_tbti[a1] = v->next_tbti[a2];
	v->active[a1] = v->active[a2];
	v->tbt_use_count[a1] = v->tbt_use_count[a2];
	v->active_snp_count[a1] = v->active_snp_count[a2];

	v->cache[a2][0][v->cdna_pos] = cells;
	v->tbtrees[a2] = tree;
	v->live_nodes[a2] = list;
	v->snp_vals[a2] = snp_vals;
	v->snp_int_vals[a2] = snp_int_val;
	v->first_snp[a2] = first_snp;
	v->snp_count[a2] = snp_count;
	v->next_tbtn[a2] = next_tbtn;
	v->next_tbti[a2] = next_tbti;
	v->active[a2] = active;
	v->tbt_use_count[a2] = use_count;
	v->active_snp_count[a2] = active_snp_count;
}

void zSNPPrintAlleles(zViterbi* v){
	int a,i,first,last;

	first = v->first_snp[0];
	last = v->first_snp[v->trellis_count-1]+
		v->snp_count[v->trellis_count-1]-1;
	if(first == -1){
		fprintf(stderr,"NO SNPS\n");		
	}
	else{
		fprintf(stderr,"%d->%d\n",first,last);
	}
	for(a = 0; a < v->trellis_count; a++){
		fprintf(stderr,"  %d (%d):",a,v->first_snp[a]);
		if(v->snp_count[a] == 0){
			fprintf(stderr,"\tnone\n");
			continue;
		}
		for(i = first; i < v->first_snp[a]; i++){
			fprintf(stderr,"\t-");
		}
		for(i = 0; i < v->snp_count[a]; i++){
			fprintf(stderr,"\t%d",v->snp_vals[a][i]);
		}
		for(i = v->first_snp[a]+v->snp_count[a]; i <= last; i++){
			fprintf(stderr,"\t-");
		}
		fprintf(stderr,"\n");
	}
}

void zSNPReleaseAllele(zViterbi* v, int a){
	int i;
	v->snp_count[a] = 0;
	v->active_snp_count[a] = 0;
	v->first_snp[a] = -1;
	for(i = 0; i < v->snp_count[a]; i++){
		v->snp_vals[a][i] = 0;
	}
	for(i = 0; i < v->hmm->states; i++){
		v->cache[a][0][v->cdna_pos][i] = NULL;
	}
	v->next_tbtn[a] = NULL;
	v->next_tbti[a] = 0;
	v->active[a] = 0;
	zResetPtrList(v->live_nodes[a]);
	zResetTBTree(v->tbtrees[a]);			
}

void zPinReleaseTrellis(zViterbi* v, int a){
	int i;
	for(i = 0; i < v->hmm->states; i++){
		v->cache[a][0][v->cdna_pos][i] = NULL;
	}
	zResetPtrList(v->live_nodes[a]);
	zResetTBTree(v->tbtrees[a]);			
}

/*
static bool zSNPDiffAllele(zViterbi* v,int a1,int a2){
	int i;
	if(v->first_snp[a1] != v->first_snp[a2] ||
		 v->snp_count[a1] != v->snp_count[a2]){
		return true;
	}
	for(i = 0; i < v->snp_count[a1]; i++){
		if(v->snp_vals[a1][i] != v->snp_vals[a2][i]){
			return true;
		}
	}
	return false;
}
*/

void zInitAlleleDecode(zAlleleDecode* ad, size_t header_size, size_t val_size){
	ad->sfv = zMalloc(sizeof(zSFVec),"zInitAlleleDecode ad->sfv");
	zInitSFVec(ad->sfv,1);
	ad->sfl = zMalloc(sizeof(zSFList),"zInitAlleleDecode ad->sfl");
	zInitSFList(ad->sfl);
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
	if(ad->header != NULL){
		zFree(ad->header);
		ad->header = NULL;
	}
	if(ad->vals != NULL){
		zFree(ad->vals);
		ad->vals = NULL;
	}
}

void zSnpViterbiMergeSFLists(zSFList* base, zSFList* ext){
	zSfeature *f;
	zSfeature *bf = zSFListMoveFirst(base);
	zSfeature *ef = zSFListMoveLast(ext);
	while(ef->end < bf->end){
		ef = zSFListMovePrev(ext);
	}
	if(ef->state == bf->state){
		bf->end = ef->end;
		ef = zSFListMovePrev(ext);
	}
	else{
		ef->start = bf->end+1;
	}
	while(ef != NULL){
		f = zSFListPrepend(base);
		zCopySfeature(ef,f);
		ef = zSFListMovePrev(ext);
	}
}

/* this subroutine traces and removes the specified allele and any other allele whose tbtree 
   links to it */
void zSNPCollapseAllele(zViterbi* v, int allele){
	char* text;

	text = zMalloc(sizeof(char)*10,"");
	sprintf(text,"%d",allele);
	zSNPTraceAllele(v,allele,true,v->tbtrees[0]->cpoint->score,NULL,text);
	zFree(text);
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

/* EVAN Heuristic cpoint search tests 

zPtrList leaves;
bool leaves_ready = false;
int node_graph[11];

void zTBTreeFillDescendantScores(zTBTreeNode* n){
	zTBTreeNode* c = n->child;
	if(c == NULL){
		n->dec_score = n->score;
		return;
	}
	n->dec_score = 0;
	while(c != NULL){
		zTBTreeFillDescendantScores(c);
		n->dec_score += c->dec_score;
		c = c->rsib;
	}
}

void zTBTreeClearSimTrims(zTBTreeNode* n){
	zTBTreeNode* c = n->child;
	n->sim_trim = 0;
	while(c != NULL){
		zTBTreeClearSimTrims(c);
		c = c->rsib;
	}
}

void zTBTreeSimKeepXPercentBest(zTBTree* tree, int n, zPtrList* leaves);

void zTBTreeSimKeepN(zTBTree* tree, int n, zPtrList* leaves);

void zTBTreeSimDropXBelowAve(zTBTree* tree, float x);
void zTBTreeSimDropXBelowMax(zTBTree* tree, float x);

void  zTBTreeSortLeafListByScore(zPtrList* leaves){
	int i;
	score_t score;
	int last = 0;
	zPtrListPos base;
	zTBTreeNode* n;
	zPtrListMoveFirst(leaves);
	zPtrListMoveNext(leaves);
	base = zPtrListGetPos(leaves);
	n = zPtrListMovePrev(leaves);
	i = 0;
	while(last != 2){
		/ move this node back to where is belongs /
		zPtrListPos this_pos = zPtrListGetPos(leaves);
		zTBTreeNode* n2 = zPtrListMovePrev(leaves);

		while(n2 != NULL && n->score > n2->score){
			n2 = zPtrListMovePrev(leaves);
		}
		zPtrListMovePosToNext(leaves,this_pos);
		zPtrListSetPos(leaves,base);
		n = zPtrListMoveNext(leaves);
		if(n == NULL){
			if(last == 0){
				last = 1;
			}
			else{
				last = 2;
			}
		}
		base = zPtrListGetPos(leaves);
		n = zPtrListMovePrev(leaves);
	}

	n = zPtrListMoveFirst(leaves);
	score = n->score;
	while(n != NULL){
		if(n->score > score){
			zDie("Sorting didn't work\n");
		}
		score = n->score;
		n = zPtrListMoveNext(leaves);
	}
}

void zTBTreeCheckHeuristicCpoints(zViterbi* v,zTBTree* tree,zTBTreeNode** cells,zPtrList* nodes){
	zHMM* hmm = v->hmm;
	zTBTreeNode* n;
	score_t total,ave,stddev,max;
	int i,count;
	
	if(!leaves_ready){
		leaves_ready = true;
		zInitPtrList(&leaves);
		for(i = 0; i <= 10; i++){
			node_graph[i] = 0;
		}
	}
	zResetPtrList(&leaves);
	for(i = 0; i < hmm->states;i++){
		if(cells[i] != NULL){
			zPtrListAddLast(&leaves,cells[i]);
		}
	}
	n = zPtrListMoveFirst(nodes);
	while(n != NULL){
		zPtrListAddLast(&leaves,n);
		n = zPtrListMoveNext(nodes);
	}

	/ calculate leaf average, stddev, max scores /
	total = 0;
	max = MIN_SCORE;
	n = zPtrListMoveFirst(&leaves);
	count = 0;
	while(n != NULL){
		if(n->score > max){
			max = n->score;
		}
		total += n->score;
		count++;
		n = zPtrListMoveNext(&leaves);
	}
	ave = total / count;
	total = 0;
	n = zPtrListMoveFirst(&leaves);
	while(n != NULL){
		total += (n->score - ave)*(n->score - ave);
		n = zPtrListMoveNext(&leaves);
	}
	stddev = sqrt(total / count);

	n = zPtrListMoveFirst(&leaves);
	while(n != NULL){
		int diff = (n->score-ave)/(stddev/2);
		diff += 5;
		if(diff < 0){
			diff = 0;
		}
		if(diff > 10){
			diff = 10;
		}
		node_graph[diff]++;
		n->ng[diff] = 1;
		n = zPtrListMoveNext(&leaves);
	}

	tree = NULL;
}

void zPrintNodeGraph(){
	int i;
	fprintf(stderr,"LEAF NODE SCORE DIST\n");
	for(i = 0; i <= 10; i++){
		int diff = i-5;
		fprintf(stderr,"%d\t%d\n",diff,node_graph[i]);
	}
}
*/

/**********************************************************************\
  zViterbi / zTrellis functions
\**********************************************************************/

static void zInitViterbi(zViterbi* v,zTrellis* trellis){
	int i,j;
	coor_t k,m;
	zHMM*  hmm = trellis->hmm;
	zSeqVariant** vars = trellis->dna->seq->variants;
	int var_count = trellis->dna->seq->var_count;

	v->trellis = trellis;
	v->hmm     = trellis->hmm;

	v->max_trellis_count = MAX_CONCURRENT_SEQUENCE_VARIANTS;
	v->max_snp_count = MAX_CONCURRENT_SNPS;
	v->max_active_snp_count = MAX_ACTIVE_SNPS;
	v->hap_threshold = HAPLOTYPE_FREQUENCY_THRESHOLD;

	if((unsigned int)v->max_snp_count > sizeof(int)*8 + 1){
		zDie("zViterbi max_snp_count %d > %d which is max for this architecture\n",
			 v->max_snp_count,sizeof(int)*8 + 1);
	}

	v->trellis_count = 0;
	v->pos = 0;
	v->cpos = 0;
	v->cdna_pos = 0;
	v->cdna_min = 0;
	v->cdna_max = 0;

	for(i = 0; i < hmm->states; i++){
		if(hmm->state[i].type == EXPLICIT){
			zDie("Explicit state currently must use the zRunViterbiAndForward fucntion for decoding"); 
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
		v->live_nodes[i] = zMalloc(sizeof(zPtrList),
								   "zInitViterbi lide_nodes[i]");
		zInitPtrList(v->live_nodes[i]);
		v->cache[i] = zMalloc(sizeof(zTBTreeNode***)*VITERBI_CACHE_LENGTH, "zInitViterbi cache[i]");
		for (k = 0; k < VITERBI_CACHE_LENGTH; k++) {
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
	for(i = 0;i < var_count; i++){
		int count = 0;
		int k;
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
		for(k = 0;k < v->trellis->dna->seq->hap_count; k++){
			zSeqHaplotype* hap;			
			hap = &v->trellis->dna->seq->haps[k];
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

	v->sfl = zMalloc(sizeof(zSFList),"zInitViterbi sfl");
	zInitSFList(v->sfl);	
	v->traceback = zMalloc(sizeof(zPtrList),
						   "zInitViterbi traceback");
	zInitPtrList(v->traceback);

	/* don't use PIN stuff */
	v->next = NULL;
	v->dead_count = -1;
	v->active_count = -1;
}

static void zInitViterbiForPin(zViterbi* v,zTrellis* trellis,coor_t start_pos){
	int            i,j,state,cur_trellis;
	coor_t         k,m,tmp_start;
	zSFList        tmp_sfl;
	zIVec*         fmap;
	zSfeature*     f;
	zPhase_t       left_phase,right_phase;
	zTBTreeNode*   live_node;
	zHMM*          hmm = trellis->hmm;

	v->active_count = 0;
	v->dead_count = 0;

	v->trellis = trellis;
	v->hmm     = trellis->hmm;

	for(i = 0; i < hmm->states; i++){
		if(hmm->state[i].type == EXPLICIT){
			zDie("Explicit state currently must use the zRunViterbiAndForward fucntion for decoding"); 
		}
	}

	if(start_pos > v->trellis->dna->length-PADDING){
		zDie("zInitViterbiForPin: requested start pos %d is beyond end of sequence",
				 start_pos-PADDING+1);
	}
	
	/* Collect long distance links crossing start_pos */
	v->sfl = zMalloc(sizeof(zSFList),"zInitViterbiForPin sfl");
	zInitSFList(v->sfl);	
	zInitSFList(&tmp_sfl);	
	for (i = 0; i < v->trellis->hmm->feature_count; i++) {
		j = start_pos;
		if(v->trellis->factory[i] != NULL){
			v->trellis->factory[i]->createi(v->trellis->factory[i], j, &tmp_sfl);
			fmap = v->trellis->fmap3[i];
			f = zSFListMoveFirst(&tmp_sfl);
			while(f != NULL){
				f->strand = '+';
				right_phase = zSFEndPhase(v->trellis->dna, f);
				left_phase = zSFBeginPhase(v->trellis->dna, f);
				for(j = 0;j < fmap->size;j++){
					state = fmap->elem[j];
					/* make sure extenral state (state) is compatible with feature */
					if(right_phase != v->trellis->hmm->state[state].phase) continue;					
					f->state = state;
					zSFListInsert(v->sfl,f);
				}
				f = zSFListMoveNext(&tmp_sfl);
			}			
			zResetSFList(&tmp_sfl);
		}
		j = v->trellis->dna->length - 1 - start_pos;
		if(v->trellis->rfactory[i] != NULL){
			v->trellis->rfactory[i]->createi(v->trellis->rfactory[i], j, &tmp_sfl);
			fmap = v->trellis->rfmap3[i];
			f = zSFListMoveFirst(&tmp_sfl);
			while(f != NULL){
				tmp_start = f->start;
				f->start  = v->trellis->dna->length - f->end - 1;
				f->end    = v->trellis->dna->length - tmp_start - 1;
				f->strand = '-';
				right_phase = zSFBeginPhase(v->trellis->rdna, f);
				left_phase = zSFEndPhase(v->trellis->rdna, f);
				for(j = 0;j < fmap->size;j++){
					state = fmap->elem[j];
					/* make sure extenral state (state) is compatible with feature */
					if(right_phase != v->trellis->hmm->state[state].phase) continue;					
					f->state = state;
					zSFListInsert(v->sfl,f);
				}
				f = zSFListMoveNext(&tmp_sfl);
			}			
			zResetSFList(&tmp_sfl);
		}		
	}	
	v->max_trellis_count = v->sfl->size;	
	
	for (i = 0; i < hmm->states; i++) {
		if(hmm->state[i].type == INTERNAL || hmm->state[i].type == GINTERNAL){
			/* these will need their own trellis */
			v->max_trellis_count++;
		}
		else if(hmm->state[i].type == EXTERNAL){
			/* handled by lond distance links */
		}
		else{
			zDie("Pin iscan can only handle INTERNAL and EXTERNAL states\n");
		}
	}
	
	/* allocate and initialize memeroy */
	v->trellis_count = v->max_trellis_count;
	v->pos = 0;
	v->cpos = 0;
	v->cdna_pos = 0;
	v->cdna_min = 0;
	v->cdna_max = 0;
	
	v->tbtrees = zMalloc(sizeof(zTBTree*)*v->max_trellis_count, "zInitViterbiForPin tbtree");
	v->live_nodes = zMalloc(sizeof(zPtrList*)*v->max_trellis_count, "zInitViterbiForPin live_nodes");
	v->cache = zMalloc(sizeof(zTBTreeNode**)*v->max_trellis_count, "zInitViterbiForPin cache");  
	v->active = zMalloc(sizeof(short)*v->max_trellis_count, "zInitViterbiForPin active");
	v->start_pos = zMalloc(sizeof(coor_t)*v->max_trellis_count, "zInitViterbiForPin active");
	v->next = zMalloc(sizeof(int)*v->max_trellis_count, "zInitViterbiForPin active");

	for(i = 0;i < v->max_trellis_count; i++){
		v->tbtrees[i] = zMalloc(sizeof(zTBTree), "zInitViterbi tbtrees[i]");
		zInitTBTree(v->tbtrees[i]);
		v->tbtrees[i]->root->pos = 0;
		v->live_nodes[i] = zMalloc(sizeof(zPtrList), "zInitViterbi lide_nodes[i]");
		zInitPtrList(v->live_nodes[i]);
		v->cache[i] = zMalloc(sizeof(zTBTreeNode***)*VITERBI_CACHE_LENGTH, "zInitViterbi cache[i]");
		for (k = 0; k < VITERBI_CACHE_LENGTH; k++) {
			v->cache[i][k] = zMalloc(sizeof(zTBTreeNode**)*(v->cdna_max - v->cdna_min + 1), "zInitViterbi cache[i][k]"); 
			for (m = v->cdna_min; m <= v->cdna_max; m++) {
				v->cache[i][k][m] = zMalloc(sizeof(zTBTreeNode*)*hmm->states, "zInitViterbi cache[i][k][m]"); 
				for(j = 0; j < hmm->states; j++){
					v->cache[i][k][m][j] = NULL;
				}
			}
		}
	}
	
	cur_trellis = 0;
	/* add internal nodes here */
	for (i = 0; i < hmm->states; i++) {
		if(hmm->state[i].type != INTERNAL && hmm->state[i].type != GINTERNAL){
			continue;
		}
		
		for(j = 0; j < hmm->states;j++){
			v->cache[cur_trellis][0][v->cdna_pos][j] = NULL;
		}
		live_node = zGetTBTreeNode(v->tbtrees[cur_trellis]);
		live_node->pos = start_pos;
		live_node->state = i;
		live_node->score = 0;
		zTBTreeSetChild(v->tbtrees[cur_trellis]->root,live_node);
		v->cache[cur_trellis][0][v->cdna_pos][i] = live_node;
		v->active[cur_trellis] = 1;
		v->active_count++;
		v->start_pos[cur_trellis] = start_pos+1;
		v->next[cur_trellis] = cur_trellis+1;
		cur_trellis++;
	}
	
	/* add ldl nodes here */
	f = zSFListMoveFirst(v->sfl);
	while(f != NULL){
		for(j = 0; j < hmm->states;j++){
			v->cache[cur_trellis][0][v->cdna_pos][j] = NULL;
		}
		live_node = zGetTBTreeNode(v->tbtrees[cur_trellis]);
		live_node->pos = f->end;
		live_node->state = f->state;
		live_node->score = 0;
		zTBTreeSetChild(v->tbtrees[cur_trellis]->root,live_node);
		v->cache[cur_trellis][0][v->cdna_pos][f->state] = live_node;
		v->active[cur_trellis] = 0;
		v->start_pos[cur_trellis] = f->end+1;
		v->next[cur_trellis] = cur_trellis+1;
		cur_trellis++;
		f = zSFListMoveNext(v->sfl);
	}	
	zResetSFList(v->sfl);
	v->next[cur_trellis-1] = -1;
	
	v->traceback = zMalloc(sizeof(zPtrList),
												 "zInitViterbiForPin traceback");
	zInitPtrList(v->traceback);
	/* don't use SNP stuff */
	v->snp_count = NULL;
	v->active_snp_count = NULL;
	v->first_snp = NULL;
	v->snp_vals = NULL;
	v->snp_int_vals = NULL;
	v->hapmap = NULL;
	v->hapmax = NULL;
	v->max_active_snp_count = 0;
	v->max_snp_count = 0;
	v->hap_threshold = 0;
}

static void zFreeViterbi(zViterbi* v){
	int i, j;
	coor_t k, m;
	/* free trellis stuff */
	for(i = 0;i < v->max_trellis_count; i++){
		for (k = 0; k < VITERBI_CACHE_LENGTH; k++) {
			for (m = v->cdna_min; m <= v->cdna_max; m++) {
				for(j = 0; j < v->trellis->hmm->states; j++){
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

	zFreeSFList(v->sfl);
	zFree(v->sfl);

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
		for(i = 0;i < v->trellis->dna->seq->var_count; i++){
			for(j = 0; j < v->trellis->dna->seq->variants[i]->variants; j++){
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
		v->trellis->dna->seq->variants[i+f]->cur = vals[j];
		v->trellis->rdna->seq->variants[i+f]->cur = vals[j]; 
		j++;
	}
	zDNASetSNPs(v->trellis->dna);
	zDNASetSNPs(v->trellis->rdna);
}

/* unsets all snps in the DNA seq for this allele */
static void zSNPUnSetSeqForAllele(zViterbi* v,int allele){
	int i,f;
	if(v->snp_count[allele] == 0){
		return;
	}
	f = v->first_snp[allele];
	for(i = 0; i < v->snp_count[allele]; i++){
		v->trellis->dna->seq->variants[i+f]->cur = 0;
		v->trellis->rdna->seq->variants[i+f]->cur = 0; 
	}
	zDNASetSNPs(v->trellis->dna);
	zDNASetSNPs(v->trellis->rdna);
}

/* copies a zTBTree object into a zSFList tracing back from a zTBTreeNode */
static void zTBTree2SFL(zSFList* sfl, zTBTree* tree, zTBTreeNode* n, zHMM* hmm, zTrellis* trellis){
	zSfeature* f = NULL;

	while(n != tree->root){
		f = zSFListAppend(sfl);
		f->name   = hmm->state[n->state].name;
		f->end    = n->pos;
		f->start  = n->parent->pos+1;
		f->strand = hmm->state[n->state].strand;
		f->state  = n->state;
		f->group  = NULL;
		/*scanner       = ('-' == f->strand) 
			? trellis->rscanner[v->trellis->hmm->state[n->state].model]
			: trellis->scanner[v->trellis->hmm->state[n->state].model];
			istate.score  = scanner->scoref(scanner, f);
		*/
		f->score  = n->score;

		zChar2FragFrame(n->frame_data,f);
		
		while(n != tree->root && n->parent->state == n->state){
			/*EVAN heuristic
			int i;
			for(i = 0;i <= 10; i++){
				if(n->ng[i]){
					n->parent->ng[i] = 1;
				}
			}
			*/
			n = n->parent;
		}
		f->start = n->parent->pos+1;
		if(n->parent != tree->root){
			f->score -= (n->parent->score + 
						 zGetTransitionScore(trellis->hmm, n->parent->state, n->state,trellis->tiso_group));
			/*EVAN now removing transition scores from output */
		}
		/*EVAN heuristic
		if(1){
			int i;
			fprintf(stderr,"TRACE");
			for(i = 0;i <= 10; i++){
				fprintf(stderr," %d",n->ng[i]);
			}
			fprintf(stderr,"\n");
		}
		*/
		n = n->parent;
	}
}

/* copies a zTBTree object into a zSFList tracing back from a zTBTreeNode */
static void zBackwardTBTree2SFL(zSFList* sfl, zTBTree* tree, zTBTreeNode* n, zHMM* hmm,zTrellis* trellis){
	zSfeature* f = NULL;

	while(n != tree->root){
		f = zSFListPrepend(sfl);
		f->name   = hmm->state[n->state].name;
		f->start    = n->pos;
		f->end  = n->parent->pos-1;
		f->strand = hmm->state[n->state].strand;
		f->state  = n->state;
		f->group  = NULL;
		/*scanner       = ('-' == f->strand) 
			? trellis->rscanner[v->trellis->hmm->state[n->state].model]
			: trellis->scanner[v->trellis->hmm->state[n->state].model];
			istate.score  = scanner->scoref(scanner, f);
		*/
		f->score  = n->score;
		zChar2FragFrame(n->frame_data,f);
		while(n != tree->root && n->parent->state == n->state){
			n = n->parent;
		}
		f->end = n->parent->pos-1;
		if(n->parent != tree->root){
			f->score -= (n->parent->score + 
						 zGetTransitionScore(trellis->hmm, n->state, n->parent->state,trellis->tiso_group));
			/*EVAN now removing transition score from output */
		}
		n = n->parent;
	}
}

/* Update live nodes for a given feature. */
static int zViterbiAddFeature(zViterbi* v, int t, int fnum, zSfeature *f){
	zTrellis*      trellis = v->trellis;
	int            i,j,state,prestate,poss;
	score_t        score;
	zIVec*         fmap;
	zIVec*         jumps;
	zPhase_t       left_phase,right_phase;
	zTBTreeNode*   live_node;
	zTBTreeNode*   old_p;
	zTBTree*       tree = v->tbtrees[t];
	zTBTreeNode**  cells = v->cache[t][0][v->cdna_pos];
	const char*    state_name;

	int new_node = 0;

	poss = ('-' != f->strand);

	if(poss){
		fmap = trellis->fmap3[fnum];
		right_phase = zSFEndPhase(trellis->dna, f);
		left_phase = zSFBeginPhase(trellis->dna, f);
	}
	else{
		fmap = trellis->rfmap3[fnum];
		right_phase = zSFBeginPhase(trellis->rdna, f);
		left_phase = zSFEndPhase(trellis->rdna, f);
	}
	/* for each (external) state which this feature can end in */
	for(i = 0;i < fmap->size;i++){
		state = fmap->elem[i];
		/* make sure extenral state (state) is compatible with feature */
		if(right_phase != trellis->hmm->state[state].phase) continue;					
		
		live_node = zFindViterbiLiveNode(v,t,f->end,state);
		if(live_node->score == MIN_SCORE){
			new_node = 1;
		}

		/* Select the best path to live_node
			 This path may go through any of the states which transition to to_pad->state */
		jumps = trellis->hmm->jmap[state];
		for(j = 0;j < jumps->size;j++){
			prestate = jumps->elem[j];
			/* make sure internal state (prestate) is compatible with feature */
			if (!zCompatPhases(poss, trellis->hmm->state[prestate].phase, left_phase)) 
				continue;
			if(cells[prestate] == NULL){
				continue;
			}
			/* this stupid check is added for compatability and should be removed
			   once we all agree that checking the last base (and only the last base)
			   of a non-exon external state during the transition computations is stupid */
			if(trellis->estseq != NULL && !trellis->est_para_mode) {
				char c = zGetEstseqSeq(trellis->estseq,f->end); 
				state_name = zStrIdx2Char(trellis->hmm->state[state].name);
				if (c == '1' && state_name[0] != 'E'){
					continue;
				}
			}
			
			/* for each possible state state/pos, if path to end state though prestate 
			   is the best we've seen so far use it */
			score = zScoreExternalState(trellis,state,f->end,f->end-f->start+1,
										f->score) + 
				zGetTransitionScore(trellis->hmm,prestate,state,trellis->tiso_group) +
				cells[prestate]->score;
			if(score > live_node->score){
				live_node->score = score;
				live_node->frame_data = zFragFrame2Char(f->lfrag,f->rfrag,f->frame);
				old_p = live_node->parent;
				zTBTreeSetChild(cells[prestate],live_node);
				if(old_p != NULL && old_p->pos < v->pos-1 && old_p->children == 0){
					/* remove old (now dead) parent */
					zReleaseDeadTBTreeNode(tree,old_p);
				}
			}
		}
		if(live_node->parent == NULL){
			zRemoveViterbiLiveNode(v,t,live_node);
		}
	}
	return new_node;
}

/* Update live nodes for a given feature. */
static int zViterbiAddFeatureBackwards(zViterbi* v, int t, int fnum, zSfeature *f){
	zTrellis*      trellis = v->trellis;
	int            i,j,state,poss,poststate;
	score_t        score;
	zIVec*         fmap;
	zIVec*         postjumps;
	zPhase_t       left_phase,right_phase;
	zTBTreeNode*   live_node;
	zTBTreeNode*   old_p;
	zTBTree*       tree = v->tbtrees[t];
	zTBTreeNode**  cells = v->cache[t][0][v->cdna_pos];
	const char*    state_name;
	
	int new_node = 0;

	poss = ('-' != f->strand);

	if(poss){
		fmap = trellis->fmap3[fnum];
		right_phase = zSFEndPhase(trellis->dna, f);
		left_phase = zSFBeginPhase(trellis->dna, f);
	}
	else{
		fmap = trellis->rfmap3[fnum];
		right_phase = zSFBeginPhase(trellis->rdna, f);
		left_phase = zSFEndPhase(trellis->rdna, f);
	}

	/*fprintf(stderr," ADD %d:  ",f->start);
	  zWriteSfeature(stderr,f);*/
	
	/* for each (external) state which this feature can end in */
	for(i = 0;i < fmap->size;i++){
		score_t base_score;
		state = fmap->elem[i];
		/* make sure extenral state (state) is compatible with feature */
		if(right_phase != trellis->hmm->state[state].phase) continue;					

		base_score = zScoreExternalState(trellis,state,f->end,f->end-f->start+1,
										 f->score);
		
		live_node = zFindViterbiLiveNode(v,t,f->start,state);

		/* for each (internal) state which this feature can transition to */
		/*  select the best path back from live_node */
		postjumps = trellis->hmm->fmap[state];
		for(j = 0;j < postjumps->size;j++){
			poststate = postjumps->elem[j];
			if(cells[poststate] == NULL){
				continue;
			}
			
			/* EVAN
			   this stupid check is added for compatability and should be removed
			   once we all agree that checking the last base (and only the last base)
			   of a non-exon external state during the transition computations is stupid */
			if(trellis->estseq != NULL && !trellis->est_para_mode) {
				char c = zGetEstseqSeq(trellis->estseq,f->end); 
				state_name = zStrIdx2Char(trellis->hmm->state[state].name);
				if (c == '1' && state_name[0] != 'E'){
					continue;
				}
			}
			score = base_score + zGetTransitionScore(trellis->hmm,state,poststate,trellis->tiso_group) +
				cells[poststate]->score;
			
			if(score > live_node->score){
				live_node->phase = left_phase;
				live_node->score = score;
				live_node->frame_data = zFragFrame2Char(f->lfrag,f->rfrag,f->frame);
				old_p = live_node->parent;
				zTBTreeSetChild(cells[poststate],live_node);
				if(old_p != NULL && old_p->pos > v->pos+1 && old_p->children == 0){
					/* remove old (now dead) parent */
					zReleaseDeadTBTreeNode(tree,old_p);
				}
			}
		}
		
		if(live_node->parent == NULL){
			zRemoveViterbiLiveNode(v,t,live_node);
			continue;
		}

	}

	return new_node;
}

/* create all features starting at pos and moving towards the end of the sequence */
static int zSNPViterbiCreateLiveNodes(zViterbi* v,int allele){
	zTrellis *trellis = v->trellis;
	coor_t l, tmp_start;
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
				new_node_count += zViterbiAddFeature(v,allele,j,f);
				f = zSFListMoveNext(sfl);
			}
			zResetSFList(sfl);
		}
		l = trellis->dna->length - v->pos - 1;
		if(trellis->rfactory[j] != NULL){
			trellis->rfactory[j]->create3(trellis->rfactory[j], l, sfl);
			f = zSFListMoveFirst(sfl);
			while(f != NULL){
				f->strand = '-';
				f->state  = j;
				tmp_start = f->start;
				f->start  = trellis->dna->length - f->end - 1;
				f->end    = trellis->dna->length - tmp_start - 1;
				new_node_count += zViterbiAddFeature(v,allele,j,f);
				f = zSFListMoveNext(sfl);
			}
			zResetSFList(sfl);
		}
	}		
	
	return new_node_count;
}

/* create all features starting at pos and moving towards the end of the sequence */
static int zPinViterbiCreateLiveNodes(zViterbi* v){
	zTrellis *trellis = v->trellis;
	coor_t l, tmp_start;
	int i,j;
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
				i = 0;
				while(i >= 0){
					if(v->active[i]){
						new_node_count += zViterbiAddFeature(v,i,j,f);
					}
					i = v->next[i];
				}
				f = zSFListMoveNext(sfl);
			}
			zResetSFList(sfl);
		}
		l = trellis->dna->length - v->pos - 1;
		if(trellis->rfactory[j] != NULL){
			trellis->rfactory[j]->create3(trellis->rfactory[j], l, sfl);
			f = zSFListMoveFirst(sfl);
			while(f != NULL){
				f->strand = '-';
				tmp_start = f->start;
				f->start  = trellis->dna->length - f->end - 1;
				f->end    = trellis->dna->length - tmp_start - 1;
				i = 0;
				while(i >= 0){
					if(v->active[i]){
						new_node_count += zViterbiAddFeature(v,i,j,f);
					}
					i = v->next[i];
				}
				f = zSFListMoveNext(sfl);
			}
			zResetSFList(sfl);
		}
	}		
	return new_node_count;
}

/* create all features starting at pos and moving towards the beginning of the seq */
static int zPinViterbiCreateLiveNodesBackwards(zViterbi* v){
	zTrellis *trellis = v->trellis;
	coor_t l, tmp_start;
	int i,j;
	zSFList* sfl = v->sfl;
	zSfeature *f;
	
	int new_node_count = 0;
	
	for (j = 0; j < trellis->hmm->feature_count; j++) {
		l = v->pos;
		if(trellis->factory[j] != NULL){
			trellis->factory[j]->create3(trellis->factory[j], l, sfl);
			f = zSFListMoveFirst(sfl);
			while(f != NULL){
				f->strand = '+';
				i = 0;
				while(i >= 0){
					if(v->active[i]){
						new_node_count += zViterbiAddFeatureBackwards(v,i,j,f);
					}
					i = v->next[i];
				}
				f = zSFListMoveNext(sfl);
			}
			zResetSFList(sfl);
		}
		l = trellis->dna->length - v->pos - 1;
		if(trellis->rfactory[j] != NULL){
			trellis->rfactory[j]->create5(trellis->rfactory[j], l, sfl);
			f = zSFListMoveFirst(sfl);
			while(f != NULL){
				f->strand = '-';
				tmp_start = f->start;
				f->start  = trellis->dna->length - f->end - 1;
				f->end    = trellis->dna->length - tmp_start - 1;
				i = 0;
				while(i >= 0){
					if(v->active[i]){
						new_node_count += zViterbiAddFeatureBackwards(v,i,j,f);
					}
					i = v->next[i];
				}
				f = zSFListMoveNext(sfl);
			}
			zResetSFList(sfl);
		}
	}		
	return new_node_count;
}

/* Find the maximum position of any postential feature starting at pos
   and in pos's frame */  
static coor_t zSNPGetMaxFeature(zViterbi* v){
	coor_t f_pos;
	/* get furthest in frame stop under any SNP value */
	int sc = zGetStopSeqMaxNextPos(v->trellis->stopseq,v->pos);
	int rsc = v->trellis->dna->length - 1 -
		zGetStopSeqMinPrevPos(v->trellis->rstopseq, v->trellis->dna->length - v->pos -1);
	f_pos = (sc > rsc) ? sc : rsc;
	/* hack to handle fixed length features (promoter and polya) */
	return zIntMax(f_pos,v->pos+PADDING);
}

/*EVAN this is not being used anywhere
static void zSNPDropSNP(zViterbi* v,int snp_idx){
	int max_snp_idx,max_snp_freq_idx,i,j,next_free;
	float max_snp_freq;
	zSeqVariant** vars = v->trellis->dna->seq->variants;

	/ pick the best SNP to drop /
	max_snp_idx = snp_idx;
	max_snp_freq_idx = 0;
	max_snp_freq = vars[max_snp_idx]->freqs[max_snp_freq_idx];
	i = max_snp_idx;
	for(j = 0; j < vars[i]->variants;j++){
		if(vars[i]->freqs[j] > max_snp_freq){
			max_snp_idx = i;
			max_snp_freq_idx = j;
			max_snp_freq = vars[i]->freqs[j];
		}			
	}
	for(i = 0;i < v->snp_count; i++){
		for(j = 0; j < vars[v->snps[i]]->variants;j++){
			if(vars[v->snps[i]]->freqs[j] > max_snp_freq){
				max_snp_idx = v->snps[i];
				max_snp_freq_idx = j;
				max_snp_freq = vars[v->snps[i]]->freqs[j];
			}
		}
	}
	if(max_snp_idx == snp_idx){
		/ we are just going to skip the next SNP so set it to the max freq
			 and move on /
		v->trellis->dna->seq->variants[snp_idx]->cur = max_snp_freq_idx;
		v->trellis->rdna->seq->variants[snp_idx]->cur = max_snp_freq_idx;
		zDNASetSNPs(v->trellis->dna);
		zDNASetSNPs(v->trellis->rdna);
		return;
	}
	else{
		/ we are removing an old SNP /
		next_free = 0;
		for(i = 0;i < v->trellis_count; i++){
			if(v->snp_vals[i][max_snp_idx] == max_snp_freq_idx){
				/ keep i /
				/ find next slot we will be dropping /
				while(v->snp_vals[next_free][max_snp_idx] == max_snp_freq_idx){
					next_free++;
				}
				if(next_free < i){
					/ remove snp from i /
					for(j = 0; j < v->snp_count;j++){
						
					}
					/ move i to next_free /
					zSNPSwapAlleles(v,i,next_free);
					next_free++;
				}
			}
		}
		for(i = next_free;i < v->trellis_count; i++){
			/ clear released allele /
			for(j = 0; j < v->snp_count; j++){
				v->snp_vals[i][j] = 0;
			}
			for(j = 0; j < v->trellis->hmm->states; j++){
				v->cells[i][j] = NULL;
			}
			zResetPtrList(v->live_nodes[i]);
			zResetTBTree(v->tbtrees[i]);			
		}
		v->snp_count = next_free;
	}	
}
*/

/* split each live node list into one for each value of the SNP at snp_idx  */
static void zSNPSplitLiveNodes(zViterbi* v,int snp_idx){
	int a,i,j, new_a, first_copy, prev_val;
	zSeqVariant** vars = v->trellis->dna->seq->variants;
	zSeqVariant*  var  = vars[snp_idx];
	
	fprintf(stderr,"SPLIT at pos: %d\tsnp(%d): %d\n",v->pos,snp_idx,var->pos);

	/* allele 0 starts at its last cpoint and we don't care about snps
	 in allele 0 before last cpoint (they are all 0) */
	v->start_pos[0] = v->tbtrees[0]->cpoint->pos;
	while(v->snp_count[0] > 0 && 
		  vars[v->first_snp[0]]->pos < v->tbtrees[0]->cpoint->pos){
		v->first_snp[0]++;
		v->snp_count[0]--;
	}

	/*zSNPPrintAlleles(v);*/

	new_a = v->trellis_count;
	/* if we currently are dealing with only the reference sequence (no "active" snps) */
	if(v->trellis_count == 1){
		v->first_snp[0] = snp_idx;
		v->snp_count[0] = 0;
	}
	/* need to create a new block of v->max_snp_count alleles */
	for(a = 0; a < v->trellis_count; a++){
		if(v->snp_count[a] < v->max_snp_count){
			/* this allele can handle another snp so add it and make a copy */
			/* make a new allele for each value of new SNP and old allele*/
			if(v->snp_count[a] >= 1){
				/* if this allele has a value for the preceeding snp */
				prev_val = v->snp_vals[a][v->snp_count[a]-1];
			}
			else{
				/* otherwise just assume the reference sequence value */
				prev_val = 0;
			}
			/* up the snp count for this allele */
			v->snp_count[a]++;
			if(snp_idx > 0 && a > 0){
				/* if we have hap info for precedding snp (not first snp) and this is not the 
				   0 allele (which is restricted to having all 0 snp values) */
				v->snp_vals[a][v->snp_count[a]-1] = v->hapmax[snp_idx-1][prev_val];
				if(v->hapmax[snp_idx-1][prev_val] != 0){
					/* if we use a non zero value increase the active snp count for this allele */
					v->active_snp_count[a]++;
				}
			}
			else{
				/* if this is the first snp then there is no hap info for preceding snp */
				v->snp_vals[a][v->snp_count[a]-1] = 0;
			}
			if(v->active_snp_count[a] >= v->max_active_snp_count){
				continue;
			}
			for(i = 0; i < var->variants; i++){
				/*fprintf(stderr,"\tFULLSPLIT %d -> %d\n",a,new_a);*/
				if(i == v->snp_vals[a][v->snp_count[a]-1]){
					/* don't need to make a copy for the snp value we used above */ 
					continue;
				}
				if(snp_idx > 0 && v->hapmap[snp_idx-1][prev_val][i] <= v->hap_threshold){
					/* if this combination of values is not common enough then skip it */
					continue;
				}
				zSNPCopyAllele(v,a,new_a);
				v->active[new_a] = 1;
				v->snp_vals[new_a][v->snp_count[new_a]-1] = i;
				v->active_snp_count[new_a]++;
				new_a++;
			}
		}
		else{
			/* this allele cannot handle another snp so we only create a copy when
			   all snps preceding snp_idx - v->max_snp_count are 0 */
			if(a > 0){
				/* if newly created allele will not overlap the snps in this allele
				   at all and this allele is not at index 0 it will contain a 1 
				   before (snp_idx - v->max_snp_count) so skip it */
				if(v->first_snp[a] + v->max_snp_count - 1 < 
				   snp_idx - v->max_snp_count + 1){
					continue;
				}
				else{
					int no_copy = 0;
					/* check that all snps preceding (snp_idx - v->max_snp_count) 
					   are 0 */ 
					for(i = 0;i <= snp_idx - v->max_snp_count - v->first_snp[a];i++){
						if(v->snp_vals[a][i] != 0){
							no_copy = 1;
							break;
						}
					}
					if(no_copy) continue;
				}
			}
			if(v->active_snp_count[a] >= v->max_active_snp_count){
				continue;
			}
			prev_val = 0;
			if(v->first_snp[a] + v->snp_count[a] == snp_idx){
				prev_val = v->snp_vals[a][v->snp_count[a]-1];
			}
			/* start at 1 since the 0 value at snp_idx is implicitly used by the snp
			   we are copying from */ 
			for(i = 1; i < var->variants; i++){
				/*fprintf(stderr,"\tPARTSPLIT %d -> %d\n",a,new_a);*/
				if(v->hapmap[snp_idx-1][prev_val][i] <= v->hap_threshold){
					/* if this combination of values is not common enough then skip it */
					continue;
				}
				/* copy allele */
				zSNPCopyAllele(v,a,new_a);
				v->active[new_a] = 1;
				/* set first snp */
				v->first_snp[new_a] = snp_idx-v->max_snp_count+1;
				v->snp_count[new_a] = v->max_snp_count;
				v->active_snp_count[new_a]++;
				/* set new snp val */
				v->snp_vals[new_a][v->max_snp_count-1] = i;
				/* copy in old snp vals in appropriate slots*/
				first_copy = v->first_snp[new_a] - v->first_snp[a];
				for(j = first_copy; j < v->max_snp_count; j++){
					v->snp_vals[new_a][j-first_copy] = v->snp_vals[a][j];
				}
				/* fill in 0 for snps which weren't copied */
				for(j = zIntMax(v->max_snp_count-first_copy,0);
					j < v->max_snp_count-1; j++){
					v->snp_vals[new_a][j] = 0;
				}
				/* next index */
				new_a++;
			}
		}
	}

	fprintf(stderr,"  ALLELECOUNT %d->%d\n",v->trellis_count,new_a);
	v->trellis_count = new_a;

	/*zSNPPrintAlleles(v);*/
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

 
static void zSNPTraceAllele(zViterbi* v, int a, bool cptrace, score_t base_score, zSFList* base_list,char* text){
	zAlleleDecode* ad;
	zTBTreeNode* tbnode;
	zSfeature* f;
	score_t best_score;
	int i,state;
	zSeqVariant* var;
	char* new_text;/*EVAN*/
	zTBTreeNode** cells = v->cache[a][0][v->cdna_pos];

	fprintf(stderr,"TRACE: pos %d\t%d",v->pos,a);
	for(i = 0; i < v->snp_count[a]; i++){
		var = v->trellis->dna->seq->variants[v->first_snp[a]+i];	
		fprintf(stderr,"\t%d=%c(%d)",var->pos,var->values[v->snp_vals[a][i]],
				v->snp_vals[a][i]);
	}
	fprintf(stderr,"\n");

	
	if(cptrace == true){
		tbnode = v->tbtrees[a]->cpoint;
	}
	else{
		tbnode = NULL;
		best_score = MIN_SCORE;
		for(state = 0;state < v->trellis->hmm->states;state++){
			if(cells[state] != NULL && 
				 cells[state]->score > best_score){
				tbnode = cells[state];
				best_score = cells[state]->score;
			}
		}
		if(tbnode == NULL){
			zDie("No live states at end of trellis.  This should never happen");
		}
	}


	ad = zMalloc(sizeof(zAlleleDecode),"zSNPCollapsePaths ad");
	zInitAlleleDecode(ad,(v->snp_count[a]+1)*100,v->snp_count[a]);
	ad->snp_count = v->snp_count[a];
	ad->first_snp = v->first_snp[a];
	ad->start = v->tbtrees[a]->root->pos+1;
	ad->stop = tbnode->pos;
	ad->score_diff = tbnode->score - base_score;
	for(i = 0; i < ad->snp_count; i++){
		ad->vals[i] = v->snp_vals[a][i];
		var = v->trellis->dna->seq->variants[i+ad->first_snp];
	}		

	/* traceback for this allele */
	zTBTree2SFL(ad->sfl,v->tbtrees[a],tbnode,v->trellis->hmm,v->trellis);
	f = zSFListMoveLast(ad->sfl);
	if(f->start < PADDING){
		f->score -= zGetInitProb(v->trellis->hmm,f->state,v->trellis->iiso_group);
		f->start = PADDING;
	}
	if(base_list != NULL){
		zSnpViterbiMergeSFLists(ad->sfl,base_list);
	}

	v->tbt_use_count[a]--;
	
	new_text = zMalloc(sizeof(char)*(strlen(text)+10),"");
	sprintf(new_text,"%s -> %d",text,a);
	for(i = 1; i < v->trellis_count; i++){
		if(v->next_tbti[i] == a){
			/* this math gives the correct score offset, the score of the last pos in i
			   minus this score will give the overall score_offset */
			zSNPTraceAllele(v,i,true,v->next_tbtn[i]->score - ad->score_diff,ad->sfl,new_text);
		}
	}
	zFree(new_text);
	
	if(v->next_tbtn[a] != NULL){
		v->tbt_use_count[v->next_tbti[a]]--;
		v->next_tbtn[a] = NULL;
		v->next_tbti[a] = 0;
	}
	v->active[a] = 0;

	if(v->tbt_use_count[a] != 0){
		for(i = 1; i < v->max_trellis_count; i++){
			if(v->next_tbti[i] == a){
				if(i >= v->trellis_count){
					fprintf(stderr,"Allele %d should have no refs but %d (>= trellis_count=%d) refs it\n",a,v->trellis_count,i);
				}
				else{
					fprintf(stderr,"Allele %d should have no refs but %d refs it\n",a,i);
				}
			}
		}
		fprintf(stderr,"well, shit. tbt_use_count[%d](%d)-- = %d\n",a,v->tbtrees[a]->id,v->tbt_use_count[a]);
		zDie("This should never happen. %d (%d)",a,v->tbt_use_count[a]);
	}

	zTranslateSFList(ad->sfl,-PADDING);
	
	/* add traceback to list */
	zPtrListAddLast(v->traceback,ad);
}
	
/*EVAN problem here when a cpoint in one tree matches a node coming before
  (closer to the root) the cpoint in another tree.  I don't know if it 
  happens but it is not being checked for.  The problem goes away at the next 
  cpoint since both trees will share it, but potentially wastes computation.
  I don't know if this is worth checking for (in terms of speeding up the 
  computation) or if it will clow things down */
static void zSNPCollapsePaths(zViterbi* v){
	coor_t cpos; /* cpoint position */
	int cstate; /* cpoint state */
	int a,i,first_snp,last_live;
	bool collapsed = false;

	zSeqVariant** vars = v->trellis->dna->seq->variants;	

	if(any_new_cpoint == 0){
		return;
	}
	cpos = any_new_cpoint;
	any_new_cpoint = 0;

	first_snp = v->first_snp[0];
	if(first_snp >= v->trellis->dna->seq->var_count){
		return;
	}

	v->cpos = cpos;/*EVAN*/

	if(cpos == (coor_t)-1 || first_snp < 0){
		/* cpoint comes before first snp */
		return;
	}
	if(cpos < vars[first_snp]->pos){
		return;
	}

	for(a = 0; a < v->trellis_count; a++){
		int a1,a2,b;
		if(v->active[a] == 0){
			continue;
		}
		if(!(zTBTreeCheckNewCpoint(v->tbtrees[a]))){
			continue;
		}
		cpos = v->tbtrees[a]->cpoint->pos;
		cstate = v->tbtrees[a]->cpoint->state;
		for(b = 0; b < v->trellis_count; b++){
			if(v->active[a] == 0){
				break;
			}
			if(b == a || v->active[b] == 0){
				continue;
			}
			if((v->tbtrees[b]->cpoint->pos == cpos) &&
			   (v->tbtrees[b]->cpoint->state == cstate)){
				bool skip = false;
				a1 = zIntMin(a,b);
				a2 = zIntMax(a,b);
				for(i = 0; i < v->snp_count[a1];i++){
					/* if we are pre-cpos it doesn't matter what value the snp has */
					if(vars[v->first_snp[a1]+i]->pos < cpos){
						continue;
					}
					/* needs to be a zero if we don't overlap a2 and come after cpos */
					if((v->first_snp[a1]+i < v->first_snp[a2]) ||
					   (v->first_snp[a1]+i >= v->first_snp[a2]+v->snp_count[a2])){
						if(v->snp_vals[a1][i] != 0){
							skip = true;
							break;
						}
					}
					/* need to match other a2 if we overlap it and come  after cpos */
					else{
						if(v->snp_vals[a1][i] != 
						   v->snp_vals[a2][i+v->first_snp[a1]-v->first_snp[a2]]){
							skip = true;
							break;
						}
					}
				}					
				if(skip){
					continue;
				}
				for(i = 0; i < v->snp_count[a2];i++){
					/* if we are pre-cpos it doesn't matter what value the snp has */
					if(vars[v->first_snp[a2]+i]->pos < cpos){
						continue;
					}
					/* needs to be a zero if we don't overlap a1 and come after cpos */
					if((v->first_snp[a2]+i < v->first_snp[a1]) ||
					   (v->first_snp[a2]+i >= v->first_snp[a1]+v->snp_count[a1])){
						if(v->snp_vals[a2][i] != 0){
							skip = true;
							break;
						}
					}
					/* need to match a1 if we overlap it and come after cpos */
					else{
						if(v->snp_vals[a2][i] != 
						   v->snp_vals[a1][i+v->first_snp[a2]-v->first_snp[a1]]){
							skip = true;
							break;
						}
					}
				}					
				if(skip){
					continue;
				}
				/*
				int pre_snp = 0;
				int diff_snp = 0;
				for(i = 0; i < v->snp_count[a2]; i++){
					if(vars[v->first_snp[a2]+i]->pos < v->pos){
						pre_snp++;
						if((v->first_snp[a2]+i < v->first_snp[a1]) ||
						   (v->first_snp[a2]+i >= v->first_snp[a1]+v->snp_count[a1])){
							if(v->snp_vals[a2][i] != 0){
								diff_snp++;
							}
						}
						else{
							if(v->snp_vals[a2][i] != 
							   v->snp_vals[a1][i+v->first_snp[a2]-v->first_snp[a1]]){
								diff_snp++;
							}
						}
					}
				}
				*/
				/* the snps after cpos all match so we can collapse a2 */
				/* we have collapsed to allele 0 so we can clean up and trace this allele */ 
				if(a1 == 0){
					/*fprintf(stderr,"COLLAPSE: %d,%d -> %d\t%d of %d, %d diff\n",v->cp[a2],v->sp[a2],v->pos,
							pre_snp,v->snp_count[a2],diff_snp);EVAN*/
					zSNPCollapseAllele(v,a2);
					collapsed = true;
				}
				/* we have collapsed to a different allele so we can link to the tbtree and
				   quit doing viterbi work for this allele */
				else{
					/*fprintf(stderr,"MATCH: %d,%d -> %d\t%d of %d, %d diff\n",v->cp[a2],v->sp[a2],v->pos,
							pre_snp,v->snp_count[a2],diff_snp);EVAN*/
					/* link to a1's tbtree, increate a1's use_count, and deactivate a2 */	
					if(v->next_tbti[a2] != 0 ||
					   v->next_tbtn[a2] != NULL){
						zDie("");
					}
					v->next_tbtn[a2] = v->tbtrees[a1]->cpoint;
					v->next_tbti[a2] = a1;
					zTBTreeLockNode(v->tbtrees[a1]->cpoint);
					v->tbt_use_count[a1]++;
					/*fprintf(stderr,"set tbt_use_count[%d](%d)++ = %d\n",a1,v->tbtrees[a1]->id,v->tbt_use_count[a1]);EVAN*/

					zClearTBTreeNodeChildren(v->tbtrees[a2],v->tbtrees[a2]->cpoint);
					v->active[a2] = 0;
					/* remove the tbtree below the cpoint in a2 and clear the cells and 
					   live_nodes data */
					zTBTreeClearSubTree(v->tbtrees[a2],v->tbtrees[a2]->cpoint);
					for(i = 0;i < v->trellis->hmm->states; i++){
						v->cache[a2][0][v->cdna_pos][i] = NULL;
					}
					v->cache[a2][0][v->cdna_pos][v->tbtrees[a2]->cpoint->state] = v->tbtrees[a2]->cpoint;
					zResetPtrList(v->live_nodes[a2]);
				}
			}
		}
	}

	if(!collapsed){
		return;
	}

	/* collapse data to lowest indices. */
	last_live = v->trellis_count-1;
	v->temp[0] = 0;
	for(a = 1; a < v->trellis_count; a++){
		v->temp[a] = a;
	}
	for(a = 1; a < v->trellis_count; a++){
		if(last_live < a){
			break;
		}
		if(v->tbt_use_count[a] == 0){
			/*fprintf(stderr,"\trelease %d\n",a);EVAN*/
			zSNPReleaseAllele(v,a);
			while(v->tbt_use_count[last_live] == 0){
				/*fprintf(stderr,"\trelease %d\n",last_live);EVAN*/
				zSNPReleaseAllele(v,last_live);
				last_live--;	
			}
			if(a < last_live){
				int b;
				/*fprintf(stderr,"\tswap %d <-> %d\n",a,last_live);EVAN*/
				zSNPSwapAlleles(v,a,last_live);
				for(b = 1; b < v->trellis_count; b++){
					if(v->next_tbti[b] == last_live){
						v->next_tbti[b] = a;
					}
				}
				v->temp[last_live] = a;
				last_live--;
			}
		}
	}
	/*fprintf(stderr,"  collapse alleles\t%d->%d\n",v->trellis_count,last_live+1);EVAN*/
	v->trellis_count = last_live+1;

	for(a = 0; a < v->trellis_count; a++){
		zTBTreeClearNewCpoint(v->tbtrees[a]);
	}
}

zSFVec* zRunViterbi(zTrellis* trellis, score_t *path_score){
	zPtrList *list;
	zAlleleDecode *ad;
	zSFVec* sfv;
	zSfeature* f;
	if(trellis->dna->seq->var_count > 0){
		zDie("Use zRunSNPViterbi not zRunViterbi on sequences containing snps.");
	}
	list = zRunSNPViterbi(trellis,path_score);
	ad = zPtrListMoveFirst(list);
	sfv = zMalloc(sizeof(zSFVec),"zRunViterbi sfv");
	zInitSFVec(sfv,ad->sfl->size);
	f = zSFListMoveFirst(ad->sfl);
	while(f != NULL){
		zPushSFVec(sfv,f);
		f = zSFListMoveNext(ad->sfl);
	}
	zFreeAlleleDecode(ad);
	zFree(ad);
	zFreePtrList(list);
	zFree(list);
	return sfv;
}
	
zPtrList* zRunSNPViterbi(zTrellis* trellis, score_t *path_score){
	zViterbi*     viterbi;
	zHMM*         hmm = trellis->hmm;
	zDNA*         dna = trellis->dna;
	zAlleleDecode* ad;
	score_t       score,best_score;
	int           state,best_state,prestate,i;
	int           j;        /* iterator for internal states */
	zIVec         *jumps;
	zTBTreeNode  *tbtn;
	coor_t        max_feature;
	coor_t        sc[3];
	coor_t        next_snp_pos;
	int           snp_idx;
	int           allele;
	zTBTree*      tree;
	zPtrList*     nodes;
	zTBTreeNode** cells;
	coor_t        cpos;
	int           cstate;
	zPtrList*     ret_list;
	zSfeature*    f;

	/* Prepare Viterbi Vars */
	viterbi = zMalloc(sizeof(zViterbi),"zRunViterbi viterbi");
	zInitViterbi(viterbi,trellis);

	/* Init allele 0 */
	allele = 0;
	viterbi->first_snp[0] = -1;
	viterbi->snp_count[0] = 0;
	viterbi->snp_int_vals[0] = 0;
	tree = viterbi->tbtrees[0];
	nodes = viterbi->live_nodes[0];
	cells = viterbi->cache[0][0][viterbi->cdna_pos];
	viterbi->trellis_count = 1;
	viterbi->tbt_use_count[0] = 1;
	viterbi->start_pos[0] = PADDING-1;
	viterbi->active[0] = 1;

	/* Init allele 0 cells */
	for(j = 0; j < hmm->states;j++){
		score = zGetInitProb(trellis->hmm, j,trellis->iiso_group);
		if(score > MIN_SCORE){
			tbtn = zGetTBTreeNode(tree);
			tbtn->pos = PADDING-1;
			tbtn->state = j;
			tbtn->score = score;
			zTBTreeSetChild(tree->root,tbtn);
			cells[j] = tbtn;
		}
		else{
			cells[j] = NULL;
		}
	}

	/* update stop codon tracking */
	max_feature = 0;
	for(i = 0;i <= 2; i++){
		viterbi->pos = PADDING+i;
		j = (PADDING+i) % 3;
		sc[j] = zSNPGetMaxFeature(viterbi);
		if(sc[j] > max_feature){
			max_feature = sc[j];
		}
	}

	/* init snp tracking */
	snp_idx = 0;
	if(viterbi->trellis->dna->seq->var_count == 0){
		next_snp_pos = viterbi->trellis->dna->length+1;
	}
	else{
		next_snp_pos = viterbi->trellis->dna->seq->variants[0]->pos;
	}
	
	cpos = -1;
	cstate = -1;

	/* walk forward through sequence */
	for(viterbi->pos = PADDING; viterbi->pos < dna->length-PADDING;viterbi->pos++){
		
		/* update feature-snp tracking */
		j = viterbi->pos % 3; 
		sc[j] = zSNPGetMaxFeature(viterbi);
		if(sc[j] > max_feature){
			max_feature = sc[j];
		}
		
		if(max_feature >= next_snp_pos){
			/* split the live node lists */
			zSNPSplitLiveNodes(viterbi,snp_idx);

			/* watch for next snp */
			snp_idx++;
			if(viterbi->trellis->dna->seq->var_count <= snp_idx){
				next_snp_pos = viterbi->trellis->dna->length + 1;
			}
			else{
				next_snp_pos = viterbi->trellis->dna->seq->variants[snp_idx]->pos;
			}
		}
		/* check for collapsing stuff */
		zSNPCollapsePaths(viterbi);

		for(allele = 0; allele < viterbi->trellis_count; allele++){
			if(viterbi->active[allele] == 0) continue;

			zSNPSetSeqForAllele(viterbi,allele);

			/* STEP 1 - create new live nodes */
			zSNPViterbiCreateLiveNodes(viterbi,allele);
			
			/* STEP 2 - step one pos forward */			
			tree = viterbi->tbtrees[allele];
			nodes = viterbi->live_nodes[allele];
			cells = viterbi->cache[allele][0][viterbi->cdna_pos];
			
			for(state = 0;state < hmm->states;state++){
				/* clear out the temp array */
				viterbi->cache[allele][1][viterbi->cdna_pos][state] = NULL;
				/* handle INTERNAL state transitions */ 
				if(hmm->state[state].type == INTERNAL || hmm->state[state].type == GINTERNAL){
					best_state = -1;
					best_score = MIN_SCORE;
					jumps = hmm->jmap[state];
					for(i = 0; i < jumps->size; i++) {
						prestate = jumps->elem[i];
						if(cells[prestate] == NULL) continue;
						
						score = zScoreInternalState(trellis,state,viterbi->pos,
													(viterbi->pos == PADDING || state != prestate));
						
						if(hmm->state[state].type == INTERNAL){
							score += zGetTransitionScore(trellis->hmm,prestate,state,trellis->tiso_group);
						}
						else{
							score += zGetIntergenicContinueScore(trellis);							
						}

						if(score != MIN_SCORE){
							score += cells[prestate]->score;
							if(score > best_score){
								best_score = score;
								best_state = prestate;
							}
						}
					}
					if(best_score != MIN_SCORE){
						/* fill in traceback stuff here */
						tbtn = zGetTBTreeNode(tree);
						viterbi->cache[allele][1][viterbi->cdna_pos][state] = tbtn;
						tbtn->pos = viterbi->pos;
						tbtn->state = state;
						tbtn->score = best_score;
						tbtn->frame_data = 0;
						zTBTreeSetChild(cells[best_state],tbtn);
					}
				}
				/* EXTERNAL states dealt with in when incorporating live nodes below */
				else if(hmm->state[state].type == EXTERNAL){
				}
				else{
					zDie("zRunViterbi only supports INTERNAL and EXTERNAL states");
				}
			}
			/* cleanup unused/unneeded old_cells */
			for(state = 0;state < hmm->states;state++){
				tbtn = cells[state];
				if(tbtn == NULL) continue;
				if(tbtn->children == 0){
					zReleaseDeadTBTreeNode(tree,tbtn);
				}
				else if((tbtn->children == 1) &&
								(tbtn->child->state == tbtn->state)){
					zReleaseRedundantTBTreeNode(tree,tbtn);
				}
			}
			
			/* STEP 3 - move tmp_cells to cells */
			for(state = 0;state < hmm->states;state++){
				cells[state] = viterbi->cache[allele][1][viterbi->cdna_pos][state];
			}
			
			/* STEP 4 - remove/incorporate old live nodes */
			zViterbiIncorporateLiveNodes(viterbi,allele);			

			zSNPUnSetSeqForAllele(viterbi,allele);

			/* EVAN heuristic checking */
			/*zTBTreeFillDescendantScores(tree->cpoint);*/
			/*zTBTreeCheckHeuristicCpoints(viterbi,tree,cells,nodes);*/
			
		}
	}
	
	/*zPrintNodeGraph();EVAN*/

	/* Viterbi trace back */
	zTrace2("traceback\n");
	*path_score = 0;

	cells = viterbi->cache[0][0][viterbi->cdna_pos];
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
	
	/* pick the best node to traceback from for each allele*/
	for(allele = 1; allele < viterbi->trellis_count; allele++){
		if(viterbi->active[allele]){
			char* text = zMalloc(sizeof(char)*10,"");
			sprintf(text,"%d",allele);
			/*EVAN 
			  if(1){
				int diff_snp = 0;
				for(i = 0; i < viterbi->snp_count[allele]; i++){
					if(viterbi->snp_vals[allele][i] != 0){
						diff_snp++;
					}
				}
				fprintf(stderr,"FINAL COLLAPSE: %d,%d -> %d\t%d of %d, %d diff\n",
						viterbi->cp[allele],viterbi->sp[allele],viterbi->pos,
						viterbi->snp_count[allele],viterbi->snp_count[allele],diff_snp);
			}
			*/
			zSNPTraceAllele(viterbi,allele,false,best_score,NULL,text);
			zFree(text);
		}
	}
	viterbi->trellis_count = 1;
	
	ad = zMalloc(sizeof(zAlleleDecode),"zRunSNPViterbi ad");
	zInitAlleleDecode(ad,25,1);
	ad->vals[0] = 0;
	
	sprintf(ad->header,"## Main SNP Sequence\n"); 
	
	zTBTree2SFL(ad->sfl,viterbi->tbtrees[0],cells[best_state],viterbi->trellis->hmm,viterbi->trellis);
	f = zSFListMoveLast(ad->sfl);
	if(f->start < PADDING){
		f->score -= zGetInitProb(viterbi->trellis->hmm,f->state,viterbi->trellis->iiso_group);
		f->start = PADDING;
	}
	zTranslateSFList(ad->sfl, -PADDING);
	zPtrListAddFirst(viterbi->traceback,ad);
	zSNPCleanUpDecodes(viterbi);
	ret_list = viterbi->traceback;
	viterbi->traceback = NULL;
	zFreeViterbi(viterbi);
	zFree(viterbi);
	return ret_list;
}

void zPinViterbiCollapsePaths(zViterbi* v){
	int i,j,c,p;
	/*EVAN problem here when a cpoint in one tree matches a node coming before
	  (closer to the root) the cpoint in another tree.  I don't know if it 
	  happens but it is not being checked for.  The problem goes away at the next 
	  cpoint since both trees will share it, but potentially wastes computation */
	
	if(any_new_cpoint == 0){
		return;
	}
	any_new_cpoint = 0;
	
	/* when we start we know no two cpoints are the same unless at least one
	   of their new cpoint flags are set */
	i = 0;
	while(i >= 0){
		if(v->active[i] == 0){
			i = v->next[i];
			continue;
		}
		if(zTBTreeClearNewCpoint(v->tbtrees[i])){
			j = 0;
			while(j >= 0){
				if(v->active[j] == 0){
					j = v->next[j];
					continue;
				}
				else if(j == i){
					j = v->next[j];
					continue;
				}
				/* if these two cpoints are the same */
				if(v->tbtrees[i]->cpoint->pos == v->tbtrees[j]->cpoint->pos &&
				   v->tbtrees[i]->cpoint->state == v->tbtrees[j]->cpoint->state){
					/* collapse the higher of the two */
					c = zIntMax(i,j);
					v->active[c] = 0;
					p = c-1;
					/* this will never get < 0 since we always collapse the higher index
					   (we never collapse index 0) */
					while(v->next[p] == -1){
						p--;
					}
					/* sanity check */
					if(v->next[p] != c){
						zDie("zViterbi next array got messed up.  This should never happen");
					}
					v->next[p] = v->next[c];
					v->next[c] = -1;
					v->active_count--;
					v->dead_count++;
					if(c == i){
						i = p;
					}
					break;
				}
				j = v->next[j];
			}
		}
		i = v->next[i];
	}
}

static void zRunBackwardViterbi(zViterbi* v, int start_state, coor_t start_pos, 
								coor_t end_pos, score_t* path_score){

	zTBTree*      tree = v->tbtrees[0];
	zPtrList*     nodes = v->live_nodes[0];
	zTBTreeNode** cells = v->cache[0][0][v->cdna_pos];
	zTBTreeNode*  live_node;
	zTBTreeNode*  tbtn;
	coor_t        cpos;
	int           cstate,state,best_state,i,poststate;
	score_t       best_score,score;
	zIVec*        jumps;
	zSfeature*    f;

	v->pos = start_pos;
	live_node = zGetTBTreeNode(tree);
	live_node->pos = start_pos;
	live_node->state = start_state;
	live_node->score = 0;
	tree->root->pos = start_pos+1;
	zTBTreeSetChild(tree->root,live_node);
	tree->cpoint = live_node;
	cells[start_state] = live_node;

	cpos = -1;
	cstate = -1;

	for(v->pos = start_pos; v->pos > PADDING; v->pos--){
		if(v->pos % 10000 == 0){
			fprintf(stderr, "B_VITERBI %d cpoint0 at %d from %d, %d live nodes\n",
					v->pos,cpos,v->pos,nodes->size);
		}
		
		if(cpos != tree->cpoint->pos || 
		   cstate != tree->cpoint->state){
			cpos = tree->cpoint->pos;
			cstate = tree->cpoint->state;
			fprintf(stderr,"B_CPOINT at %d (%d) from %d\n",cpos,cstate,v->pos);
		}
		
		if(cpos < end_pos){
			break;
		}

		/* STEP 1 - create new live nodes */
		zPinViterbiCreateLiveNodesBackwards(v);
		
		/* STEP 2 - step one pos back */
		for(state = 0;state < v->trellis->hmm->states;state++){
			if(v->trellis->hmm->state[state].type != INTERNAL &&
			   v->trellis->hmm->state[state].type != GINTERNAL) continue;
			/* clear out the temp array */
			v->cache[0][1][v->cdna_pos][state] = NULL;
			best_state = -1;
			best_score = MIN_SCORE;
			jumps = v->trellis->hmm->fmap[state];
			for(i = 0; i < jumps->size; i++) {
				poststate = jumps->elem[i];
				if(cells[poststate] == NULL){
					continue;
				}
				if(cells[poststate]->score == MIN_SCORE) continue;
				
				if(v->trellis->hmm->state[poststate].type == EXTERNAL){
					/*EVAN is this check redundant?  Don't all states have a fixed phase?
					  If so, why would transitions which break the phase have a non-zero
					  probability? */
					if (!zCompatPhases(v->trellis->hmm->state[poststate].strand == '+'?1:0, 
									   v->trellis->hmm->state[state].phase,
									   cells[poststate]->phase)) 
						continue;
				}
				
				if(v->trellis->hmm->state[poststate].type == GINTERNAL){
					/*EVAN this ignores the actual transition score for the first 
					  base but it is this way for compatability with zGInternalTrans
					  in zTransition.c
					*/
					score = zGetIntergenicContinueScore(v->trellis);
				}
				else{
					score = zGetTransitionScore(v->trellis->hmm,state,poststate,v->trellis->tiso_group);
				}
				if(score == MIN_SCORE) continue;

				score += zScoreInternalState(v->trellis,state,v->pos,(state != poststate));
				if(score == MIN_SCORE) continue;
				
				score += cells[poststate]->score;
				if(score > best_score){
					best_score = score;
					best_state = poststate;
				}
			}
			if(best_score != MIN_SCORE){
				/* fill in traceback stuff here */
				tbtn = zGetTBTreeNode(tree);
				v->cache[0][1][v->cdna_pos][state] = tbtn;
				tbtn->pos = v->pos;
				tbtn->state = state;
				tbtn->score = best_score;
				tbtn->frame_data = 0;
				zTBTreeSetChild(cells[best_state],tbtn);
			}
		}
		/* cleanup unused/unneeded old_cells */
		for(state = 0;state < v->trellis->hmm->states;state++){
			tbtn = cells[state];
			if(tbtn == NULL) continue;
			if(tbtn->children == 0){
				zReleaseDeadTBTreeNode(tree,tbtn);
			}
			else if((tbtn->children == 1) &&
							(tbtn->child->state == tbtn->state)){
				zReleaseRedundantTBTreeNode(tree,tbtn);
			}
		}
		
		/* STEP 3 - move tmp_cells to cells */
		for(state = 0;state < v->trellis->hmm->states;state++){
			cells[state] = v->cache[0][1][v->cdna_pos][state];
		}
		
		/* STEP 4 - remove/incorporate old live nodes */	
		zBackwardViterbiIncorporateLiveNodes(v,0);			
	}

	/* trace and return path */
	tbtn = tree->cpoint;
	if(end_pos < tbtn->pos){
		best_state = -1;
		best_score = MIN_SCORE;
		for(state = 0;state < v->trellis->hmm->states;state++){
			if(cells[state] != NULL && 
			   cells[state]->score + zGetInitProb(v->trellis->hmm, state, v->trellis->iiso_group) > best_score){
				best_state = state;
				best_score = cells[state]->score + zGetInitProb(v->trellis->hmm, state,v->trellis->iiso_group);
			}
		}
		tbtn = cells[best_state];
	}
	*path_score = tbtn->score;
	
	zResetSFList(v->sfl);
	zBackwardTBTree2SFL(v->sfl,tree,tbtn,v->trellis->hmm,v->trellis);
		
	f = zSFListMoveLast(v->sfl);
	while(f != NULL && f->end < end_pos){
		zSFListRemoveLast(v->sfl);
		f = zSFListMoveLast(v->sfl);
	}
}

zSFVec* zRunPinViterbi(zTrellis* trellis, score_t *path_score,
					   coor_t start_pos,coor_t end_pos){
	zViterbi*     viterbi;
	zHMM*         hmm = trellis->hmm;
	zDNA*         dna = trellis->dna;
	score_t       score,best_score;
	int           state,best_state,prestate,i;
	int           pin_state = -1;
	zIVec        *jumps;
	zTBTreeNode  *tbtn;
	int           t;
	zTBTree*      tree;
	zPtrList*     nodes;
	zTBTreeNode** cells;
	coor_t        cpos;
	int           cstate;
	zSFVec*       sfv;
	coor_t        pin_pos = (coor_t)-1;
	score_t       pin_score = MIN_SCORE;
	zSfeature*    f;
	zSFList*      temp_sfl;
	score_t       pre_path_score = 0;
	clock_t start, end;
	double cpu_time_used;
	start = clock();

	if(start_pos == 0){
		start_pos = 1;
	}

	/*adjust for padding, 0 index */
	start_pos += PADDING-1;
	end_pos += PADDING-1;

	/* Prepare Viterbi Vars */
	viterbi = zMalloc(sizeof(zViterbi),"zRunViterbi viterbi");
	zInitViterbiForPin(viterbi,trellis,start_pos);
	
	cpos = -1;
	cstate = -1;
	
	/* walk forward through sequence */
	fprintf(stderr, "VITERBI %d, %d live nodes, %d active, %d dead\n",
			viterbi->pos,viterbi->live_nodes[0]->size,
			viterbi->active_count,viterbi->dead_count);
	for(viterbi->pos = start_pos+1; viterbi->pos < dna->length-PADDING;viterbi->pos++){
		
		if(viterbi->pos % 10000 == 0){
			fprintf(stderr, 
					"VITERBI %d cpoint0 at %d from %d, %d live nodes, %d active, %d dead\n",
					viterbi->pos,cpos,viterbi->pos,viterbi->live_nodes[0]->size,
					viterbi->active_count, viterbi->dead_count);
		}
		
		if(cpos != viterbi->tbtrees[0]->cpoint->pos || 
		   cstate != viterbi->tbtrees[0]->cpoint->state){
			cpos = viterbi->tbtrees[0]->cpoint->pos;
			cstate = viterbi->tbtrees[0]->cpoint->state;
			fprintf(stderr,"CPOINT at %d (%d) from %d\n",cpos,cstate,viterbi->pos);
		}

		/* check for collapsing paths */
		if(viterbi->active_count > 1){
			zPinViterbiCollapsePaths(viterbi);
			if(viterbi->active_count == 1  && viterbi->dead_count == viterbi->trellis_count-1){
				pin_pos = viterbi->tbtrees[0]->cpoint->pos;
				pin_state = viterbi->tbtrees[0]->cpoint->state;
				pin_score = viterbi->tbtrees[0]->cpoint->score;
				end = clock();
				cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
				fprintf(stderr,"FOUND PIN %d %d after %f cpu seconds\n",
						pin_pos,viterbi->tbtrees[0]->cpoint->state,cpu_time_used);
			}
		}

		/* if we are down to 1 path and it has a cpoint after end_pos */
		if(viterbi->active_count == 1){
			if(cpos > end_pos){
				fprintf(stderr,"STOP FORWARD cpos = %d > %d = end_pos\n",cpos,end_pos);
				break;
			}
		}
		
		/* STEP 1 - create new live nodes */
		zPinViterbiCreateLiveNodes(viterbi);

		/* move through the list of non-collapsed trellises (see end of loop) */
		t = 0;
		while(t >= 0){
			if(viterbi->active[t] == 0){
				if(viterbi->start_pos[t] == viterbi->pos){
					/*fprintf(stderr,"Activate %d (%d active) %d\n",
					  t,viterbi->active_count,viterbi->pos);*/
					viterbi->active[t] = 1;
					viterbi->active_count++;
				}
				else{
					t = viterbi->next[t];
					continue;
				}
			}

			/* STEP 2 - step one pos forward */			
			tree = viterbi->tbtrees[t];
			nodes = viterbi->live_nodes[t];
			cells = viterbi->cache[t][0][viterbi->cdna_pos];
			
			for(state = 0;state < hmm->states;state++){
				/* clear out the temp array */
				viterbi->cache[t][1][viterbi->cdna_pos][state] = NULL;
				/* handle INTERNAL state transitions */
				if(hmm->state[state].type == INTERNAL || 
				   hmm->state[state].type == GINTERNAL){
					best_state = -1;
					best_score = MIN_SCORE;
					jumps = hmm->jmap[state];
					for(i = 0; i < jumps->size; i++) {
						prestate = jumps->elem[i];
						if(cells[prestate] == NULL) continue;
						/*length 1 for duration score. it works properly since length 0 
						  only needed for special case of pos=first position in sequence. 
						  see funtion in zTransition */
						score = zScoreInternalState(trellis,state,viterbi->pos,
													(viterbi->pos == PADDING || state != prestate));
						if(score == MIN_SCORE) continue;
						if(hmm->state[state].type == INTERNAL){
							score += zGetTransitionScore(trellis->hmm,prestate,state,trellis->tiso_group);
						}
						else{
							score += zGetIntergenicContinueScore(trellis);
						}
						score += cells[prestate]->score;
						if(score > best_score){
							best_score = score;		
							best_state = prestate;
						}
					}
					if(best_score != MIN_SCORE){
						/* fill in traceback stuff here */
						tbtn = zGetTBTreeNode(tree);
						viterbi->cache[t][1][viterbi->cdna_pos][state] = tbtn;
						tbtn->pos = viterbi->pos;
						tbtn->state = state;
						tbtn->score = best_score;
						tbtn->frame_data = 0;
						zTBTreeSetChild(cells[best_state],tbtn);
						/*printf("%d trace %d: %d -> %d\n",t,viterbi->pos,state,best_state);*/
					}
				}
				/* EXTERNAL states dealt with in when incorporating live nodes below */
				else if(hmm->state[state].type == EXTERNAL){
				}
				else{
					zDie("zRunViterbi only supports INTERNAL and EXTERNAL states");
				}
			}
			/* cleanup unused/unneeded old_cells */
			for(state = 0;state < hmm->states;state++){
				tbtn = cells[state];
				if(tbtn == NULL) continue;
				if(tbtn->children == 0){
					zReleaseDeadTBTreeNode(tree,tbtn);
				}
				else if((tbtn->children == 1) &&
						(tbtn->child->state == tbtn->state)){
					zReleaseRedundantTBTreeNode(tree,tbtn);
				}
			}
			
			/* STEP 3 - move tmp_cells to cells */
			for(state = 0;state < hmm->states;state++){
				cells[state] = viterbi->cache[t][1][viterbi->cdna_pos][state];
			}
			
			/* STEP 4 - remove/incorporate old live nodes */
			zViterbiIncorporateLiveNodes(viterbi,t);			
			
			/* get the next non-collapsed trellis index */
			t = viterbi->next[t];
		}
	}

	if(pin_pos == (coor_t)-1){
		zDie("No Pin Found!\n");
	}
	zResetSFList(viterbi->sfl);
	/* if the pin is inside the decode range */
	temp_sfl = zMalloc(sizeof(zSFList), "zRunBackwardViterbi sfl");
	zInitSFList(temp_sfl);
	if(pin_pos < end_pos){
		tbtn = viterbi->tbtrees[0]->cpoint;
		/* if don't have a cpoint after the end_pos trace from best path*/
		if(end_pos > tbtn->pos){
			cells = viterbi->cache[0][0][viterbi->cdna_pos];
			best_state = -1;
			best_score = MIN_SCORE;
			for(state = 0;state < hmm->states;state++){
				if(cells[state] != NULL && 
					 cells[state]->score > best_score){
					best_state = state;
					best_score = cells[state]->score - pin_score;
				}
			}
			fprintf(stderr,"Traceback from end\n");
			*path_score = cells[best_state]->score - pin_score;
			zTBTree2SFL(temp_sfl,viterbi->tbtrees[0],cells[best_state],
						viterbi->trellis->hmm,viterbi->trellis);
		}
		/* else trace from cpoint */
		else{			
			fprintf(stderr,"Traceback from cpoint at %d\n",tbtn->pos);
			*path_score = tbtn->score - pin_score;
			zTBTree2SFL(temp_sfl,viterbi->tbtrees[0],tbtn,viterbi->trellis->hmm,viterbi->trellis);
		}
		f = zSFListMoveLast(temp_sfl);
		if(f->start < PADDING){
			f->score -= zGetInitProb(viterbi->trellis->hmm,f->state,viterbi->trellis->iiso_group);
			f->score = PADDING;
		}
		/* trim any features coming before pin */
		f = zSFListMoveLast(temp_sfl);
		while(f != NULL && f->start < pin_pos){
			fprintf(stderr,"PREPINKILL %d->%d %d < %d\n",f->start, f->end, f->state, pin_pos);
			zSFListRemoveLast(temp_sfl);
			f = zSFListMoveLast(temp_sfl);
		}
		/* trim any feature coming after end_pos */
		f = zSFListMoveFirst(temp_sfl);
		while(f != NULL && f->start > end_pos){
			fprintf(stderr,"ENDKILL %d->%d %d > %d\n",f->start, f->end, f->state, end_pos);
			zSFListRemoveFirst(temp_sfl);
			f = zSFListMoveFirst(temp_sfl);
		}
	}
	else{
		*path_score = 0;
	}
	
	fprintf(stderr,"forward pass complete, %d features after trimming\n",temp_sfl->size);

	/* trace back from pin to start_pos */
	zPinReleaseTrellis(viterbi,0);
	zRunBackwardViterbi(viterbi,pin_state,pin_pos,start_pos,&pre_path_score);
	*path_score += pre_path_score;
	if(pin_pos > end_pos){
		/* trim any feature coming after end_pos */
		f = zSFListMoveFirst(viterbi->sfl);
		while(f != NULL && f->end < end_pos){
			fprintf(stderr,"ENDKILL %d->%d %d > %d\n",f->start, f->end, f->state, end_pos);
			zSFListRemoveFirst(viterbi->sfl);
			f = zSFListMoveFirst(viterbi->sfl);
		}
	}
	
	f = zSFListMoveFirst(viterbi->sfl);
	while(f != NULL){
		zSFListMoveLast(temp_sfl);
		zSFListInsertNext(temp_sfl,f);
		f = zSFListMoveNext(viterbi->sfl);
	}
	
	zTranslateSFList(temp_sfl, -PADDING);
	sfv = zMalloc(sizeof(zSFVec),"zRunPinViterbi sfv");
	zInitSFVec(sfv,temp_sfl->size);
	
	zSFList2SFVec(temp_sfl,sfv);
	zFreeSFList(temp_sfl);
	zFree(temp_sfl);

	zFreeViterbi(viterbi);

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stderr,"finished after %f cpu seconds\n",cpu_time_used);
	return sfv;
}

