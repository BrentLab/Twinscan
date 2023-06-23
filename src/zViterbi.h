#include "zSfeature.h"
#include "zAlnFeature.h"
#include "zHMM.h"
#include "zTrellis.h"
#include "zPairTrellis.h"

/**********************************************************************\
  zViterbi

  This holds the data needed for the memory optimized viterbi run.
  The same structure is used by both GHMM and GPAIRHMM. zViterbi.c
  implements all the functions that are used in both, and some that are
  used in GHMM exclusively, and zPairViterbi.c implements all the 
  functions that are used in GPAIRHMM exclusively.

  See also:
	zViterbi.c, zPairViterbi.c

\**********************************************************************/

/* zViterbi structure holds the trellis, traceback tree and tbnodes */ 
struct zViterbi{

	/* The trellis pointers: only one of these can be NON-NULL at any time */

	zTrellis*      trellis;
	zPairTrellis*  pair_trellis;

	zHMM*          hmm;

	/* For Evan's SNP stuff */
	short*         snp_count;
	short*         active_snp_count;
	int*           first_snp;
	int*           temp;
	short**        snp_vals;
	int*           snp_int_vals;

	zTBTreeNode**  next_tbtn;
	int*           next_tbti;
	zTBTree**      tbtrees;
	short*         tbt_use_count;
	zPtrList**     live_nodes;

	/* These are added specifically for GPAIRHMM, and both GHMM and GPAIRHMM viterbis have been modified to use the same zViterbi structure */

	zTBTreeNode***** cache;                    /* cache[allele][cache_level][cdna_pos][state] is a zTBTreeNode*. cdna_pos=0 for regular GHMM */
	int            cache_length;               /* How many previous positions should be cached? *
						    * For geometric only, you just need 2, if you   *
						    * have fixed length states then you need that   *
						    * length + 1, if you have explicit states, then *
						    * you need the max_explicit_length + 1          */
	zSFList*       sfl;
	zAFList*       afl;
	int            max_trellis_count;
	short          max_snp_count;
	short          max_active_snp_count;
	float          hap_threshold;
	int            trellis_count;
	zPtrList*      traceback;
	coor_t         pos;
	coor_t         cpos;

	/* These are added for compatibility with PairHMMs */

	coor_t         cdna_pos;
	coor_t         cdna_min;
	coor_t         cdna_max;
	
	int            active_count;
	int            dead_count;
	short*         active;
	coor_t*        start_pos;
	int*           next;

	float***       hapmap;
	short**        hapmax;

};
typedef struct zViterbi zViterbi;

zTBTreeNode* zFindViterbiLiveNode(zViterbi *v,int allele, coor_t pos,int state);
void         zRemoveViterbiLiveNode(zViterbi *v,int t, zTBTreeNode* rem);
int          zViterbiIncorporateLiveNodes(zViterbi* v, int allele);
int          zBackwardViterbiIncorporateLiveNodes(zViterbi* v, int allele);
void         zSNPCopyAllele(zViterbi* v, int a1, int a2);
void         zSNPSwapAlleles(zViterbi* v, int a1,int a2);
void         zSNPReleaseAllele(zViterbi* v, int a);
void         zPinReleaseTrellis(zViterbi* v, int a);
void         zSnpViterbiMergeSFLists(zSFList* base, zSFList* ext);
void         zSNPCollapseAllele(zViterbi* v, int allele);

struct zAlleleDecode{
	char*   header;
	zSFVec* sfv;
	zAFVec* afv;
	zSFList* sfl;
	zAFList* afl;
	coor_t  start;
	coor_t  stop;
	coor_t  diff_start;
	coor_t  diff_stop;
	short*  vals;
	int     snp_count;
	int     first_snp;
	score_t score_diff;
};
typedef struct zAlleleDecode zAlleleDecode;

void         zFreeAlleleDecode(zAlleleDecode*);
void         zSNPFillAlleleHeader(zViterbi* v,zAlleleDecode* ad);
