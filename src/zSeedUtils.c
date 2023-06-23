#include <assert.h>
#include "zSeedUtils.h"

static void zAddNewBlock(zVec *range, const zHSP *next); 

int zReadHSP(FILE *stream, zHSP* hsp) {
 	coor_t g_start, c_start, g_end, c_end;      /*  pin regions */
	char line[255];

	(void)fgets(line, sizeof(line), stream);
	if (sscanf(line, "(%u, %u) (%u, %u)", &g_start,  &c_start,  &g_end, &c_end) != 4) {
		zDie("Error reading HSP: %s\n", line);
	}
	/* Seed alignment's HSPs should be ungapped blocks */
	if ((c_end - c_start) != (g_end - g_start)) {
		zDie("Seed alignment needs ungapped HSPs: %s", line);
	}
	hsp->g_start = g_start; 
	hsp->g_end   = g_end;
	hsp->c_start = c_start; 
	hsp->c_end   = c_end; 
	return 1;
}

void zWriteHSP(FILE *stream, zHSP* hsp) {
	fprintf(stream, "(%u, %u) (%u, %u)\n", hsp->g_start,  hsp->c_start,  hsp->g_end, hsp->c_end);
	return;
}

int zTranslateHSPGenomic(zHSP* hsp, coor_t offset) {
	hsp->g_start += offset; 
	hsp->g_end   += offset;
	return 1;
}

int zTranslateHSPCDna(zHSP* hsp, coor_t offset) {
	hsp->c_start += offset; 
	hsp->c_end   += offset; 
	return 1;
}

int zCopyHSP(zHSP* from, zHSP* to) {
	to->g_start = from->g_start;
	to->c_start = from->c_start;
	to->g_end   = from->g_end;
	to->c_end   = from->c_end;
	return 1;
}

int zHSPCmp(const zHSP x, const zHSP y) {
	int result = -100;

	if (((result = zCoortCmp((void*)&x.g_start, (void*)&y.g_start)) == 0) && ((result = zCoortCmp((void*)&x.g_end, (void*)&y.g_end)) == 0)) {
		fprintf(stderr, "zHSPCmp for:\n");
		zWriteHSP(stderr, (zHSP*) &x);
		zWriteHSP(stderr, (zHSP*) &y);
		zDie("Merge: Cannot sort. Tell Mani.");
	}
	return result;
}

int zHSPPtrCmp(const void* a, const void* b) {
	return zHSPCmp(*(zHSP*)a, *(zHSP*)b);
}

static void zSortHSPVec (zVec *vec) {
	int i;
	int size = vec->size;
	zHSP* elements = zMalloc(sizeof(zHSP)*size, "zSortVec: elements");
	for (i = 0; i < size; i++) {
		memcpy(elements+i, vec->elem[i], sizeof(zHSP));
		zFree(vec->elem[i]);
	}
	zFreeVec(vec);
	zInitVec(vec, 2);
	qsort(elements, size, sizeof(zHSP), zHSPPtrCmp);
	for (i = 0; i < size; i++) {
		zHSP *hsp = zMalloc(sizeof(zHSP), "zSortHSPVec: hsp");
		memcpy(hsp, elements+i, sizeof(zHSP));
		zPushVec(vec, hsp);
	}
	zFree(elements);
}

/* Returns: number of HSPs (>=0) if there are any, -1 if EOF is reached */
int zReadSeedAlignment (FILE *stream, zSeedAlignment *seed) {
	int    hsp_count = 0;                       /* hsp_count from pin file */

	coor_t boundary_start, boundary_end;        /* genomic boundaries */
	char   defline[4096];                       /* Fasta def line */
	int    i;
 	char   c;
	char   strand[4];

	/* Validate */
	c = fgetc(stream);
	if (c == EOF) return -1;
	if (c != '>') {
		seed = NULL;
		zWarn("zReadPins > not found");
		return -1;
	}
	(void)ungetc(c, stream);

	/* read the def line */
	(void)fgets(defline, sizeof(defline), stream);
	seed->def = zMalloc(strlen(defline) + 1, "zReadPins defline");
	(void)strcpy(seed->def, defline);
	seed->def[strlen(seed->def) -1] = '\0'; /* remove newline */
	
	/* Get a line. This could be optional boundary and mandatory count */
	(void)fgets(defline, sizeof(defline), stream);
	if (sscanf(defline, "genomic_boundary_start=%u genomic_boundary_end=%u strand=%s", &boundary_start, &boundary_end, strand) == 3) {
		seed->strand                 = zText2Strand(strand);
		seed->gb_start = boundary_start;
		seed->gb_end   = boundary_end;
		(void)fgets(defline, sizeof(defline), stream);
	} else {
		fprintf(stdout, "#Missing genomic boundary: Will scan the whole genomic sequence\n");
		seed->strand                 = zText2Strand(".");
		seed->gb_start = 0;
		seed->gb_end   = 0;
	}

	if (sscanf(defline, "count=%d", &hsp_count) != 1) {
		zWarn("Error reading seed alignment count\n");
		return -1;
	}

	/* Empty seed alignment set */

	if (hsp_count == 0) {
		seed->hsps = 0;
		seed->hsp  = NULL;
		return hsp_count;
	}

	seed->hsp = (zHSP*) zMalloc( hsp_count*sizeof(zHSP), "zReadSeedAlignment: hsp");

	for (i = 0; i < hsp_count; i++) {
		zReadHSP(stream, &seed->hsp[i]);
	}
	seed->hsps = hsp_count; 

	if (seed->hsps > 0 && seed->gb_end == 0) {
		seed->gb_start = MAX(seed->hsp[0].g_start, 10000)  - 10000;
		seed->gb_end   =     seed->hsp[seed->hsps-1].g_end + 10000;
	}
	return hsp_count;
}

void zWriteSeedAlignment(FILE *stream, zSeedAlignment *seed) {
	int i;
	if (seed == NULL) {
		fprintf(stream, ">empty seed alignment\n");
		return;
	}
	fprintf(stream, "%s\ngenomic_boundary_start=%u genomic_boundary_end=%u strand=%c\n", seed->def, seed->gb_start, seed->gb_end, (seed->strand == UNDEFINED_STRAND)?'.':seed->strand);
	fprintf(stream, "count=%d\n", seed->hsps);
	for (i = 0; i < seed->hsps; i++) {
		zWriteHSP(stream, &seed->hsp[i]);
	}
	return;
}

void zCopySeedAlignment(zSeedAlignment *seed, zSeedAlignment *copy) {
	int i;
	copy->strand = seed->strand;
	copy->gb_start = seed->gb_start;
	copy->gb_end = seed->gb_end;
	copy->def  = zMalloc(sizeof(char)*(strlen(seed->def)+1), "zCopySeedAlignment: trellis->copy->def");
	strcpy(copy->def, seed->def);
	copy->hsps = seed->hsps;
	copy->hsp  = (zHSP*) zMalloc(copy->hsps*sizeof(zHSP), "zCopySeedAlignment: copy->hsp");
	for (i = 0; i < copy->hsps; i++) {
		zCopyHSP(&seed->hsp[i], &copy->hsp[i]);
	}
}

void zInitSeedAlignment(zSeedAlignment *seed) {
	seed->hsp = NULL;
	seed->def = NULL;
}

void zFreeSeedAlignment(zSeedAlignment *seed) {
	if (seed == NULL) return;
	zFree(seed->hsp);
	zFree(seed->def);
}

int zTranslateSeedAlignmentGenomic(zSeedAlignment *seed, coor_t offset) {
	int i;
	for (i = 0; i < seed->hsps; i++) {
		zTranslateHSPGenomic(&seed->hsp[i], offset);
	}
	seed->gb_start += offset;
	seed->gb_end   += offset;
	return 1;
}

int zTranslateSeedAlignmentCDna(zSeedAlignment *seed, coor_t offset) {
	int i;
	for (i = 0; i < seed->hsps; i++) {
		zTranslateHSPCDna(&seed->hsp[i], offset);
	}
	return 1;
}

int zTranslateSeedAlignment(zSeedAlignment *seed, coor_t offset) {
	return zTranslateSeedAlignmentGenomic(seed, offset) &&
		zTranslateSeedAlignmentCDna(seed, offset);
}

int zGetAlignmentBlock(zSeedAlignment *mem_blocks, coor_t gpos, coor_t cpos) {
	int k;
	for (k = 0; k < mem_blocks->hsps; k++) {
		zHSP* r = &mem_blocks->hsp[k];
		if (gpos <= r->g_end && cpos <= r->c_end) {
			break;
		}
	}
	if (k == mem_blocks->hsps) {
		zDie("Point (%u, %u) lies outside alignment blocks", gpos, cpos);
	}
	return k;
}

int zReadMultipleSeedAlignments(FILE* stream, zVec *seeds) {
	zSeedAlignment *entry = (zSeedAlignment*) zMalloc(sizeof(zSeedAlignment), "zReadMultipleSeedAlignments: entry");
	if (seeds == NULL) {
		zDie("NULL pointer passed in zReadMultipleSeedAlignments");
	}
	while(zReadSeedAlignment(stream, entry) != -1) {
		zPushVec(seeds, (void*)entry);
		entry = (zSeedAlignment*) zMalloc(sizeof(zSeedAlignment), "zReadMultipleSeedAlignments: entry");
	}
	zFree(entry);
	return seeds->size;
}

/*
Prune HSP of the terminal <prune> positions that
are contiguous matches. Input seeds are ungapped. 
If there is not such a region, the hsp is set to
(0, 0) (0, 0)
If the region has length < 2*prune, the HSP is
set at the middle position of that region with
length of HSP = 1

This prune should be at least BLOCK_OVERLAP to 
make sure that there are at least BLOCK_OVERLAP
matching bases from the edge of the HSP in either
direction

*/

int zPruneSeedAlignment(zDNA *genomic, zDNA *cdna, zSeedAlignment *seed, coor_t prune) {
	int i;
	assert(prune >= BLOCK_OVERLAP);

	if (seed->gb_end > genomic->length) {
		seed->gb_end = genomic->length;
	}
	for (i = 0; i < seed->hsps; i++) {
		zHSP* hsp = &seed->hsp[i];
		zHSP* last = NULL;
		if (i > 0) {
			last = &seed->hsp[i-1];
			/* If any HSP crosses in both coordinates over the previous one, merge them */
			if (last->g_end > hsp->g_start && last->c_end > hsp->c_start) {
				last->g_end = hsp->g_end;
				last->c_end = hsp->c_end;
				hsp->g_start = hsp->g_end = hsp->c_start = hsp->c_end = 0;
			}
		}
	}
	for (i = 0; i < seed->hsps; i++) {
		zHSP* hsp = &seed->hsp[i];
		coor_t length = hsp->g_end - hsp->g_start + 1;
		coor_t index = 0;
		coor_t count = 0;

		/* Check 5' end */
		for (index = 0; index < length; index++) {
			if (zGetDNAS5(genomic, hsp->g_start + index) == zGetDNAS5(cdna, hsp->c_start + index)) {
				count++;
			} else {
				count = 0;
			}
			if (count > prune) {
				hsp->g_start += index;
				hsp->c_start += index;
				break;
			}
		}
		/* If we do not have continuous matches, ignore this HSP */
		if (count < prune) {
			hsp->g_start = hsp->g_end = hsp->c_start = hsp->c_end = 0;
			continue;
		}

		/* Check 3' end */
		count = 0;
		for (index = 0; index < length; index++) {
			if (zGetDNAS5(genomic, hsp->g_end - index) == zGetDNAS5(cdna, hsp->c_end - index)) {
				count++;
			} else {
				count = 0;
			}
			if (count > prune) {
				hsp->g_end -= index;
				hsp->c_end -= index;
				break;
			}
		}

		/* If the HSP is too short, the points will cross each other, so get a mid-point */
		if (hsp->g_end < hsp->g_start) {
			hsp->g_start = (hsp->g_start + hsp->g_end)/2;
			hsp->c_start = (hsp->c_start + hsp->c_end)/2;
			hsp->g_end = hsp->g_start;
			hsp->c_end = hsp->c_start;
		}

		/* If any coordinate is very close to the edge, skip that HSP */

		if (MIN(hsp->g_start, hsp->c_start) < BLOCK_OVERLAP + PADDING ||
		    hsp->g_end > genomic->length - BLOCK_OVERLAP - PADDING ||
		    hsp->c_end > cdna->length - BLOCK_OVERLAP - PADDING) {
			hsp->g_start = hsp->g_end = hsp->c_start = hsp->c_end = 0;
		}
	}
	return 1;
}

/*
seed should have been passed through zPruneSeedAlignment by now 
coordinates here are 1 based.
*/
int zSeedAlignment2AlignmentBlocks(zDNA *genomic, zDNA *cdna, zSeedAlignment *seed, zSeedAlignment *blocks) {
	zVec* buffer = (zVec*) zMalloc(sizeof(zVec), "zSeedAlignment2AlignmentBlocks: buffer");
	int i;
	zHSP* block = (zHSP*) zMalloc(sizeof(zHSP), "zSeedAlignment2AlignmentBlocks: block"); /* Will not be freed at the end of function */
	zHSP* hsp;
	if (blocks == NULL) {
		zDie("NULL Pointer passed for blocks");
	}

	/* To Do: Check for compatibility */

	zInitVec(buffer, 2);

	/* Convert the pruned HSPs into blocks */
	block->c_start = PADDING-1;
	block->g_start = MAX(PADDING-1, seed->gb_start);
	assert(block->c_start > 0);
	assert(block->g_start > 0);
	for (i = 0; i < seed->hsps; i++) {
		hsp = &seed->hsp[i];
		if (hsp->g_start == 0) continue;
		block->g_end = hsp->g_start + BLOCK_OVERLAP;
		block->c_end = hsp->c_start + BLOCK_OVERLAP;
		zPushVec(buffer, block);
		block = (zHSP*) zMalloc(sizeof(zHSP), "zSeedAlignment2AlignmentBlocks: block");
		if (hsp->g_start != hsp->g_end) {
			block->g_start = hsp->g_start - BLOCK_OVERLAP + 1;
			block->c_start = hsp->c_start - BLOCK_OVERLAP;
			block->g_end = hsp->g_end + BLOCK_OVERLAP;
			block->c_end = hsp->c_end + BLOCK_OVERLAP;
			zPushVec(buffer, block);
			block = (zHSP*) zMalloc(sizeof(zHSP), "zSeedAlignment2AlignmentBlocks: block");
		}
		block->g_start = hsp->g_end - BLOCK_OVERLAP + 1;
		block->c_start = hsp->c_end - BLOCK_OVERLAP;
	}
	block->g_end = MIN(genomic->length - 1, seed->gb_end);
	block->c_end = cdna->length - 1;
	zPushVec(buffer, block);

	blocks->strand = seed->strand;
	blocks->gb_start = seed->gb_start;
	blocks->gb_end = seed->gb_end;
	blocks->def  = zMalloc(sizeof(char)*(strlen(seed->def)+1), "zInitPairTrellis: trellis->blocks->def");
	strcpy(blocks->def, seed->def);
	blocks->hsps = buffer->size;
	blocks->hsp  = (zHSP*) zMalloc(blocks->hsps*sizeof(zHSP), "zSeedAlignment2Blocks: blocks->hsp");
	for (i = 0; i < buffer->size; i++) {
		zCopyHSP((zHSP*)buffer->elem[i], &blocks->hsp[i]);
		zFree((zHSP*)buffer->elem[i]);
	}
	zFreeVec(buffer);
	zFree(buffer);
	return blocks->hsps;
}

void zAlignmentBlocks2MemoryBlocks(zSeedAlignment *blocks, zSeedAlignment *mem_blocks) {
	int i;
	zHSP *current, *previous;
	zHSP *r = (zHSP*) zMalloc(sizeof(zHSP), "zAllocViterbiVars: r");
	zVec *range = (zVec*) zMalloc(sizeof(zVec), "zAllocViterbiVars: range");
	zVec *copy = (zVec*) zMalloc(sizeof(zVec), "zAllocViterbiVars: copy");

	zInitVec(range, 2);

	current = &blocks->hsp[0];
	zCopyHSP(current, r);
	r->g_start = current->g_start;
	r->g_end   = current->g_end - 2*BLOCK_OVERLAP - PADDING;
	r->c_start = current->c_start;
	r->c_end   = current->c_end;

	for (i = 1; i < blocks->hsps; i++) {
		zAddNewBlock(range, r);
		current  = &blocks->hsp[i];
		previous = &blocks->hsp[i-1];

		r->g_start =  current->g_start - PADDING; /* Non overlapping genomic coordinates */
		r->g_end   = previous->g_end;
		r->c_start = (previous->c_start < PADDING)?previous->c_start:(previous->c_start - PADDING); /* for hsp[0].c_start = 0 */
		r->c_end   =  current->c_end;
		zAddNewBlock(range, r);

		r->g_start = previous->g_end;
		r->g_end   =  current->g_end - 2*BLOCK_OVERLAP - PADDING;
		r->c_start =  current->c_start - PADDING;
		r->c_end   =  current->c_end;
	}

	r->g_end += (2*BLOCK_OVERLAP + PADDING); /* The last one (or the only one if hsps==1) needs to reach the end of sequence */
	zAddNewBlock(range, r);

	/* Clean Up the blocks - it's a heuristic method, so there are cases where two blocks have the same cdna stretch */

	zInitVec(copy, 2);
	for (i = 0; i < range->size; i++) {
		zHSP* current = (zHSP*) range->elem[i];
		zHSP* last = (zHSP*) copy->last;
		if (copy->size == 0) { /* First entry */
			zPushVec(copy, current);
		} else {
			if (last->c_start == current->c_start && last->c_end == current->c_end) { /* Blocks to be merged */
				last->g_end = current->g_end;
				zFree(current);
			} else { /* As is */
				zPushVec(copy, current);
			}
		}
	}

	
	mem_blocks->strand = blocks->strand;
	mem_blocks->gb_start = blocks->gb_start;
	mem_blocks->gb_end = blocks->gb_end;
	mem_blocks->def  = zMalloc(sizeof(char)*(strlen(blocks->def)+1), "zInitPairTrellis: trellis->blocks->def");
	strcpy(mem_blocks->def, blocks->def);
	mem_blocks->hsps = copy->size;
	mem_blocks->hsp  = (zHSP*) zMalloc(sizeof(zHSP)*mem_blocks->hsps, "zCalculateMemoryBlocks: mem_blocks->hsp");
	for (i = 0; i < mem_blocks->hsps; i++) {
		zCopyHSP((zHSP*) copy->elem[i], &mem_blocks->hsp[i]);
		zFree(copy->elem[i]);
	}
	zFreeVec(copy);
	zFree(copy);
	zFreeVec(range);
	zFree(range);
	zFree(r);
}

static void zAddNewBlock(zVec *range, const zHSP *next) {
	zHSP* last;
	zHSP* current = (zHSP*) zMalloc(sizeof(zHSP), "zAddNewBlock: current");
	
/*
printf("In: %u--%u, %u--%u\n", next->g_start, next->g_end, next->c_start, next->c_end);
*/
	if (next->g_start >= next->g_end) {
		zFree(current);
		return;
	}

	memcpy(current, next, sizeof(zHSP));
	if (range->size == 0) {
		zPushVec(range, current);
		return;
	}
	last = range->last;


	if (next->g_start == last->g_end + 1) { /* No problem -- continuous */
		current->c_end = MAX(current->c_end, last->c_end); /* Just in case! */
		zPushVec(range, current);
	} else if (next->g_start > last->g_end + 1) {
		zWriteHSP(stderr, last);
		zDie("Discontiguous ranges. Tell Mani!");
	} else { /* Crossing over to an earlier block */
		int k, fit;
		zHSP *r;

		/* Find the block the start falls into */
		for (k = 0; k < range->size; k++) {
			r = (zHSP*) range->elem[k];
			if (next->g_start <= r->g_end) {
				break;
			}
		}
		if (k == range->size) {
			zDie("Point (%u, %u) cannot be added", next->g_start, next->c_start);
		}
		fit = k;

		/* Expand every block after that */
		for (k = fit + 1; k < range->size; k++) {
			r = (zHSP*) range->elem[k];
			r->c_start = MIN(current->c_start, r->c_start);
			r->c_end   = MAX(current->c_end,   r->c_end);
		}

		/* Split the <fit> block into two */
		r = (zHSP*) range->elem[fit];
		current->g_end = r->g_end;
		r->g_end = current->g_start - 1;
		current->c_start = MIN(current->c_start, r->c_start);
		current->c_end   = MAX(current->c_end,   r->c_end);
		zPushVec(range, current);
		/* sort it */
		zSortHSPVec(range);
		last = range->last;
		current = (zHSP*) zMalloc(sizeof(zHSP), "zAddNewBlock: current()");
		memcpy(current, next, sizeof(zHSP));
		current->g_start = last->g_end + 1;
		zPushVec(range, current);
	}
/*
{int k; for (k = 0; k < range->size; k++) zWriteHSP(stdout, (zHSP*)range->elem[k]); }
*/
}
