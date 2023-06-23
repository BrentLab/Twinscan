#include "zTools.h"
#include "zDNA.h"

#define BLOCK_OVERLAP 15

struct zHSP {
	coor_t g_start; 
	coor_t g_end; 
	coor_t c_start;
	coor_t c_end;
};
typedef struct zHSP zHSP; 

void zWriteHSP(FILE *stream, zHSP *hsp);

struct zSeedAlignment {
	int       hsps; 
	coor_t    gb_start; 
	coor_t    gb_end; 
	zHSP     *hsp; 
	char     *def;
	strand_t  strand; 
}; 
typedef struct zSeedAlignment zSeedAlignment; 

void zInitSeedAlignment(zSeedAlignment *seed);
void zFreeSeedAlignment(zSeedAlignment *seed);
int zReadSeedAlignment (FILE *stream, zSeedAlignment *seed); 
void zWriteSeedAlignment(FILE *stream, zSeedAlignment *seed);
void zCopySeedAlignment(zSeedAlignment *seed, zSeedAlignment *copy);
int zReadMultipleSeedAlignments(FILE* stream, zVec *seeds);
int zTranslateSeedAlignment(zSeedAlignment *seed, coor_t offset); 
int zGetAlignmentBlock(zSeedAlignment *mem_blocks, coor_t gpos, coor_t cpos); 
int zSeedAlignment2AlignmentBlocks(zDNA *genomic, zDNA *cdna, zSeedAlignment *seed, zSeedAlignment *blocks); 
void zAlignmentBlocks2MemoryBlocks(zSeedAlignment *blocks, zSeedAlignment *mem_blocks); 
int zPruneSeedAlignment(zDNA *genomic, zDNA *cdna, zSeedAlignment *seed, coor_t prune); 
