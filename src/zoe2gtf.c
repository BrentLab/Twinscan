/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */  
/*****************************************************************************\
zoe2gtf2.c

Convert a zoe list of features to a GTF file;
\*****************************************************************************/

#include <stdio.h>
#include "ZOE.h"

static const char* usage =
    "Usage: \n\
    zoe2gtf hmm_file zoe_file -s=seqname\n\
seqname is a string that will be used for the first column of the GTF file";


void openHMM(zHMM *hmm, const char* hmm_file , 	bool  est_para_mode_bool ) {
	FILE *hmm_stream;

	if ((hmm_stream = fopen(hmm_file, "r")) == NULL)
		zDie("can't open hmm file (%s)", zOption("c"));
	if (!zReadHMM(hmm_stream, hmm, est_para_mode_bool )) zDie("error reading hmm");
	fclose(hmm_stream);
}

int main(int argc, char** argv) {
    zHMM      hmm;
    char* seqname;
    char* hmm_file;
    char* zoe_file;
    FILE* zoe_stream;

    bool est_para_mode_bool = false; 
    zSFVec  *sfv;
    zGTFVec *gtfvec;

    zSetProgramName(argv[0]);
    zParseOptions(&argc, argv);
    seqname = zOption("s");

    if (argc != 3)
        zDie(usage);
    if (NULL == seqname) 
        zDie("must specify a sequence name (use -s=)\n%s", usage);

    hmm_file = argv[1];
    zoe_file = argv[2];
    
	if(zOption("-pe")) {
        est_para_mode_bool = true;
	}

    openHMM(&hmm, hmm_file, est_para_mode_bool);
		
    if ((zoe_stream = fopen(zoe_file, "r")) == NULL) {
        zDie("Couldn't open zoe file %s", zoe_stream);
    }
    sfv = zReadSFVec(zoe_stream);
    fclose(zoe_stream);
    if (sfv->size == 0) zWarn("Did not read any features");
		
    gtfvec = zMalloc(sizeof(zGTFVec), "GTFVec");
    zInitGTFVec(gtfvec, sfv->size);
    zSFVec2GTFVec(&hmm, sfv, gtfvec, seqname);
    zSortGTFVec(gtfvec);
    zWriteGTFVec(stdout, gtfvec);

    zFreeGTFVec(gtfvec);
    zFree(gtfvec);
    zFreeSFVec(sfv);
    zFree(sfv);

    zFreeHMM(&hmm);

    return 0;
}
