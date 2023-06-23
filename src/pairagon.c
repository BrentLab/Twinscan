/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
 iscan.c

A zoe gene prediction test program

\*****************************************************************************/

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>
#include <sys/resource.h>
#include "ZOE.h"
#include "zHardCoding.h"

/* Auxiliary functions and structures, used only inside this file */

void    zDescribeHMM (zHMM*);
void    zWriteGlobalHeaders(FILE* outfile, char* full_command_line,char* parameter_file_name, char* time_string);
void    zWriteLocalHeaders(FILE* stream, zDNA* genomic, zDNA* cdna, int amode, int smode, double time, score_t score);
zIVec*  zParseInts(char*);
zSFVec* zGetRanges(const char*, coor_t default_length);
zSFVec* zParseRanges(const char *);

#define zGetAlignmentModeString(a) ((a==FORWARD)?"forward":"reversed")
#define zGetSpliceModeString(a)    ((a==FORWARD)?"forward":"REVERSED")

extern int optind;                 /* from <unistd.h> */
extern int FORWARD, REVERSE, BOTH; /* 01, 10, 11 from zAlnFeature.c */

/* Specific to sparc architecture - Solaris machines have a different definition of ctime_r */
#if defined (__sparc)
	extern char *ctime_r(const time_t *timep, char *buf, int buflen);
#else
	extern char *ctime_r(const time_t *timep, char *buf);
#endif

void print_help() {
	puts("Pairagon - PairHMM based program for cDNA to genome alignment");

#ifdef BUILD
/* Program Version information, if available */
	puts("");
	printf("%s\n", BUILD);
#endif

	puts("");
	puts("Usage:");
	puts("    pairagon [--alignment_mode={forward|reverse|both}] [--splice_mode={forward|reverse|both}] [--seed=file] [-i] [--nonull] hmm_file cdna_file genomic_file");
	puts("");
	printf("Arguments:\n%s\n%s\n%s\n",
		"    hmm_file          - file containing the pairHMM model specification",
		"    cdna_file         - input cDNA sequence in Fasta format",
		"    genomic_file      - input genomic sequence in Fasta format");

	/* Options */
	puts("");
	printf("Options:\n%s\n%s\n%s\n%s\n%s\n",
		"	--alignment_mode - direction of the cDNA sequence to consider (default:both)",
		"	--splice_mode    - direction of sense in the genomic sequence (default:both)",
		"	 -i              - output the state sequence rather than the est_genome style output (default:false)",
		"	--nonull         - do not use the null model, if present in the parameter file (default:false)",
	        "	--seed=<file>    - use the seed alignments defined in <file> for Stepping Stone algorithm");
	/* Pin file format */
	puts("");
	printf("%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
		"Format of a seed alignment file: ", 
		"",
		">cdna_seq_description", 
		"genomic_boundary_start=<start_position> genomic_boundary_end=<end_position> strand={+|-}",
		"count=<number of seed HSPs>",
		"(hsp1_genomic_start, hsp1_cdna_start) (hsp1_genomic_end, hsp1_cdna_end)",
		"(hsp2_genomic_start, hsp2_cdna_start) (hsp2_genomic_end, hsp2_cdna_end)",
		"...",
		"",
		"Example:",
		">NM_012345",
		"genomic_boundary_start=1000001 genomic_boundary_end=1100000 strand=+",
		"count=3",
		"(1008001,   1) (1008451,  451)",
		"(1052255, 448) (1052557,  750)",
		"(1083431, 761) (1083755, 1085)");
}

/* exit gracefully */
void usage(const char * error) {
        if (NULL != error) {
                fprintf(stderr, "Pairagon error: %s\n", error);
        }
	print_help();
	exit(0);
}

/* exit with ERRNO = -1 */
void dieusage(const char * error) {
        if (NULL != error) {
                fprintf(stderr, "Pairagon error: %s\n", error);
        }
	print_help();
        exit(-1);
}

/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char *argv[]) {

	/* File Reading / Writing */

	FILE*           stream;

	/* HMM and Trellis */

	zHMM            hmm;        /* HMM model from the parameter file */
	zPairTrellis    trellis;    /* Trellis to store all the information for decoding */

	/* Program Modes */

	bool            snp_mode;
	bool            optimized_mode;
	bool            standard_mode;

	/* Alignment specific modes */

	int             alignment_mode = BOTH; /* Strand of cDNA to align against + strand of genomic */

	/* Input Files */

	char*           parameter_file_name;   /* Input zHMM parameter file */
	char*           cdna_file;             /* Input cDNA fasta file, can be multiple sequences */
	char*           genomic_file;          /* Input genomic fasta file, MUST BE SINGLE SEQUENCE */
	char*           snp_file;              /* SNP file, doesn't work yet. */

	/* Sequences */

	zDNA*           cdna;      /* cDNA sequence */
	zDNA*           genomic;   /* Genomic sequence */

	/* Helper Objects to store Sequence Information */

	int             cdna_entries;          /* No. of cDNA entries in the input cDNA file */
	zVec*           multi_cdna_vec;        /* zVec that stores all the entries in the input cDNA file. Used when reading in from file only */
	zDNA**          multi_cdna;            /* Array to store the cDNA entries in the zVec above */
	zVec*           multi_seed_vec = NULL; /* Seed alignment for each cDNA entry in the cDNA fasta file */
	zSeedAlignment *seed = NULL;           /* The working Seed alignment. Do not free here, will be freed by zInitPairTrellis */

	/* Output Helper Objects */

	int             verbose;
	int             size;
	char*           full_command_line;     /* To print the first line of output */
	int             lib_verbosity = 0;
	score_t         score;                 /* To print the "Score: yadda yadda" */
	time_t          stop_time;             /* To print the "Date: yadda yadda" */
	struct rusage   ru;                    /* System resource usage - time */
	double          previous_usage, current_usage;

	/* General Iterator */

	int             i;

	/*******************************************************************
	 * Handle the command line
	 *******************************************************************/

	size = 0;
	for(i = 0; i < argc; i++){
		size += strlen(argv[i])+1;
	}
	full_command_line = zMalloc(sizeof(char)*size,"main full_command_line");
	sprintf(full_command_line,"%s",argv[0]);
	size = strlen(argv[0]);
	for(i = 1; i < argc; i++){
		sprintf(&(full_command_line[size])," %s",argv[i]);
		size += strlen(argv[i])+1;
	}

	/* set the program name */
	zSetProgramName(argv[0]);
	zParseOptions(&argc, argv);
	
	if (zOption("h") || zOption("-help")) {
		usage(NULL);
	}
	if (argc < 4) {
		dieusage("incorrect number of parameters\n");
	}
	cdna_file = argv[optind+1];
	genomic_file = argv[optind+2];
	snp_file = zOption("s");
	
	if(snp_file != NULL){
		snp_mode = true;
	}
	else{
		snp_mode = false;
	}

	if(zOption("o")){
		optimized_mode = true;
	}
	else{
		optimized_mode = false;
	}

	if(snp_mode || optimized_mode){
		standard_mode = false;
	}
	else{
		standard_mode = true;
	}

	/* Set alignment mode */

	if (zOption("-alignment_mode") != NULL && strcmp(zOption("-alignment_mode"), "forward") == 0) {
		alignment_mode = FORWARD;
	} else if (zOption("-alignment_mode") != NULL && strcmp(zOption("-alignment_mode"), "reverse") == 0) {
		alignment_mode = REVERSE;
	}

	/* Get verbosity */
	verbose = (zOption("d") ? 1 : 0);
	lib_verbosity = (zOption("v") ? 6 : 0);
	zSetVerbosityLevel(lib_verbosity);
	
	/*******************************************************************
	 * Read the files in
	 *******************************************************************/

	/* read in HMM from file */
	if ((stream = fopen(argv[optind], "r")) == NULL) {
		zDie("hmm file error (%s)", argv[1]);
	} 

	parameter_file_name = argv[optind];

	if (!zReadHMM(stream, &hmm, GPAIRHMM)) zDie("error reading hmm");
	fclose(stream);
	if (zOption("h")) zDescribeHMM(&hmm);

	/* Add the NULL model */
	if (zOption("-nonull") == NULL) {
		if (!zNullifyHMM(&hmm)) {
			zDie("Cannot process NULL model. Try without --nonull");
		}
	}

	/* cDNA */
	if ((stream = fopen(cdna_file, "r")) == NULL) {
		zDie("cdna file error (%s)", cdna_file);
	}
	fclose(stream);

	/* genomic */
	if ((stream = fopen(genomic_file, "r")) == NULL) {
		zDie("genomic file error (%s)", genomic_file);
	}
	fclose(stream);

	cdna    = zMalloc(sizeof(zDNA), "main: cdna");
	genomic = zMalloc(sizeof(zDNA), "main: genomic");
	zInitDNA(cdna);
	zInitDNA(genomic);

	zLoadDNAFromFasta(genomic,genomic_file,snp_file);

	/* Read in a cDNA fasta file with multiple sequences */
	multi_cdna_vec = (zVec*) zMalloc(sizeof(zVec), "main: multi_cdna_vec");
	zInitVec(multi_cdna_vec, 2);
	cdna_entries = zLoadMultiDNAFromMultiFasta(multi_cdna_vec, cdna_file, NULL);

	multi_cdna = (zDNA**) zMalloc(cdna_entries*sizeof(zDNA*), "main: multi_cdna");
	for (i = 0; i < cdna_entries; i++) {
		multi_cdna[i] = (zDNA*) multi_cdna_vec->elem[i];
	}

	/* Read the seed alignments */

	if (zOption("-seed") != NULL) {
		if ((stream = fopen(zOption("-seed"), "r")) == NULL) {
			zDie("seed alignment file error (%s)", zOption("-seed"));
		}
		multi_seed_vec = (zVec*) zMalloc(sizeof(zVec), "main: multi_cdna_vec");
		zInitVec(multi_seed_vec, 2);
		if (zReadMultipleSeedAlignments(stream, multi_seed_vec) != cdna_entries) {
			zDie("Incorrect number of seed alignments");
		}
		fclose(stream);
	}

	/*******************************************************************
	 *  Run the alignment algorithm
	 *******************************************************************/

	/* Align each cDNA file */

	time(&stop_time);
	zWriteGlobalHeaders(stdout,full_command_line,parameter_file_name,ctime(&stop_time));
	getrusage(RUSAGE_SELF,&ru);
	previous_usage = ru.ru_utime.tv_sec + ru.ru_stime.tv_sec;
	for (i = 0; i < cdna_entries; i++) {
		score_t best_score          = MIN_SCORE;
		int     best_alignment_mode = -1;
		int     best_splice_mode    = -1;
		int     smode, amode;
		int j, k;
		zIVec  *avec         = zMalloc(sizeof(zIVec), "main: avec");
		zIVec  *svec         = zMalloc(sizeof(zIVec), "main: svec");
		zAFVec *best_afv     = zMalloc(sizeof(zAFVec), "main: best_afv");
		zDNA   *best_genomic = zMalloc(sizeof(zDNA), "main: best_genomic");
		zDNA   *best_cdna    = zMalloc(sizeof(zDNA), "main: best_cdna");

		if (multi_seed_vec != NULL) {
			seed = (zSeedAlignment*) multi_seed_vec->elem[i];
			if (seed->strand != UNDEFINED_STRAND) {
				/* Once you are here, you have seed alignments, so you dont need the original value of alignment_mode */
				alignment_mode = (seed->strand == '-')?REVERSE:FORWARD;
				fprintf(stdout, "# Seed alignment found in %c strand of cDNA. %s alignment_mode enforced\n", seed->strand, (seed->strand=='-')?"reverse":"forward");
			}
			if (seed->gb_end == 0) {
				zWarn("# Empty seed alignment found. Skipping this cDNA");
				zWriteLocalHeaders(stdout, genomic, multi_cdna[i], FORWARD, FORWARD, 0, (score_t)0);
				zFree(best_afv); zFree(best_genomic); zFree(best_cdna);
				continue;
			}
		} else {
			seed = NULL;
		}

		zInitIVec(avec, 1);
		if (alignment_mode == FORWARD) {
			zPushIVec(avec, FORWARD);
		} else if (alignment_mode == REVERSE) {
			zPushIVec(avec, REVERSE);
		} else {
			zPushIVec(avec, FORWARD);
			zPushIVec(avec, REVERSE);
		}

		zInitAFVec(best_afv, 2);
		zInitDNA(best_genomic);
		zInitDNA(best_cdna);

		for (j = 0; j < avec->size; j++) { /* Over all alignment modes */
			amode = avec->elem[j];
			zInitIVec(svec, 1);
			if (zOption("-splice_mode") != NULL) {
				       if (strcmp(zOption("-splice_mode"), "forward") == 0) {
					zPushIVec(svec, FORWARD);
				} else if (strcmp(zOption("-splice_mode"), "reverse") == 0) {
					zPushIVec(svec, REVERSE);
				} else if (strcmp(zOption("-splice_mode"), "both") == 0) {
					zPushIVec(svec, FORWARD);
					zPushIVec(svec, REVERSE);
				} else if (strcmp(zOption("-splice_mode"), "cdna") == 0) {
					zPushIVec(svec, amode);
				}
			} else { /* This is default */
				zPushIVec(svec, amode);
			}
			for (k = 0; k < svec->size; k++) { /* Over all splice modes */
				smode = svec->elem[k];

				fprintf(stderr, "# Running alignment_mode=%s, splice_mode=%s\n", zGetAlignmentModeString(amode), zGetSpliceModeString(smode));
				zCopyDNA(multi_cdna[i], cdna);
				if (smode == REVERSE) { /* splice_mode = reverse */
					zSetHMMStrand(&hmm, '-');
				} else {
					zSetHMMStrand(&hmm, '+');
				}
				if (amode == REVERSE) { /* alignment_mode = reverse */
					zAntiDNA(cdna);
				}
				/* Run Pairagon using the options */
				{
					zAFVec    *afv;
					zInitPairTrellis(&trellis, seed, genomic, cdna, &hmm);

					if(optimized_mode){
						afv = zRunPairViterbi(&trellis,&score);
					}
					else
					{
						afv = zRunPairViterbiAndForward(&trellis,&score);
					}
					
					if (score > best_score) {
						int m;
						best_score = score;
						best_splice_mode = smode;
						best_alignment_mode = amode;
						zFreeDNA(best_genomic);
						zInitDNA(best_genomic);
						zCopyDNA(trellis.genomic, best_genomic);
						zFreeDNA(best_cdna);
						zInitDNA(best_cdna);
						zCopyDNA(trellis.cdna, best_cdna);
						zFreeAFVec(best_afv);
						zInitAFVec(best_afv, 2);
						zCopyAFVec(afv, best_afv);
						for (m = 0; m < best_afv->size; m++) {
							best_afv->elem[m].genomic = best_genomic;
							best_afv->elem[m].cdna = best_cdna;
						}
					}
					zFreeAFVec(afv);
					zFree(afv);

					zFreePairTrellis(&trellis);
					zFreeDNA(cdna);
				}
			}
			zFreeIVec(svec);
		}
		zFreeIVec(avec);

		getrusage(RUSAGE_SELF,&ru);
		current_usage = ru.ru_utime.tv_sec + ru.ru_stime.tv_sec;
		zWriteLocalHeaders(stdout, genomic, multi_cdna[i], best_alignment_mode, best_splice_mode, current_usage-previous_usage, best_score);
		previous_usage = current_usage;
		if (zOption("i") != NULL) {
			zWriteAFVec(stdout, best_afv, 0, 0);
		} else {
			zWriteAlignment(stdout, best_afv, 0);
		}
		fflush(stdout);
		zFree(avec);
		zFree(svec);
		zFreeAFVec(best_afv);
		zFree(best_afv);
		zFreeDNA(best_cdna);
		zFree(best_cdna);
		zFreeDNA(best_genomic);
		zFree(best_genomic);
	}

	/* Multiple cDNA and multiple seed alignment stuff */
	for (i = 0; i < cdna_entries; i++) {
		zFreeDNA((zDNA*)multi_cdna_vec->elem[i]);
		zFree(multi_cdna_vec->elem[i]);
	}
	zFreeVec(multi_cdna_vec);
	zFree(multi_cdna_vec);
	zFree(multi_cdna);
	if (multi_seed_vec != NULL) {
		zFreeVec(multi_seed_vec);
		zFree(multi_seed_vec);
	}

	zFreeDNA(genomic);
	zFree(genomic);
	zFree(cdna);
	zFreeHMM(&hmm);
	zFreeOptions();
	zFree(full_command_line);
	zFreeVerbosityGlobalVariable();
	zStringPoolFree(); /* Don't add anything between this line and "return 0;"*/
	return 0;
}

void zDescribeHMM (zHMM *hmm) {
	int i, j, k;
	char strand[2];
	
	printf("\n===== HMM =====\n\n");
	
	printf("state[]\n");
	for (i = 0; i< hmm->states; i++) {
		zStrand2Text(hmm->state[i].strand, strand);
		printf("%i: %s\t%i %i %s %g", i, zStrIdx2Char(hmm->state[i].name), 
									hmm->state[i].type,
									hmm->state[i].phase,
									strand,
									(float)zScore2Float(hmm->state[i].init[0]));
		for (j = 1; j < hmm->iso_states; j++) {
			printf(" %g", (float) zScore2Float(hmm->state[i].init[j]));
		}
		printf("\n");
	}		
	
	printf("tmap\n");
	for (i = 0; i < hmm->states; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (hmm->tmap[i][j][0] == MIN_SCORE) continue;
			printf("\t%s\t%s\t%g", zStrIdx2Char(hmm->state[i].name), 
					zStrIdx2Char(hmm->state[j].name), hmm->tmap[i][j][0]);
			for (k = 1; k < hmm->iso_transitions; k++) {
				printf("\t%g", hmm->tmap[i][j][k]);
			}
			printf("\n");
		}
	}
	
	printf("jmap\n");
	for (i = 0; i < hmm->states; i++) {
		if (hmm->jmap[i] == NULL) continue;
		printf("\t%s:\t", zStrIdx2Char(hmm->state[i].name));
		for (j = 0; j < hmm->jmap[i]->size; j++) {
			printf("%s ", zStrIdx2Char(hmm->state[hmm->jmap[i]->elem[j]].name));
		}
		printf("\n");
	}
	
	printf("initial score\n");
	for (i = 0; i < hmm->states; i++) {
		printf("\t%s\t%g", zStrIdx2Char(hmm->state[i].name), hmm->state[i].init[0]);
		for (j = 1; j < hmm->iso_states; j++) {
			printf("\t%g", hmm->state[i].init[j]);
		}
	}
	
	
	printf("somap\n"); /* !!! should be somap */
	for (i = 0; i < hmm->feature_count; i++) {
		if (hmm->somap[i] == -1) continue;
		printf("\t%s\t%s\n", zStrIdx2Char(i), 
			zStrIdx2Char(hmm->orig_state[hmm->somap[i]].name));
	}
	
	printf("simap\n");
	for (i = 0; i < hmm->orig_states; i++) {
		printf ("\t%s:", zStrIdx2Char(hmm->orig_state[i].name));
		for (j = 0; j < hmm->simap[i].size; j++) {
			printf(" %s", zStrIdx2Char(hmm->state[hmm->simap[i].elem[j]].name));
		}
		printf("\n");
	}

	printf("dmap\n");
	for (i = 0; i < hmm->feature_count; i++) {
		if (hmm->dmap[i] == NULL) continue;
		printf("\t%s\n", zStrIdx2Char(i));
	}
	
	printf("mmap\n");
	for (i = 0; i < hmm->feature_count; i++) {
		if (hmm->mmap[i] == NULL) continue;
		printf("\t%s\n", zStrIdx2Char(i));
	}
}

/* Write the commented header */
/* This writes the header for the total run */
void zWriteGlobalHeaders(FILE* stream, char* full_command_line,char* parameter_file_name, 
				   char* time_string) {
	fprintf(stream,"# %s\n", full_command_line);
	fprintf(stream,"# Date: %s", time_string);
#ifdef BUILD
	fprintf(stream, "# %s\n", BUILD);
#else
	fprintf(stream, "# %s\n", zGetProgramName());
#endif
	fprintf(stream, "# PairHMM Parameters: %s\n", parameter_file_name);
}

/* This writes the header for the alignment of each cDNA sequence */

void zWriteLocalHeaders(FILE* stream, zDNA* genomic, zDNA* cdna, int amode, int smode, double time, score_t score) {
		fprintf(stream, "# Spliced alignment without quality value for cDNA\n");
		fprintf(stream, "# Genomic Sequence: %s, %ubp\n", genomic->def, genomic->length);
		fprintf(stream, "# cDNA    Sequence: %s, %ubp\n", cdna->def, cdna->length);
		fprintf(stream, "# Note Best alignment is between %s est and forward genome, %s splice sites imply %s gene\n", zGetAlignmentModeString(amode), (smode==REVERSE)?"but":"and", zGetSpliceModeString(smode));
		fprintf(stream, "# Completed in %.2f CPU seconds\n", time); 
		fprintf(stream, "# Optimal score: %f\n", score);
}

zIVec* zParseInts(char* ints) {
	char  *num = ints;
	char  *pos = ints;
	char  *err = NULL;
	int    i;
	zIVec *result = NULL;

	if (NULL == ints) return NULL;

	result = zMalloc(sizeof(zIVec), "zParseInts");
	zInitIVec(result, 0);

	while ('\0' != *pos) {
		if (',' == *pos) {
			*pos = '\0';
			err = num;
			i = strtol(num, &err, 10);
			if (0 == i) {
				/* !!! report error! */
				zWarn("the string \"%s\" in an invalid coor., ignoring", num);
			}
			else {
				zPushIVec(result, i);
			}
			num = pos+1;
		}
		++pos;	
	}
	if (pos > num+1) {
		*pos = '\0';
		err = num;
		i = strtol(num, &err, 10);
		if (0 == i) {
			/* !!! report error! */
			zWarn("0 is an impossible position, ignoring");
		}
		else {
			zPushIVec(result, i);
		}
		num = pos+1;
	}		

	return result;
}

zSFVec* zParseRanges(const char* l) {
	char* list;

	char* start; /* start of a range */
	char* div;   /* the range divider, a colon */
	char* end;   /* a comma, or \0, the end of the range */
	char *err1, *err2;

	int start_num, end_num;

	zSFVec* result = NULL;
	zSfeature cur;

	/* allocate the vector for the result */
	result = zMalloc(sizeof(zSFVec), "zParseRanges: result vector");

	zInitSFVec(result, 1);
	if (l == NULL) {
		zFree(result->elem);
		result->elem=NULL; 
		return result; 
	}

	cur.name = -1; cur.state = 0; cur.score = 0;
	cur.lfrag = cur.rfrag = cur.frame = 0;
	cur.group = NULL;

	/* copy the argument into a working copy */
	/* list = zCalloc(strlen(l)+2, sizeof(char), "zParseRanges: copy arg"); */
	list = zMalloc((strlen(l)+2)*sizeof(char), "zParseRanges: copy arg");
	strcpy(list, l); *(list+strlen(l)+1) = '\0';

	start = div = end = list;
	while (*start != '\0') {
		/* locate the end of the range*/
		while (*end != ',' && *end != '\0') ++end;
		*end = '\0';
		/* find the range divider, if there is one */
		while (*div != ':' && *div != '\0') ++div;
		*div = '\0';

		/* parse the numbers */
		cur.strand = '+';
		start_num = strtol(start, &err1, 10);
		if (err1 != div || start == div) {
			zDie("Couldn't read \"%s\" as a number", start);
		}
		if (end == div || end == div+1) {
			end_num = start_num;
		} else {
			end_num = strtol(div+1, &err2, 10);
			if (err2 != end) {
				zDie("Couldn't read \"%s\" as an upperlimit", div+1);
			}
		}
		if (start_num == 0 || end_num == 0) {
			zDie("Zero is an illegal coordinate");
		}

		/* fill out the struct */
		if (start_num < 0) {
			cur.start = -start_num - 1;
			cur.strand = '-';
		} else {
			cur.start = start_num - 1;
		}
		if (end_num < 0) {
			cur.end = -end_num - 1;
			cur.strand = '-';
		} else {
			cur.end = end_num - 1;
		}
		zPushSFVec (result, &cur);

		start = div = end = end+1;
	}

	zFree(list);
	return result;
}

zSFVec* zGetRanges(const char* line, coor_t default_length) {
	zSFVec* result;
	zSfeature f;
	
	result = zParseRanges(line);
	if (result->size == 0) {
		f.name = -1; 
		f.state = 0; 
		f.start = 0; 
		f.end = default_length - 1;
		f.score = 0;
		f.lfrag = f.rfrag = f.frame = 0;
		f.group = NULL;

		f.strand = '+';
		zPushSFVec(result, &f);
		f.strand = '-';
		zPushSFVec(result, &f);
	}
	return result;
}
