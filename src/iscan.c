/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/*****************************************************************************\
 iscan.c

A zoe gene prediction test program

\*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <libgen.h>
#include <time.h>
#include "ZOE.h"

/* auxilliary */
void    zDescribeHMM (zHMM*);
void    zDescribeTrellis (zTrellis*);
void    zOutputSFVec(zTrellis*,zSFVec*, char*, char*);
void    zWriteHeaders(zTrellis* trellis, char* full_command_line,char* parameter_file_name, 
					  char* time_string, char* snp_file, score_t score);
void    zWriteGTFTrace(zDNA*, zHMM*, zSFVec*);
void    zScoreModel(zTrellis*, char*, zSFVec*);
void    zReportFeatures(zTrellis*, char*, zSFVec*);
void    zScorePath(zTrellis*, zSFVec*);
zIVec*  zParseInts(char*);
zSFVec* zGetRanges(const char*, coor_t default_length);
zSFVec* zParseRanges(const char *);
void    zDumpTrellis(zTrellis* trellis, coor_t start, coor_t end);

struct range {
	char   strand;
	coor_t begin;
	coor_t end;
};

typedef enum { NONCONS=1, CONS=2, ALIGN=3 } model_t;

extern int optind;

/* output file, currently always stdout */
FILE *outfile;

/* Specific to sparc architecture - Solaris machines have a different definition of ctime_r */
#if defined (__sparc)
	extern char *ctime_r(const time_t *timep, char *buf, int buflen);
#else
	extern char *ctime_r(const time_t *timep, char *buf);
#endif

void print_help(const char *error) {
	if (NULL != error) {
		fprintf(stderr, "%s\n", error);
	}
	fprintf(stdout, "Usage:\n");
	fprintf(stdout, " Gene Prediction:\n");
	fprintf(stdout, "     iscan [options] hmm_file seq_file [-c=conseq_file | -a=align file] [-e=estseq_file] [-s=snp_file] [-p=start,stop]\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "            Parameters:\n");
	fprintf(stdout, "             Required hmm_file contains the model paramters for initial states,\n");
	fprintf(stdout, "             state transitions, length distributions, DNA models, and\n");
	fprintf(stdout, "             conseq/phylogenetic models\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "            DNA sequence:\n");
	fprintf(stdout, "              Required DNA sequence in Fasta format for which gene predictions\n");
	fprintf(stdout, "              will be generated\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "           Options:\n");
	fprintf(stdout, "             -i  output in zoe format\n");
	fprintf(stdout, "             -pe  use parameterized method for est sequence\n");
	fprintf(stdout, "             -o  run memory optimized code (automatically used with snp file\n");
	fprintf(stdout, "                 or pin mode, explicit states not allowed)\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "            Mode:\n");
	fprintf(stdout, "              Optional conseq_file or align_file  allows iscan to be run in\n");
	fprintf(stdout, "              either TWINSCAN mode (by choosing -c=conseq_file) or in NSCAN\n");
	fprintf(stdout, "              mode (by choosing -a=align_file)\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "            EST sequence:\n");
	fprintf(stdout, "              Optional EST sequence file (EST option invoked by presence of\n");
	fprintf(stdout, "              -e=estseq_file) \n");
	fprintf(stdout, "\n");
	fprintf(stdout, "            SNP file:\n");
	fprintf(stdout, "              Optional SNP file (-s=snp_file) will run in SNP mode creating\n");
	fprintf(stdout, "              predictions for multiple sequence variants\n");
	fprintf(stdout, "\n");
	fprintf(stdout, "            Parallel mode:\n");
	fprintf(stdout, "              Using parallel mode (-p=start,stop) will calclate and return\n");
	fprintf(stdout, "              only the portion of the overall gene predictions between start\n");
	fprintf(stdout, "              and stop (both start and stop must be positive integers)\n");
	fprintf(stdout, "\n");
	fprintf(stdout, " Restricted Gene Prediction: \n");
	fprintf(stdout, "        iscan hmm_file seq_file [-c=conseq_file] -t=zoe_trace \n");
	fprintf(stdout, "        iscan hmm_file seq_file [-c=conseq_file] -T=gtf_trace \n");
	fprintf(stdout, "\n");
	fprintf(stdout, "   With the options below, output can be restricted to specific ranges \n");
	fprintf(stdout, "   using the -r option.  For example -r=1:100,200,-1:-50 \n");
	fprintf(stdout, "   will score a model from coordinates 1 to 100 and 200 on the plus \n");
	fprintf(stdout, "   strand, and -1 to -50 on the minus strand.  On the minus strand, \n");
	fprintf(stdout, "   coordinates start from the 3' end, i.e. 1 and -1 are a base pair. \n");
	fprintf(stdout, "\n");
	fprintf(stdout, "   Score a Model:  \n");
	fprintf(stdout, "        iscan hmm_file seq_file [-c=conseq_file] -m=model_name \n");
	fprintf(stdout, "   Score all Features (Exons, PolyA's, etc.):  \n");
	fprintf(stdout, "        iscan hmm_file seq_file [-c=conseq_file] -f=feature_name \n");
}

/* exit gracefully */
void usage(const char * error) {
	print_help(error);
	exit(0);
}

/* exit with ERRNO = -1 */
void dieusage(const char * error) {
	print_help(error);
	exit(-1);
}

/*****************************************************************************\
 Main Program
\*****************************************************************************/

int main (int argc, char *argv[]) {
	FILE*      stream;
	zHMM       hmm;
	zTrellis   trellis;
	int        verbose;

	zSFVec    *sfv;
	zPtrList  *ptrl = NULL;
	zAlleleDecode *ad = NULL;

	model_t    model = NONCONS;  /* flag that determines conservation model */
	char*      est_mode;
	bool       est_para_mode_bool; 
	bool       snp_mode;
	bool       optimized_mode;
	bool       pin_mode;
	bool       standard_mode;

	char*      hmm_file;
	char*      dna_file;
	char*      unmasked_dna_file;
	char*      conseq_file;
	char*      estseq_file;
	char*      align_file;
	char*      snp_file;
	coor_t     pin_start;
	coor_t     pin_stop;
	int        lib_verbosity = 0;
	score_t    score;
	time_t     stop_time;
	char       *parameter_file_name;
	int        i,size;
	char*      full_command_line;

	outfile = stdout;

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
	if (argc < 3) {
		dieusage("incorrect number of parameters\n");
	}
	hmm_file = argv[optind];
	dna_file = argv[optind+1];
	unmasked_dna_file = zOption("u");
	conseq_file = zOption("c");
	estseq_file = zOption("e");
	align_file = zOption("a");
	snp_file = zOption("s");
	
	if(snp_file != NULL){
		snp_mode = true;
	}
	else{
		snp_mode = false;
	}

	if(zOption("p")){
		const char* input = zOption("p");
		int div = -1;
		if(snp_mode){
			zDie("cannot use both SNP mode (-s) and parallel mode (-p)");
		}
		pin_mode = true;
		for(i = 0; i < (int)strlen(input); i++){
			if(!isdigit(input[i])){
				if(input[i] == ',' && div == -1){
					div = i;
				}
				else{
					zDie("Bad range for -p option.  Must have form \"-p=1234,5678\"\n");
				}
			}
		}
		if(div <= 0 || div == (int)strlen(input)-1){
			zDie("Bad range for -p option.  Must have form \"-p=1234,5678\"\n");
		}
		pin_start = atol(input);
		pin_stop = atol(&(input[div+1]));
	}
	else{
		pin_mode = false;
		pin_start = (coor_t)-1;
		pin_stop = (coor_t)-1;
	}

	if(zOption("o")){
		optimized_mode = true;
	}
	else{
		optimized_mode = false;
	}

	if(snp_mode || pin_mode || optimized_mode){
		standard_mode = false;
	}
	else{
		standard_mode = true;
	}

	/* Set conservation model flag */
	if (zOption("c") != NULL && zOption("a") != NULL) {
		zDie ("Error, cannot use both conseq and alignment\n");
	}
	else if (zOption("c") != NULL) {
		model = CONS;
	}
	else if (zOption("a") != NULL) {
		model = ALIGN;
	}
	else {
		model = NONCONS;
	}

	/* Set estseq */
	est_mode = zOption("e");
	/* Set est_para_mode_bool */
	if(zOption("pe")) {
		est_para_mode_bool = true;
	}
	else {
		est_para_mode_bool = false;
	}

	/* Get verbosity */
	verbose = (zOption("d") ? 1 : 0);
	lib_verbosity = (zOption("v") ? 6 : 0);
	zSetVerbosityLevel(lib_verbosity);
	
	/* read in HMM from file */
	if ((stream = fopen(argv[optind], "r")) == NULL) {
		zDie("hmm file error (%s)", argv[1]);
	} 

	parameter_file_name = argv[optind];

	if (!zReadHMM(stream, &hmm, GHMM)) zDie("error reading hmm");
	fclose(stream);
	if (zOption("h")) zDescribeHMM(&hmm);

	/* DNA */
	if ((stream = fopen(dna_file, "r")) == NULL) {
		zDie("dna file error (%s)", dna_file);
	}
	fclose(stream);

	/* Unmasked DNA */
	if(unmasked_dna_file != NULL){
		if ((stream = fopen(unmasked_dna_file, "r")) == NULL) {
			zDie("unmasked dna file error (%s)", unmasked_dna_file);
		}
		fclose(stream);
	}

	/* EST seq */
	if(est_mode != NULL) {
		if ((stream = fopen(zOption("e"), "r")) == NULL) {
			zDie("Est file error (%s)", zOption("e"));
		}
		fclose(stream);
	}
	else {
		if(est_para_mode_bool) zDie("Est sequence needed\n");
	}

	/* if using conseq read in CONSEQ from file */
	if(model == CONS) {
		if ((stream = fopen(conseq_file, "r")) == NULL) {
			zDie("conseq file error (%s)", zOption("c"));
		}
		fclose(stream);
	}
	
	/* if using an alignment, read it in */
	if (model == ALIGN) {
		if ((stream = fopen(zOption("a"), "r")) == NULL) {
			zDie("alignment file error (%s)", zOption("a"));
		}
		fclose(stream);
	}
	

	/* initilize Trellis object */
	zInitTrellis(&trellis, dna_file, &hmm, conseq_file, estseq_file, align_file,
				 est_para_mode_bool, unmasked_dna_file, snp_file);

	/* Score a model */
	if (zOption("m")) {
		if(!standard_mode){
			zWarn("Option -m supersedes -p, -s, and -o opotions");
		}
		if (zOption("r") == NULL) {
			zDie("Invalid or missing value for -r. Use iscan -h for more detailed information.");
		}
		sfv = zGetRanges(zOption("r"), trellis.dna->length);
		zScoreModel(&trellis, zOption("m"), sfv);
		zFreeSFVec(sfv);
		zFree(sfv);
		exit(0);
	}
	
	/* Run a factory */
	if (zOption("f")) {
		if(!standard_mode){
			zWarn("Option -f supersedes -p, -s, and -o opotions");
		}
		if (zOption("r") == NULL) {
			zDie("Invalid or missing value for -r. Use iscan -h for more detailed information.");
		}
		sfv = zGetRanges(zOption("r"), trellis.dna->length);
		zReportFeatures(&trellis, zOption("f"), sfv);
		zFreeSFVec(sfv);
		zFree(sfv);
		exit(0);
	}
	
	/* Trace a path */
	if (zOption("t")) {
		if(!standard_mode){
			zWarn("Option -t supersedes -p, -s, and -o opotions");
		}
		if ((stream = fopen(zOption("t"), "r")) == NULL) {
			zDie("Could not read features from file \"%s\"", zOption("t"));
		}
		sfv = zReadSFVec(stream);
 		zScorePath(&trellis, sfv); 
		exit(0);
	}

	if(snp_mode){
		ptrl = zRunSNPViterbi(&trellis,&score);
		ad = zPtrListMoveFirst(ptrl);
		sfv = ad->sfv;
	}
	else if(pin_mode){
		sfv = zRunPinViterbi(&trellis,&score,pin_start,pin_stop);
	}
	else if(optimized_mode){
		sfv = zRunViterbi(&trellis,&score);
	}
	else{
		sfv = zRunViterbiAndForward(&trellis,&score);
		if (zOption("dump")) {
			zDumpTrellis(&trellis, atoi(zOption("r")), atoi(zOption("P")));
			zFreeSFVec(sfv);
			zFree(sfv);
			exit(0);
		}
	}
	
	time(&stop_time);
	zWriteHeaders(&trellis,full_command_line,parameter_file_name,ctime(&stop_time),snp_file,score);
	zOutputSFVec(&trellis,sfv,basename(dna_file),NULL);

	if(snp_mode){
		/*output extra decodes here */
		zFreeAlleleDecode(ad);
		zFree(ad);
		ad = zPtrListMoveNext(ptrl);
		while(ad != NULL){
			if(!(ad->score_diff < 0.01 && ad->score_diff > -0.01)){
				zOutputSFVec(&trellis,ad->sfv,basename(dna_file),ad->header);
			}
			zFreeAlleleDecode(ad);
			zFree(ad);
			ad = zPtrListMoveNext(ptrl);
		}
		zFreePtrList(ptrl);
		zFree(ptrl);
	}
	else{
		zFreeSFVec(sfv);
		zFree(sfv);
	}
	
	zFreeTrellis(&trellis);
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

void zDescribeTrellis(zTrellis* t) {
	int   i;
	coor_t j;
	zHMM* hmm = t->hmm;

	printf("===== TRELLIS =====\n\n");

	printf("Externals\n");
	for (i = 0; i < hmm->states; i++) {
		if (t->external[i] == NULL) continue;
		printf("\t%s\n", zStrIdx2Char(hmm->state[i].name));
                
		for (j = 0; j < t->dna->length; j++) {
			if (t->external[i][j].size > 0) {
				printf("%u:%d ", j, t->external[i][j].size);
			}
		}
		printf("\n");
	}
	
	printf("Scanners\n");
	for (i = 0; i < hmm->feature_count; i++) {
		if (t->scanner[i] == NULL) continue;
		printf("\t%s\n", zStrIdx2Char(i));
	}
	
	printf("Factories\n");
	for (i = 0; i < hmm->states; i++) {
		if (t->factory[i] == NULL) continue;
		printf("\t%s\n", zStrIdx2Char(hmm->state[i].name));
	}
	
	printf("Cells\n");
	for (i = 0; i < hmm->states; i++) {
		if (t->cell[i] == NULL) continue;
		printf("\t%s\n", zStrIdx2Char(hmm->state[i].name));
	}
}

void zOutputSFVec(zTrellis* trellis, zSFVec* sfv, char* filename, char* header){
	zGTFVec gtfvec;
	
	zFlipSFVec(sfv);
	if(header != NULL){
		fprintf(outfile, "#%s\n", header);
	}
	else{
		fprintf(outfile, "#\n");
	}
	if (zOption("i")) {
		zWriteSFVec(outfile, sfv);
	} else {
		zInitGTFVec(&gtfvec, sfv->size);
		zSFVec2GTFVec(trellis->hmm, sfv, &gtfvec, filename);
		zSortGTFVec(&gtfvec);
		zWriteGTFVec(outfile, &gtfvec);
		zFreeGTFVec(&gtfvec);
		/*zWriteGTFTrace(trellis->dna, trellis->hmm, sfv); */
	}

}

/* Write the commented header */
void zWriteHeaders(zTrellis* trellis, char* full_command_line,char* parameter_file_name, 
				   char* time_string, char* snp_file, score_t score){	
	fprintf(outfile,"# %s\n", full_command_line);
	fprintf(outfile,"# Date: %s", time_string);
#ifdef BUILD
	fprintf(outfile, "# %s\n", BUILD);
#else
	fprintf(outfile, "# %s\n", zGetProgramName());
#endif
	fprintf(outfile, "# Genome Parameters: %s\n", parameter_file_name);
	if (trellis->conseq != NULL) {
		fprintf(outfile, "# Conservation Parameters: %s\n", parameter_file_name);
	}
	fprintf(outfile, "# Target Sequence: %s\n", trellis->dna->def);
	fprintf(outfile, "# Target Sequence Read... %ubp C+G = %.4f%%\n", trellis->dna->length - 2*PADDING, 100*trellis->dna->gc);
	if (trellis->conseq != NULL) {
		fprintf(outfile, "# Conservation Sequence: %s\n", trellis->conseq->def);
	}
	if(snp_file != NULL){
		fprintf(outfile, "# SNP file: %s\n", snp_file);
	}
	fprintf(outfile, "# Completed in %.2f CPU seconds\n",(double)clock()/CLOCKS_PER_SEC); 
	fprintf(outfile, "# Score: %g\n", score);
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

void    zScoreModel(zTrellis* trellis, char* m, zSFVec* vec) {
	int        i;
	coor_t     pos;
	coor_t     coor;
	coor_t     end;
	coor_t     delta;
	zStrIdx    model;
	zScanner  *s;
	zScanner  *r;

	zScanner  *scanner;
	score_t    score;
	zSfeature *f;
	char       score_str[32];

	model = zChar2StrIdx(m);
	if (model >= trellis->hmm->feature_count) {
		zDie("Name \"%s\" not found anywhere in parameter file\n", m);
	}

	s = trellis->scanner[model];
	r = trellis->rscanner[model];
   
	if (NULL == s || NULL == r) {
		zDie("Model \"%s\" not found in parameter file\n", m);
	}

	if (vec) {
		zTranslateSFVec(vec, trellis->padding);
		for (i = 0; i < vec->size; ++i) {
			f = &vec->elem[i];
			delta = (f->start < f->end) ? 1 : -1;
			end   = f->end + delta;
			for (pos = f->start; pos != end; pos += delta) {
				scanner = (f->strand == '-') ? r   : s;
				coor    = (f->strand == '-') ? trellis->dna->length - pos - 1
					: pos;
				score = scanner->score(scanner,coor);
				zScore2Text(score, score_str);
				printf("%c%u: %s\n", f->strand, pos - trellis->padding + 1, score_str);
			}
		}
		zTranslateSFVec(vec, -trellis->padding);
	}			

}


void zReportFeatures(zTrellis* trellis, char* fac_string, zSFVec* vec) {
	int              i;
	int              delta;
	coor_t           pos;
	coor_t           end,cpos,tmp_start;
	zStrIdx          fac_idx;
	zSFVec           sfv;
	zSfeature       *range,*f;
	zSFList          sfl;

	zInitSFList(&sfl);
	zInitSFVec(&sfv,10);
	
	fac_idx = zChar2StrIdx(fac_string);
	if (fac_idx >= trellis->hmm->feature_count) {
		zDie("Name \"%s\" not found anywhere in parameter file\n", fac_string);
	}

	if (vec) {
		for (i = 0; i < vec->size; ++i) {
			range = &vec->elem[i];
			delta = (range->start < range->end) ? 1 : -1;
			end   = range->end + delta;
			for (pos = range->start; pos != end; pos += delta) {
				if (range->strand == '-') {
					cpos = trellis->dna->length - (pos + trellis->padding) - 1;
					if(trellis->rfactory[fac_idx] != NULL){
						zResetSFList(trellis->rexternal[fac_idx]);		
						trellis->rfactory[fac_idx]->create5(trellis->rfactory[fac_idx], cpos, &sfl);
						f = zSFListMoveFirst(&sfl);
						while(f != NULL){
							f->strand = '-';
							tmp_start = f->start;
							f->start  = trellis->dna->length - f->end - 1;
							f->end    = trellis->dna->length - tmp_start - 1;
							f = zSFListMoveNext(&sfl);
						}
					}
				} else {
					/* create external links ending at pos */
					cpos = pos + trellis->padding;
					if(trellis->factory[fac_idx] != NULL){
						zResetSFList(trellis->external[fac_idx]);		
						trellis->factory[fac_idx]->create3(trellis->factory[fac_idx], cpos, &sfl);
						f = zSFListMoveFirst(&sfl);
						while(f != NULL){
							f->strand = '+';	
							f = zSFListMoveNext(&sfl);
						}
					}
				}
				zResetSFVec(&sfv);
				zSFList2SFVec(&sfl,&sfv);
				zResetSFList(&sfl);
				zTranslateSFVec(&sfv, -trellis->padding);
				zWriteSFVec(stdout, &sfv);
				zTranslateSFVec(&sfv, trellis->padding);
			}
		}
	}
	zFreeSFList(&sfl);
	zFreeSFVec(&sfv);
}

void zScorePath(zTrellis* trellis, zSFVec* sfv) {
	int i, stateidx;
	zSfeature *f;
	zSFVec result;

	zInitSFVec(&result,0);

	for (i = 0; i < sfv->size; ++i) {
		f = &sfv->elem[i];
		stateidx = zGetStateIdxFromStrIdx(trellis->hmm, f->name);
		f->state = stateidx;
/* 		printf("%d %s\n", stateidx,  */
/* 			   zStrIdx2Char(trellis->hmm->state[stateidx].name)); */
	}
	zTranslateSFVec(sfv, trellis->padding);
	printf("## Score of path: %g\n", zGetPathScore(trellis, sfv, &result));
	zTranslateSFVec(&result, -trellis->padding);
	zWriteSFVec(stdout, &result);
}


void zDumpTrellis(zTrellis* trellis, coor_t start, coor_t end) {
	coor_t pos;
	int    state;
	int i = 0;

	printf("dump score in (%u, %u)\n", start, end);

	for (pos = start; pos<= end; pos++) {
		printf("pos %u\n", pos);
		if(i == 10){ i = 0;}
		if(i==0) { 
			printf("\n");

			for (state = 0; state < trellis->hmm->states; state++) {
				printf("%10s\t", zStrIdx2Char(trellis->hmm->state[state].name) );
			}
			printf("\n");
		}

		i ++;
		for (state = 0; state < trellis->hmm->states; state++) {
			printf("%10g\t", trellis->cell[state][pos].score);
		}
		printf("\n");
	}

}
