/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: nil -*- */  
/******************************************************************************\
zDuration.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_DURATION_C
#define ZOE_DURATION_C

#include "zDuration.h"

float zReadSubDuration (FILE *stream, zDuration *dm) {
	float ubound;
	char  buf[17];
	int distributions;
	int j;

	if (fscanf(stream, "%16s %f %d", buf, &ubound, &distributions) == 3 && strcmp(buf, "ISO") == 0) {  /* ISO 75 2 */
		dm->distributions = distributions;
		dm->distribution = zMalloc(distributions * sizeof(zDistribution),
								   "zReadSubDuration: zDuration header");

		/* distributions */
		dm->min = (coor_t)-1;
		dm->max = 0;
		for (j = 0; j < distributions; j++) {
			if(!zReadDistribution(stream, &dm->distribution[j])) {
				zWarn("zReadDuration distribution");
				return 0;
			}
			if (dm->distribution[j].start < dm->min) {
				dm->min = dm->distribution[j].start;
			}
			if (dm->distribution[j].end > dm->max) {
				dm->max = dm->distribution[j].end;
			}
		}
	}else{
		zWarn("zReadSubDuration: zDuration Header error");
		return (float)(-1);
	}
	return ubound;
}

int zReadDurationGroup (FILE *stream, zDurationGroup *dm) {
	char          duration_name[17];
	int           distributions, i, j;
	char          line[256];
	char          buf[17];

	int           duration_groups = 1;
	float         ubound;

	while(isspace((int)(line[0]=(char)getc(stream))));                                                   /* skip white space */
	if(fgets(line+1, sizeof(line) - 1, stream) == NULL) return 0;                           /* get first line */

	if (sscanf(line, "%16s %16s %d", duration_name, buf, &duration_groups) == 3 && strcmp(buf, "ISO") == 0){    /* <name> ISO <num> */
		/*              printf("Durations: ISO header found\n"); */
		dm->name      = zChar2StrIdx(duration_name); /* name only for zDurationGroup */
		dm->durations = duration_groups;
		dm->iso_bound = zCalloc(duration_groups, sizeof(float), "zReadDurationGroup(1): zCalloc for iso_bound");
		dm->duration  = zCalloc(duration_groups, sizeof(zDuration), "zReadDurationGroup(1): zCalloc for duration");
		for( i = 0; i < dm->durations; i++){
			if( (ubound = zReadSubDuration(stream, &dm->duration[i])) == (float)(-1)){
				zWarn("zReadDurationGroup: zReadSubDuration error");
				return 0;
			}else{
				dm->iso_bound[i] = ubound / 100.0;
			}
		}
	}else if ( sscanf(line, "%16s %d", duration_name, &distributions) == 2) {
        /* compatibility for old model - single isochore */
		/*              printf("Durations: ISO header NOT found\n"); */

		/* hard coding groups = 1 */
		duration_groups = 1;      /* can match ONLY ONCE */
		i = 0;                    /* hardcoding for the ONLY isochore group */

		/* name */
		dm->name      = zChar2StrIdx(duration_name);
		dm->durations = duration_groups;
		dm->iso_bound = zCalloc(duration_groups, sizeof(float), "zReadDurationGroup(2): zCalloc for iso_bound");
		dm->duration  = zCalloc(duration_groups, sizeof(zDuration), "zReadDurationGroup(2): zCalloc for duration");
		dm->iso_bound[i] = 1.0;

		dm->duration[i].distributions = distributions;
		dm->duration[i].distribution = zMalloc(distributions * sizeof(zDistribution),
											   "zReadDurationGroup header");


		/* distributions */
		dm->duration[i].min = (coor_t)-1;
		dm->duration[i].max = 0;
		for (j = 0; j < distributions; j++) {
			if(!zReadDistribution(stream, &dm->duration[i].distribution[j])) {
				zWarn("zReadDuration distribution");
				return 0;
			}
            if (dm->duration[i].distribution[j].start <= dm->duration[i].max){
                zWarn("%s (%d <= %d)", 
                      "Distributions must be in order and nonoverlapping (%d)",
                      dm->duration[i].distribution[j].start,
                      dm->duration[i].max);
                return 0;
            }
			if (dm->duration[i].distribution[j].start < dm->duration[i].min) {
				dm->duration[i].min = dm->duration[i].distribution[j].start;
			}
			if (dm->duration[i].distribution[j].end > dm->duration[i].max) {
				dm->duration[i].max = dm->duration[i].distribution[j].end;
			}
		}
	}else {
		zWarn("zDurationGroup header error");
		return 0;
	}

	return 1;
}

void zWriteDurationGroup (FILE *stream, const zDurationGroup *dm) {
	int i, j;

	(void) fprintf(stream, "%s ISO %d levels\n", zStrIdx2Char(dm->name), dm->durations);
	for (i = 0; i < dm->durations; i++){
		(void) fprintf(stream, "ISO %f %d\n", dm->iso_bound[i], dm->duration[i].distributions);
		for (j = 0; j < dm->duration[i].distributions; j++) {
			zWriteDistribution(stream, &dm->duration[i].distribution[j]);
		}
	}
}

void zFreeDuration(zDuration *dm) {
	int i;
	
	for (i = 0; i < dm->distributions; i++) {
		zFreeDistribution(&dm->distribution[i]);
	}
	zFree(dm->distribution);
}

void zFreeDurationGroup(zDurationGroup *dm) {
	int i;

	for (i = 0; i < dm->durations; i++) {
		zFreeDuration(&dm->duration[i]);
	}
	zFree(dm->duration);  dm->duration  = NULL;
    zFree(dm->iso_bound); dm->iso_bound = NULL;
}

score_t zScoreDurationGroup(const zDurationGroup *dm, coor_t pos, float gc) {
	int iso_group;
	zDuration *duration;

	if (pos < 1) {
		zWarn("zScoreDurationGroup positive, non-zero integers only");
		return 0;
	}

	iso_group = zGetDurationIsochoreGroup(dm, gc);
	if ( iso_group < 0 || iso_group >= dm->durations ) {
		zWarn("zDurationGroup error: group index %d not found in %s", iso_group, zStrIdx2Char(dm->name) );
	}
	duration = &dm->duration[iso_group];

	return zScoreDuration(duration, pos); 
}

score_t zScoreDuration(const zDuration *dm, coor_t pos) {
	int i, found;

	if (pos < dm->min || pos > dm->max) {
		return MIN_SCORE;
	}
        
	found = -1;
	for (i = 0; i < dm->distributions; i++) {
		if (pos >= dm->distribution[i].start && pos <= dm->distribution[i].end) {
			found = i;
			break;
		} else if (pos >= dm->distribution[i].start && dm->distribution[i].end == 0) {
			found = i;
			break;
		} else if (pos <= dm->distribution[i].end && dm->distribution[i].start == 0) {
			found = i;
			break;
		}
	}
        
	if (found == -1) zDie("zScoreDuration out of bounds %d", pos);
	return zScoreDistribution(&dm->distribution[found], pos);
}

int zGetDurationIsochoreGroup(const zDurationGroup *dm, float gc){
	/*static*/ int iso_group=-1; /* This is static since we are using global GC content now      */
	/* If we need to use progressive scanning, this static variable */
	/* should be removed and group should be checked every time     */
	int i;

	if(iso_group == -1){
		iso_group = 0;
		if(dm->durations > 1){
			for(i = 0; i < dm->durations; i++){
				if(gc < dm->iso_bound[i]){
					iso_group=i;
					break;
				}
			}
		}
	}

	return iso_group;
}

#endif
