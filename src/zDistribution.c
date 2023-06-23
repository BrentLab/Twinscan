/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zDuration.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_DISTRIBUTION_C
#define ZOE_DISTRIBUTION_C

#include "zDistribution.h"

int zReadDistribution (FILE *stream, zDistribution *dist) {
	char  func_name[33];
	int   start, end, params, k;
	long  prev_pos = -1;
	float param;

	if (fscanf(stream, "%32s %d %d", func_name, &start, &end) != 3) {
		zWarn("fscanf zReadDistribution func");
		return 0;
	}
	if (0 == start) {
		zWarn("Distributions can not start at 0 (start at 1 instead)");
		return 0;
	}
	dist->start = start;  /* Note the implicit conversion to unsigned */
	dist->end   = end;    /* This makes the value -1 "infinite" */
				
	params = -1; /* impossible value */
	if (strcmp(func_name, "DEFINED") == 0) {
		dist->type = DEFINED;
		params = end - start + 1;
	} else if (strcmp(func_name, "GEOMETRIC") == 0) {
		dist->type = GEOMETRIC;
		params = 2;
	} else if (strcmp(func_name, "POISSON") == 0) {
		 dist->type = POISSON;
		 params = 1;
	} else if (strcmp(func_name, "CONSTANT") == 0) {
		dist->type = CONSTANT;
		params = 1;
        } else if (strcmp(func_name, "POLYSCORE") == 0) {
		dist->type = POLYSCORE;
		if (fscanf(stream, "%d", &params) != 1) {
			zWarn("zReadDistribution params");
			return 0;
		}
		params += 3;
	} else {
		zWarn("zReadDistribution unknown func_name (%s)", func_name);
		return 0;
	}
	 					 
	dist->params = params;
	dist->param  = zMalloc(params * sizeof(float), "zReadDistribution param");
	
	for (k = 0; k < params; k++) {
		if (fscanf(stream, "%g", &param) != 1) {
			if ((dist->type == GEOMETRIC) && (k == 1))
			{
				/* If the pre-factor of geometric distribution */
				/* is not specified, assume it is 1 and return */
				/* file read pointer to previous location      */
				param = 1.;   
				if (fseek(stream, prev_pos, SEEK_SET) == -1) {
					zDie("Cannot fseek on stream");
				}
			}
			else
			{
				zWarn("fscanf zReadDistribution params");
				return 0;
			}
		}
		dist->param[k] = param;
		prev_pos = ftell(stream);
	}
	 return 1;
}

void zWriteDistribution (FILE *stream, const zDistribution *dist) {
	int k;

	switch (dist->type) {
		case DEFINED:       (void)fprintf(stream, "\tDEFINED");   break;
		case GEOMETRIC:     (void)fprintf(stream, "\tGEOMETRIC"); break;
		case POISSON:       (void)fprintf(stream, "\tPOISSON");   break;
		case CONSTANT:      (void)fprintf(stream, "\tCONSTANT");  break;
                case POLYSCORE:     (void)fprintf(stream, "\tPOLYSCORE"); break;
	}
	(void)fprintf(stream, " %d %d\t", dist->start, dist->end);
        if (dist->type == POLYSCORE) (void)fprintf(stream, "%d\t", (dist->params-3));
	for (k = 0; k < dist->params; k++) {
		if ((k % 5) == 0) (void)fprintf(stream, "\n\t\t");
		(void)fprintf(stream, "%.0f\t", dist->param[k]);
	}
	(void)fprintf(stream, "\n");
}

void zFreeDistribution(zDistribution *dist) {
	zFree(dist->param);
}

score_t zScoreDistribution(const zDistribution *dist, coor_t pos) {
	switch (dist->type) {
		case GEOMETRIC:
			return zScoreGeometric(((double)dist->param[0]),
                                               ((double)dist->param[1]), 
		                               ((double)pos));
		case POISSON:
			return zScorePoisson((double)dist->param[0], 
								 (double)(pos - dist->start));
		case CONSTANT:
			return dist->param[0];
		case DEFINED:
			if (pos >= dist->start && pos <= dist->end)
				return dist->param[pos - dist->start];
			else zDie("zScoreDistribution out of bounds");
                case POLYSCORE:
                        return zScorePolyScore(dist->param, (dist->params-3), ((double) pos));
		default: zDie("zScoreDistribution unknown type");
	}
	return MIN_SCORE; /* shush compiler warning */
}

#endif
