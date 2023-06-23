/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zDuration.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_DURATION_H
#define ZOE_DURATION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "zDistribution.h"
#include "zSfeature.h"
#include "zTools.h"

/******************************************************************************\
 zDuration

A duration is an ordered collection of some number of zDistributions. For
example, it may be some defined values from 0..100, and then an exponential
distribution out to infinity. Right now this is pretty basic and doesn't do the
smart thing of making sure the values of the various distributions sum to 1.
zDurations must be given a name, which must be one of the legal state types.
zDurations automatically read the distributions of which they are composed.

	zDurationGroup d;
	score_t s;
	if (!zReadDurationGroup(stream, &d) error_handler());
	s = zScoreDuration(&d, 20);
	zWriteDuration(stream, &d);
	zFreeDuration(&d);

\******************************************************************************/

struct zDuration {
	coor_t          min;
	coor_t          max;
	int             distributions;
	zDistribution  *distribution;
};
typedef struct zDuration zDuration;

struct zDurationGroup {
        zStrIdx         name;      /* Name of the group - stripped off from individual zDuration structs */
        int             durations; /* Number of durations */
        zDuration      *duration;  /* Array of durations in the isochore group */
        float          *iso_bound; /* Upper bounds of the durations in "duration" */
};
typedef struct zDurationGroup zDurationGroup;

void    zFreeDurationGroup (zDurationGroup*);
void    zFreeDuration (zDuration*);

float   zReadSubDuration          (FILE *stream, zDuration *dm);
int     zReadDurationGroup        (FILE *stream, zDurationGroup *dg);
void    zWriteDurationGroup       (FILE *stream, const zDurationGroup *dm);
score_t zScoreDuration            (const zDuration *dm, coor_t pos); 
score_t zScoreDurationGroup       (const zDurationGroup *dm, coor_t pos, float gc);  /* Score the duration relevant to the gc value */
int     zGetDurationIsochoreGroup (const zDurationGroup *dm, float gc);     /* Get the index of the duration element using gc */

#endif
