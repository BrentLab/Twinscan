/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zDistribution.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_DISTRIBUTION_H
#define ZOE_DISTRIBUTION_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "zMath.h"
#include "zTools.h"

/******************************************************************************\
 zDistribution

A zDistribution is a mathematical function that can be of several defined types
{DEFINED, GEOMETRIC, NORMAL, LINEAR, POISSON, CONSTANT, POLYSCORE}. 
A zDistribution will return a value in score_t given an argument that is a 
non-negative integer. If you would rather have the value in float, you can 
use the zScore2Float function in zMath.

	zDistribution d;
	score_t s;
	if (!zReadDistribution(stream, &d) error_handler();
	s = zScoreDistribution(&d, 20);
	zWriteDistribution(stream, &d);
	zFreeDistribution(&d);

IO Format:

	TYPE START END
		PARAMS

Examples: The following examples can be used for values of x from 0..9.

Geometric with a mean m = 2 and pre-factor N = 0.1,
P(x) = N*(1/m)*(1-(1/m))^(x-1), for regular geometric distrib., set N = 1:

	GEOMETRIC 0 9
		2
                0.1
For N = 1, the second parameter can be omitted

Poisson with a mean of 4:

	POISSON 0 9
		4

Constant with the same score (1) at all positions:

	CONSTANT 0 9
		1

Defined with specific scores for specific positions:

	DEFINED 0 9
		10 20 30 40 50
		60 70 80 90 100

Exponentiated polynomial of N-th degree (PolyScore),

        10*log2(P(x)) = sum_{k=0}^{N} a_k*((x-xc)/xm)^k,

with N = 3, a_0 = 0.1, a_1 = 0.2, a_2 = 0.3, a_3 = 0.4, xc = 100, xm = 50:

        POLYSCORE 0 9   
                3  
                0.1
                0.2
                0.3
                0.4
                100.0
                50.0

\******************************************************************************/

typedef enum {
	DEFINED,
	GEOMETRIC,
	POISSON,
	CONSTANT,
        POLYSCORE
} zDistributionType;

struct zDistribution {
	zDistributionType type;   /* the type of function (see above) */
	coor_t        start;      /* starting coordinate */
	coor_t        end;        /* ending coordinate */
	int           params;     /* number of parameters to function */
	float        *param;      /* parameters of function */
};
typedef struct zDistribution zDistribution;

void    zFreeDistribution (zDistribution*);
int     zReadDistribution (FILE*, zDistribution*);
void    zWriteDistribution (FILE*, const zDistribution*);
score_t zScoreDistribution (const zDistribution*, coor_t);

#endif
