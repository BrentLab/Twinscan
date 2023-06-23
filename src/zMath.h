/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zMath.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_MATH_H
#define ZOE_MATH_H

#include <math.h>
#include "zTools.h"

/******************************************************************************\
 zMath

zMath provides some of the more common math operations. There are conversion
functions for floats and scores as well as some common distributions.
	
	score_t s;
	s = zScoreGeometric(mean, x);
	s = zScorePoisson(mean, x);
	s = zScoreNormal(mean, std_dev, x);
	s = zScoreLinear(slope, intercept, x);

\******************************************************************************/

score_t zFloat2Score (double);
double  zScore2Float (score_t);
score_t zFloatwiseScoreAdd(score_t, score_t);
double  zLog2 (double);
double  zLnFactorial (int);
score_t zScoreGeometric (double, double, double);
score_t zScorePoisson (double, double);
score_t zScorePolyScore(float *, int, double);
void    zDecToBase (int, int, char*);
int     zBaseToDec (int, const char*);

int     zIntMax (int, int);
int     zIntMin (int, int);
coor_t  zCoorMax (coor_t, coor_t);
coor_t  zCoorMin (coor_t, coor_t);
score_t zScoreMax (score_t, score_t);
score_t zScoreMin (score_t, score_t);

extern const int zPOWER[6][11];
#endif
