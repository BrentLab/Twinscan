/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zModel.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_BNTREEMODEL_H
#define ZOE_BNTREEMODEL_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "bntree/bntree.h"
#include "zMath.h"
#include "zTools.h"

/******************************************************************************\
 zBNTreeModel

zBNTreeModel represents a probability distribution for a base or sequence
of bases defined by a Bayesian network with tree structure (i.e., exactly
one parent for each node except the root)

\******************************************************************************/

typedef enum {
	BNTREE,        /* Single network */
	BNTREE_ARRAY,  /* Array of networks */
	BNTREE_CDS     /* Frame-dependent model */
} zBNTreeModelType;

struct zBNTreeModel {
	zBNTreeModelType type;          /* enumerated above */
	char*                 name;        
	coor_t                length;   /* length of the model */
	coor_t                focus;    /* which pos gets which score */
	/* for BNTreeModels, focus should not take into account model order,
	   unlike DNA and Conseq WWAMs */

	int                   order;
	int                   submodels;/* number of sub-models */
	struct zBNTreeModel  *submodel; /* sub-models for BNTREE_ARRAY */
	BNTree               *bntree;   /* the Bayesian network used for scoring */
};
typedef struct zBNTreeModel zBNTreeModel;

void zFreeBNTreeModel (zBNTreeModel*);
int  zReadBNTreeModel (FILE*, zBNTreeModel*);
void zBNTreeModelFlushCache (zBNTreeModel *);

#endif
