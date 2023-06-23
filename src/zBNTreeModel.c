/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zModel.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_BNTREEMODEL_C
#define ZOE_BNTREEMODEL_C

#include "zBNTreeModel.h"

#define CACHE_SIZE 30e6  /* 30 MB cache per node */

int zReadBNTreeModel (FILE *stream, zBNTreeModel *model) {
	int i;
	char type[64];
	char name[64];
	int length, focus, order;
	
	/* parse the model header */
	if (fscanf(stream, "%s %s %d", name, type, &order) != 3) {
		zWarn("zReadBNTreeModel header failed");
		return 0;
	}

	/* set model attributes */
	model->order = order;
	model->name = zMalloc(strlen(name) + 1, "zModel name");
	(void)strcpy(model->name, name);
		
	if (strcmp(type, "BNTREE_ARRAY")  == 0) model->type = BNTREE_ARRAY;
	else if (strcmp(type, "BNTREE")  == 0) model->type = BNTREE;
	else if (strcmp(type, "BNTREE_CDS") == 0) model->type = BNTREE_CDS;
	else {
		zWarn("zReadBNTreeModel sequence model unknown (%s)", type);
		return 0;
	}
		 
	/* parse the model data */
	if (model->type == BNTREE_ARRAY) {
		/* read length and focus */
		if (fscanf(stream, "%d %d", &length, &focus) != 2) {
			zWarn("zReadBNTreeModel could not read length and focus for BNTREE_ARRAY");
			return 0;
		}
		model->length = length;
		model->focus = focus;
		model->submodels = length;
		model->submodel = zMalloc(length * sizeof(zBNTreeModel), "zBNTreeModel allocate BNTREE_ARRAY submodels");

		/* read in all submodels */
		for (i=0; i<model->submodels; i++) {
			zReadBNTreeModel (stream, &(model->submodel[i]));
		}
	}
	else if (model->type == BNTREE_CDS) {
		model->length = 1;
		model->submodels = 3;
		model->submodel = zMalloc(3 * sizeof(zBNTreeModel), "zBNTree Model allocat BNTREE_CDS submodels");

		/* read in all submodels */
		for (i=0; i<model->submodels; i++) {
			zReadBNTreeModel (stream, &(model->submodel[i]));
		}
	}
	else if (model->type == BNTREE) {
		model->length = 1;
		model->submodels = 0;
		model->submodel = NULL;

		/* read in the BNTree starting from header */
		model->bntree = NewBNTreeFromFP (stream, name, order, CACHE_SIZE); 
	}

	return 1;
}

void zFreeBNTreeModel (zBNTreeModel *model) {
	int i;

	for (i = 0; i < model->submodels; i++) zFreeBNTreeModel(&model->submodel[i]);
        if (model->submodel != NULL) zFree(model->submodel);
	zFree(model->name);
	/* free_bntree(model->bntree);      Need to write this! */
}

void zBNTreeModelFlushCache (zBNTreeModel *model) {
	int i;

	if (model->type == BNTREE_ARRAY ||
		model->type == BNTREE_CDS) {
		for (i=0; i<model->submodels; i++)
			zBNTreeModelFlushCache (&model->submodel[i]);
	}
	else {
		free_all_caches (model->bntree);
	}
}
#endif
