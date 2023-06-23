/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zHMM_State.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_HMM_STATE_C
#define ZOE_HMM_STATE_C

#include "zHMM_State.h"

/* Reads a string, and sets a state's type and transition function */
/* for now, just internal, external, and phase states */
static zHMM_StateType zGetStateType(char* state_type) {
	     if (strcmp(state_type, "Internal") == 0)    return INTERNAL;
	else if (strcmp(state_type, "GInternal") == 0)   return GINTERNAL;
	else if (strcmp(state_type, "External") == 0)    return EXTERNAL;
	else if (strcmp(state_type, "Explicit") == 0)    return EXPLICIT;
	else { 
		zWarn("Unrecognized state type, defaulting to INTERNAL");
		return INTERNAL;
	}
}

typedef const char* str;
static str zGetStateTypeString(zHMM_StateType type) {
	switch (type) {
	case INTERNAL:	        return "Internal";
	case GINTERNAL:         return "GInternal";
	case EXTERNAL:          return "External";
	case EXPLICIT:          return "Explicit";
	default:        return "UNKNOWN_TYPE";
	}
}


int zReadHMM_State (FILE *stream, zHMM_State *state, int iso_groups) {

	/*
	  This function reads a line of the parameter file that
	  describes a state in the HMM and fills in a zHMM_State
	  object with information from the line.  It also associates
	  external states with feature factory initilization functions
	  and converts intitial probabilities to scores.
	 */

	char  state_name[33];
	char  state_type[33];
	char  strand[2];
	char  dur[33];
	char  model[33];
	char  phase[5];

	float temp;

	score_t *init;
	int i;
	
	init = zCalloc((size_t) iso_groups, (size_t) sizeof(score_t), "zReadHMM_State ISO");
	state->init = zCalloc((size_t)iso_groups, (size_t) sizeof(score_t), "zReadHMM_State ISO");

	/* get state information and initial prob of first isochore */
	if (fscanf(stream, "%32s %32s %1s %4s %32s %32s %f", 
			   state_name, state_type, strand, phase, dur, model, &temp) != 7) {
		zWarn("zReadHMM_State header");
		zFree(init);
		return 0;
	}
	init[0] = temp;
	
	/* read in initial probs of other isochores */
	for (i = 1; i < iso_groups; i++){
		if(fscanf(stream, "%f", &temp) != 1){
			zWarn("zReadHMM_State header ISO");
			zFree(init);
			return 0;
		}
		init[i] = temp;
	}
	
	/* convert strings read in from paramter file to internal representations and store them */
	state->type = zGetStateType(state_type);
	if (zStrIdxExists(state_name)) {
		zWarn("duplicate state name: %s", state_name);
		return 0;
	}
	state->name   = zChar2StrIdx(state_name);
	state->strand = zText2Strand(strand);
	state->phase  = zText2Phase(phase);	
	state->duration = zChar2StrIdx(dur);
	state->model    = zChar2StrIdx(model);

	/* Associate feature factory initilization function with each external state */
	state->ffactory = zGetFactory(model);

	/* external states must have a feature factory initilization function defined to compute their scores */
	if (state->ffactory == NULL && state->type == EXTERNAL) {
		zWarn("Cannot find a factory for state \"%s\" using model \"%s\"",
			  state_name, model);
	}

	/* convert initial probabilities to scores */
	for (i = 0; i < iso_groups; i++) {
		if (init[i] > 0) {
			state->init[i]  = zFloat2Score((double)init[i]);
		} else {
			state->init[i] = MIN_SCORE;
		}
	}
	zFree(init);
	
	return 1;
}

void zWriteHMM_State (FILE *stream, const zHMM_State *state, int iso_groups) {
	char strand[8];
	char phase[8];
	int i;

	zStrand2Text(state->strand, strand);
	zPhase2Text (state->phase,  phase);

	(void)fprintf(stream, "%16s %16s %1s %4s %16s %32s %f",   
				  zStrIdx2Char(state->name),
				  zGetStateTypeString(state->type),
				  strand,
				  phase,
				  zStrIdx2Char(state->duration),
				  zStrIdx2Char(state->model),
				  (float)zScore2Float(state->init[0]));
	for (i = 1; i < iso_groups; i++){
		(void) fprintf(stream, " %f", zScore2Float(state->init[i]));
	}
	(void) fprintf(stream, "\n");
}

#endif
