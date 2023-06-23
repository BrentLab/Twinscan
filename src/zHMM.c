/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zHMM.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_HMM_C
#define ZOE_HMM_C

#include "zHMM.h"
#include "zHardCoding.h"

#define KNOWN_TAGS 9

int zCheckValidTag(const char* tag) {
	int i;
	char tags[KNOWN_TAGS][33] = {   "<STATES>",
					"<STATE_TRANSITIONS>",
					"<STATE_DURATIONS>",
					"<SEQUENCE_MODELS>",
					"<CONSEQ_MODELS>",
					"<ESTSEQ_MODELS>",
					"<PHYLOGENETIC_MODELS>",
					"<STATE_INCREMENTS>",
					"<GTF_CONVERSION>"
				    };
	for (i = 0; i < KNOWN_TAGS; i++) {
		if (strcmp(tags[i], tag) == 0) {
			return 1;
		}
	}
	return 0;
}

static int zReadTag(FILE *stream, char *tag) {
	char line[5000];

	/* Skip over comment lines before tag */
	do {
		int read = fscanf(stream, "%32s", tag);
		if (read == EOF) {
			return 0;
		} else if (read != 1) {
			zDie("tag_error %s", tag);
			return 0;
		}
		if (tag[0] == '#') fgets(line, sizeof(line), stream);
		
	} while (tag[0] == '#');

	if (zCheckValidTag(tag)) {
		return 1;
	} else {
		zDie("Unknown tag: %s", tag);
		return 0;
	}
}

static int zAbortHMM (zHMM *hmm, const char *task, int num) {
	zWarn("zReadHMM aborting (%s, %d, %s)", task, num, hmm->name);
	return 0;
}

int zMapHMM (zHMM *hmm) {

	/*
	  The internal representation of the HMM has a lot of "maps".  These are
	  used to associate the names of things (states, sequence models, 
	  transistions, etc...) represented as string indxes to the to the things
	  themselves.  This function sets up these associations.
	 */
	
	int              i, j, k, m, from, to;
	zStrIdx          name;
	score_t         *score;
	zHMM_State      *ostate, *state;
	zDurationGroup  *dur;

	ostate = hmm->orig_state;
	state  = hmm->state;
	
	/* somap - state original map.  Maps string index (name of the state) to a unique number identifying
	 a state in the parameter file representation of the HMM */
	/* reverse_somap - reverse state original map.  Maps unique number identifying a state in the parameter
	 file representation of the HMM to a string index (name of the state) */
	hmm->somap = zMalloc(hmm->feature_count * sizeof(int), "zHMM: somap");
	hmm->reverse_somap = zMalloc(hmm->feature_count * sizeof(zStrIdx), "zHMM: reverse_somap");

	for (i = 0; i < hmm->feature_count; i++)  hmm->somap[i] = -1; /* initilization */

	for (i = 0; i < hmm->orig_states; i++) {
		hmm->somap[ostate[i].name] = i;
		hmm->reverse_somap[i] = ostate[i].name;
	}
	
	/* dmap - duration map.  Maps the name of the duration as a string index to the actual duration. */
	hmm->dmap = zCalloc(hmm->feature_count, sizeof(zDurationGroup*), "zHMM: dmap");
	for (i = 0; i < hmm->durations; i++) { 
		hmm->dmap[hmm->duration[i].name] = &hmm->duration[i];
	}
	
	/* mmap - model map.  Map the name of the model as a string index to the actual model */
	hmm->mmap = zCalloc(hmm->feature_count, sizeof(zModel*), "zHMM: mmap");
	for (i = 0; i < hmm->models; i++) {
		name = zChar2StrIdx(hmm->model[i].name);
		hmm->mmap[name] = &hmm->model[i];
	}
	
	/* cmmap - conservation model map.  Map the name of the conseq model as a string index to the actual conseq model */
	hmm->cmmap = zCalloc(hmm->feature_count, sizeof(zModel*), "zHMM: cmmap");
	for (i = 0; i < hmm->cons_models; i++) {
		name = zChar2StrIdx(hmm->cons_model[i].name);
		hmm->cmmap[name] = &hmm->cons_model[i];
	}

	/* emmap - est model map.  Map the name of the estseq model as a string index to the actual estseq model */
	if (hmm->est_model != NULL) {
		hmm->emmap = zCalloc(hmm->feature_count, sizeof(zModel*), "zHMM: cmmap");
		for (i = 0; i < hmm->est_models; i++) {
			name = zChar2StrIdx(hmm->est_model[i].name);
			hmm->emmap[name] = &hmm->est_model[i];
		}
	}

	/* ammap - alignment model map.  Map the name of the phylogenetic model as a string index to the actual alignment model */
	hmm->ammap = zCalloc(hmm->feature_count, sizeof(zBNTreeModel*), "zHMM: ammap");
	for (i=0; i<hmm->phylo_models; i++) {
		name = zChar2StrIdx(hmm->phylo_model[i].name);
		hmm->ammap[name] = &hmm->phylo_model[i];
	}

	/* tmap - transitions.  This one's more complex.  Here allocate memory and initilize the transistion matrix. */
	hmm->tmap = zCalloc(hmm->states, sizeof(score_t**), "zHMM tmap level 1");
	for (i = 0; i < hmm->states; i++) {
		hmm->tmap[i] = zMalloc(hmm->states * sizeof(score_t*), "zHMM tmap level 2");
		for (j = 0; j < hmm->states; j++) {
			hmm->tmap[i][j] = zCalloc(hmm->iso_transitions, sizeof(score_t), "zHMM: tmap level 3");
			for(k=0;k<hmm->iso_transitions;k++) {
				hmm->tmap[i][j][k]=MIN_SCORE;
			}
		}
	}

	score = zCalloc(hmm->iso_transitions, sizeof(score_t), "zMapHMM : zTransition score");

	for (i = 0; i < hmm->transitions; i++) { /* process each transistion given in parameter file */
		from  = hmm->somap[hmm->transition[i].from];
		to    = hmm->somap[hmm->transition[i].to];

		/* get transistion score for each isochore */
		for (j = 0; j < hmm->iso_transitions; j++) {
			score[j] = hmm->transition[i].score[j];
		}
		
		/* error checking */
		if (from < 0 || to < 0) {
			zWarn("Unrecognized state in transition #%d.", i);
			return 0;
		}

		/*
		  Since we done blow up states inside the program anymore these
		  simaps should only contain one element.  This just fills in the
		  transistion matrix, its a 3-d matrix because each transistion
		  can have different scores in different isochores.
		 */

		for (j = 0; j < hmm->simap[from].size; j++) {
			for (k = 0; k < hmm->simap[to].size; k++) {
				for(m=0;m<hmm->iso_transitions;m++){
					hmm->tmap[hmm->simap[from].elem[j]]
						[hmm->simap[to].elem[k]][m] = score[m];
				}
			}
		}
	} /* for - process each trans... */	
	
	/* set the intergenic continue scores, different isochores have different mean intergenic lengths */
	hmm->inter_continue = zCalloc(sizeof(score_t), hmm->iso_transitions, "Intergenic continue scores");
	name = zChar2StrIdx("Inter");

	if (name < hmm->feature_count) {
		to = (hmm->simap[hmm->somap[name]]).elem[0];
		dur = hmm->dmap[hmm->state[to].duration];
		for (i = 0; i < hmm->iso_transitions; i++) {
			hmm->inter_continue[i] = zScoreDurationGroup(dur, dur->duration[0].min+2, hmm->iso_transition[i]-0.001) 
				- zScoreDurationGroup(dur, dur->duration[0].min+1, hmm->iso_transition[i]-0.001);
		}
	}

	/* pick up the transitions for internal states */
	for (i = 0; i < hmm->states; i++) {
		dur = hmm->dmap[hmm->state[i].duration];

		/* only deal with twinscan style D-states with geometric distributions here */
		if (hmm->state[i].type != INTERNAL && hmm->state[i].type != GINTERNAL) continue;
		if (dur->duration[0].distribution[0].type != GEOMETRIC) continue;

		for(j = 0; j < hmm->iso_transitions; j++) {
			score[j] = zScoreDurationGroup(dur, dur->duration[0].min+3,hmm->iso_transition[j]-0.001)
				- zScoreDurationGroup(dur, dur->duration[0].min+2,hmm->iso_transition[j]-0.001);
			hmm->tmap[i][i][j] = score[j];
		}
	}
	
	/* jmap - jump map.  Stores all of the states that can transistion into a state.
	 This gets used later to iterate over all the states that can jump into a particular state
	when filling in the trellis. */
	hmm->jmap = zCalloc(hmm->states, sizeof(zIVec*), "zHMM jmap");
	for (i = 0; i < hmm->states; i++) {
		hmm->jmap[i] = zMalloc(sizeof(zIVec), "jmap[]");
		zInitIVec(hmm->jmap[i], 1);
		for (j = 0; j < hmm->states; j++) {
			if (hmm->tmap[j][i][0] != MIN_SCORE) zPushIVec(hmm->jmap[i], j); /* hard coding tmap[][][0] is OK since we
											    are only checking for presence of 
											    transitions */
											 /* check this assumption : 05/21/03 */
		}
	}

	/* fmap - not sure why its called 'fmap'.  Similar to the jmap but stores all states that can be reached
	 from a particular state. */
	hmm->fmap = zCalloc(hmm->states, sizeof(zIVec*), "zHMM fmap");
	for (i = 0; i < hmm->states; i++) {
		hmm->fmap[i] = zMalloc(sizeof(zIVec), "fmap[]");
		zInitIVec(hmm->fmap[i], 1);
		for (j = 0; j < hmm->states; j++) {
			if (hmm->tmap[i][j][0] != MIN_SCORE) zPushIVec(hmm->fmap[i], j); /* see jmap init comments */
		}
	}
	zFree(score);
	return 1;
}

int zReadHMM (FILE *stream, zHMM *hmm, hmm_mode_t hmm_mode) {
	char type[17], name[257], tag[33], buf[17], line[257];
	int i, j, states, transitions, durations, models, cons_models,
		iso_states, iso_transitions, est_models=0, phylo_models = 0;
	float *iso_state, *iso_transition;
	int *somap = NULL;
	int value; /* placeholder for zMapHMM's return value */
	
	hmm->mode = hmm_mode;

	/* Initialize HMM object's pointers */
	hmm->name       = NULL;
	hmm->orig_state = NULL;
	hmm->state      = NULL;
	hmm->transition = NULL;
	hmm->duration   = NULL;
	hmm->model      = NULL;
	hmm->gtf_conv   = NULL;
	hmm->feature_count = 0;
	hmm->dmap  = NULL;
	hmm->mmap  = NULL;
	hmm->somap = NULL;
	hmm->simap = NULL;
	hmm->tmap = NULL;
	hmm->jmap = NULL;
	hmm->fmap = NULL;
	hmm->increments = NULL;

	/* parse HMM header located at top of .zhmm file */
	
	i = fscanf(stream, "%16s name=%256s states=%d transitions=%d durations=%d seq_models=%d conseq_models=%d", type, name, &states,
				&transitions, &durations, &models, &cons_models) ;
	if(i != 7) 
		return zAbortHMM(hmm, "header", 0);


	/* _EST mode? */
	if (fscanf(stream, " estseq_models=%d", &est_models) != 1) {
		est_models = 0;
	}

	/* N-SCAN? */
	if (fscanf(stream, " phylo_models=%d", &phylo_models) != 1) {
		phylo_models = 0;
	}
	
	/* set attributes from header */
	if (strcmp(type, "zHMM") != 0) return zAbortHMM(hmm, "zHMM", 0);
	hmm->orig_states  = states;
	hmm->transitions  = transitions;
	hmm->durations    = durations;
	hmm->models       = models;
	hmm->cons_models  = cons_models;
	hmm->est_models   = est_models;
	hmm->phylo_models = phylo_models;

	hmm->name        = zMalloc(strlen(name) + 1, "zReadHMM name");
	(void)strcpy(hmm->name, name);
	
	/* Allocate memory */
	hmm->orig_state = zMalloc(states * sizeof(zHMM_State), "zReadHMM states");
	hmm->transition = zMalloc(transitions*sizeof(zTransition),"zReadHMM tran");
	hmm->model      = zMalloc(models * sizeof(zModel), "zReadHMM models");
	hmm->cons_model = zMalloc(cons_models * sizeof(zModel), "zReadHMM cons_models");
	hmm->est_model = NULL;
	if(est_models > 0) {
		hmm->est_model = zMalloc(est_models * sizeof(zModel), "zReadHMM est_models");
	}
	hmm->phylo_model = NULL;
	if (phylo_models > 0) {
		hmm->phylo_model = zMalloc(phylo_models * sizeof(zBNTreeModel), "zReadHMM phylo_models");
	}

	hmm->duration   = zMalloc(durations * sizeof(zDurationGroup),"zReadHMM duration");
	hmm->simap      = zMalloc(states * sizeof(zIVec), "zReadHMM simap");

	/* First Parse the states into temp buffer */
	/* This is done independent of the other models since the states need to be set up for everything else to work properly */

	if (!zReadTag(stream, tag)) {
		zDie("Tag error reading STATES");
	}
	while (isspace((int)(line[0] = (char)getc(stream)))); /* skip white space */
	if(fgets(line+1, sizeof(line) - 1, stream) == NULL) return 0;

	if (strcmp(tag, "<STATES>") == 0) {

		/* ISO tag present */
		if(sscanf(line, "%16s %d levels:", buf, &iso_states) == 2 && strcmp(buf,"ISO") == 0 ) { 
			hmm->iso_states = iso_states;
			iso_state = zCalloc(hmm->iso_states, sizeof(float), "zReadHMM iso_states");

			/* read in and store each isochore division */
			for(i=0;i<iso_states;i++) {
				if(fscanf(stream, " %f ", iso_state+i) != 1) {
					zWarn("ISO tag error");
					return 0;
				}
				iso_state[i] /= 100.0;
			} /* for -- read in and store... */

			hmm->iso_state = iso_state;
		} /* if -- ISO tag present */ 
		/* ISO tag absent, set up one big isochore */
		else { 
			hmm->iso_states = 1;
			hmm->iso_state = zCalloc( hmm->iso_states, sizeof(float), "zReadHMM zCalloc for iso_state");
			hmm->iso_state[0] = 1.0;
			for(i=strlen(line)-1;i>=0;i--){
				ungetc(line[i], stream);
			}
		} /* else -- ISO tag absent... */

		/* Read in each HMM state */
		for (i = 0; i < states; i++) {
			if (!zReadHMM_State(stream, &hmm->orig_state[i], hmm->iso_states))

				return zAbortHMM(hmm, "states", i);
		}

		/*
		  NOTE: the simap is an artifact of INTERNAL TRACKING and
		  UNDEFINED STRAND states which are no longer used.  The 
		  idea behind INTERNAL TRACKING states was to make it 
		  possible to define the intron state once in the parameter
		  file and "blow it up" into 6 states inside the program.  
		  This made adding sequence models like the branch point 
		  which interfere with phase checking very difficult.  The
		  decision was made instead to explicily define all phases 
		  of introns in the parameter file.  Similar argument for
		  UNDEFINED STRAND states.  Blowing up states within the
		  program is just a bad idea!  Someday this code will be
		  removed.
		 */

		somap = (int*) zMalloc(100*sizeof(int), "zReadHMM: somap");
		for (i = 0; i < hmm->orig_states; i++) somap[hmm->orig_state[i].name] = i;

		/* initialize the simap */
		for (i = 0; i < hmm->orig_states; i++) zInitIVec(&hmm->simap[i], 2);

		/* figure out how many states will be in the internal HMM representation */
		hmm->states = 0;
		for (i = 0; i < hmm->orig_states; i++) {
			j = 1;
			if (UNDEFINED_STRAND == hmm->orig_state[i].strand) j *= 2;  /* these blow up into 2 states */
			hmm->states += j;
		}

		/* copy over states from parameter file representation to internal representation */
		j = 0; /* j := number of states copied so far */
		hmm->state = zMalloc(hmm->states * sizeof(zHMM_State), "zReadHMM states");

		/* foreach state in the parameter file repersentation... */
		for (i = 0; i < hmm->orig_states; i++) {
			memcpy(&hmm->state[j], &hmm->orig_state[i], sizeof(zHMM_State));
			zPushIVec(&hmm->simap[i], j);
			j++;

			/* ...blow up UNDEFINED STRAND states */
			if (UNDEFINED_STRAND == hmm->state[j-1].strand) {
				memcpy(&hmm->state[j], &hmm->state[j-i], sizeof(zHMM_State));
				zPushIVec(&hmm->simap[i], j);
				hmm->state[j-1].strand = '+';
				hmm->state[j].strand = '-';
				j++;
			}
		}
	} else {
		zDie("Missing Parameters: STATES");
	}

	/* Now we can read the other parameters in whichever order we choose, since they dont depend on each other. */

	while (zReadTag(stream, tag)) {

		if (strcmp(tag, "<STATE_TRANSITIONS>") == 0) {
			
			while (isspace((int)(line[0] = (char)getc(stream)))); /* skip white space */
			if(fgets(line+1, sizeof(line) - 1, stream) == NULL) return 0;

			/* Parse state transitions */

			if (sscanf(line, "%16s %d levels:", buf, &iso_transitions) == 2 && strcmp(buf,"ISO")==0 ) { /* ISO tag present */
				hmm->iso_transitions = iso_transitions;
				iso_transition = zCalloc(hmm->iso_transitions, sizeof(float), "zReadHMM iso_transitions");
				for(i=0;i<iso_transitions;i++){
					if(fscanf(stream, "%f", iso_transition+i) != 1)	{
						zWarn("ISO tag error");
					}
					iso_transition[i] /= 100.0;
				}
				hmm->iso_transition = iso_transition;

			} /* if -- ISO tag present*/
			else { /* ISO tag absent, set up one big isochore */
				hmm->iso_transitions = 1;
				hmm->iso_transition = zCalloc(hmm->iso_transitions, sizeof(float), "zReadHMM zCalloc for iso_transition");
				hmm->iso_transition[0] = 100;
				for(i=strlen(line)-1;i>=0;i--){
					ungetc(line[i], stream);
				}
			} /* else -- ISO tag absent...*/

			/* read in and parse each transistion */
			for (i = 0; i < hmm->transitions; i++) {
				if (!zReadTransition(stream, &hmm->transition[i], hmm->iso_transitions))
					return zAbortHMM(hmm, "transitions", i);
			}
		}
		else if (strcmp(tag, "<STATE_DURATIONS>") == 0) {
			
			/* parse state duration models */

			/* read in and parse each feature duration model */
			for (i = 0; i < hmm->durations; i++) {
				if (!zReadDurationGroup(stream, &hmm->duration[i]))
					return zAbortHMM(hmm, "durations", 1);
			}
		}
		else if (strcmp(tag, "<SEQUENCE_MODELS>") == 0) {
			
			/* parse and ambiguate sequence models */

			for (i = 0; i < models; i++) {

				if (!zReadModel(stream, &hmm->model[i])) return zAbortHMM(hmm, "models", i);
				
				/* set up model to score sequences containing Ns */
				if(hmm->model[i].symbols == 4){
					/* if the sequence model has only 4 symbols convert it to a 5 symbol model */
					zAmbiguateModel(&hmm->model[i], -1); /* the second parameter is really meaningless */
				}
				else if(hmm->model[i].symbols != 5){
					/* if the model does not have 4 or 5 parameters something is wrong */
					zDie("Wrong number of symbols in model %s (should be 4 or 5)",hmm->model[i].name);
				}

				/* add feature to string pool -- if the model is required
				 * by a particular zFeatureFactory but not named by a state
				 * e.g. Coding and Acceptor */
				zChar2StrIdx(hmm->model[i].name); 
			}
		}
		else if (strcmp(tag, "<CONSEQ_MODELS>") == 0) {

			/* parse conservation sequence models */

			if (cons_models > 0) {
				for (i = 0; i < cons_models; i++) {
					if (!zReadModel(stream, &hmm->cons_model[i])) return zAbortHMM(hmm, "cons_models", i);
					/* add feature to string pool -- if the model is required
					 * by a particular zFeatureFactory but not named by a state
					 * e.g. Coding and Acceptor */
					zChar2StrIdx(hmm->cons_model[i].name); 
				}
			}	
		}
		else if (strcmp(tag, "<ESTSEQ_MODELS>") == 0) {

			/* parse est sequence models */

			if (est_models > 0) {
				for (i = 0; i < est_models; i++) {
					if (!zReadModel(stream, &hmm->est_model[i])) return zAbortHMM(hmm, "est_models", i);
				
					/* add feature to string pool -- if the model is required
					 * by a particular zFeatureFactory but not named by a state
					 * e.g. Coding and Acceptor */
					zChar2StrIdx(hmm->est_model[i].name); 
				}
			}
		}
		else if (strcmp(tag, "<PHYLOGENETIC_MODELS>") == 0) {

			/* parse phylogenetic models */

			if (phylo_models > 0) {
				for (i = 0; i < phylo_models; i++) {
					if (!zReadBNTreeModel(stream, &hmm->phylo_model[i])) 
						return zAbortHMM(hmm, "phylo_models", i);
					zChar2StrIdx(hmm->phylo_model[i].name);
				}
			}
		}
		else if (strcmp(tag, "<STATE_INCREMENTS>") == 0) {
			char state_name[33];
			int gincrement, cincrement, state_index;
			hmm->increments = (int**) zMalloc(sizeof(int*)*hmm->orig_states, "zReadHMM: increments");
			for (i = 0; i < hmm->orig_states; i++) {
				if (fscanf(stream, "%s %d %d", state_name, &gincrement, &cincrement) != 3) {
					zDie("Tag error reading %s for line %d from the zhmm file\n", tag, i);
				}
				if (!zStrIdxExists(state_name)) {
					zWarn("zReadHMM from unrecognized state (%s)", state_name);
					return 0;
				}
				state_index = somap[zChar2StrIdx(state_name)];
				hmm->increments[state_index] = (int*) zMalloc(sizeof(int)*2, "zReadHMM: increments[i]");
				hmm->increments[state_index][0] = gincrement;
				hmm->increments[state_index][1] = cincrement;
			}
		}
		else if (strcmp(tag, "<GTF_CONVERSION>") == 0) {

			/* Read GTF Info, used to map internal names to standard GTF features for output */

			hmm->gtf_conv = zMalloc(sizeof(zGTFConvInfo), "zReadHMM GTFConv");
			zReadGTFConvInfo(stream, hmm->gtf_conv);
		}
	}

	/* Freeze the hmm's features */
	hmm->feature_count = zStringPoolCount();

	/* Check for valid parameters */

	if (hmm->transition == NULL) {
		zDie("Missing parameters: STATE TRANSITIONS");
	}
	if (hmm->duration == NULL) {
		zDie("Missing parameters: STATE DURATIONS");
	}
	if (hmm->model == NULL) {
		zDie("Missing parameters: SEQUENCE MODELS");
	}
	if (cons_models > 0 && hmm->cons_model == NULL) {
		zDie("Missing parameters: CONSEQ MODELS");
	}
	if (est_models > 0 && hmm->est_model == NULL) {
		zDie("Missing parameters: ESTSEQ MODELS");
	}
	if (phylo_models > 0 && hmm->phylo_model == NULL) {
		zDie("Missing parameters: PHYLOGENETIC MODELS");
	}
	if (hmm->mode == GPAIRHMM && hmm->increments == NULL) {
		zDie("Missing parameters: STATE INCREMENTS");
	}
	if (hmm->mode == GHMM && hmm->gtf_conv == NULL) {
		zDie("Missing parameters: GTF CONVERSION");
	}

	zFree(somap);
	if ((value = zMapHMM(hmm)) != 1) {
		return value;
	}
	
	/* Fix the transitions for GPAIRHMM to speed Viterbi up */
	if (hmm->mode == GPAIRHMM) {
		zFixInternalTransitions(hmm);
	}
	return value;
}


void zWriteHMM (FILE *stream, const zHMM *hmm) {
	int i;

	(void)fprintf(stream, "zHMM\tname=%s\tstates=%d\ttransitions=%d\tdurations=%d\tseq_models=%d\tconseq_models=%d\n", 
				  hmm->name, hmm->states, hmm->transitions, hmm->durations, 
				  hmm->models, hmm->cons_models);
	(void)fprintf(stream, "\n<STATES>\n");
        (void)fprintf(stream, "ISO %d\n%f", hmm->iso_states, hmm->iso_state[0]);
        for(i = 1; i < hmm->iso_states; i++){
                (void) fprintf(stream, "\t%f", hmm->iso_state[i]);
        }
        (void)fprintf(stream, "\n");
	for (i = 0; i < hmm->orig_states; i++) {
		zWriteHMM_State(stream, &hmm->orig_state[i], hmm->iso_states);
	}
	
	(void)fprintf(stream, "\n<STATE_TRANSITIONS>\n");
        (void)fprintf(stream, "ISO %d\n%f", hmm->iso_transitions, hmm->iso_transition[0]);
        for(i = 1; i < hmm->iso_transitions; i++){
                (void) fprintf(stream, "\t%f", hmm->iso_transition[i]);
        }
        (void)fprintf(stream, "\n");
	for (i = 0; i < hmm->transitions; i++) {
		zWriteTransition(stream, &hmm->transition[i], hmm->iso_transitions);
	}
	
	(void)fprintf(stream, "\n<STATE_DURATIONS>\n");
	for (i = 0; i < hmm->durations; i++) {
		zWriteDurationGroup(stream, &hmm->duration[i]);
	}

	(void)fprintf(stream, "\n<SEQUENCE_MODELS>\n");
	for (i = 0; i < hmm->models; i++) {
		zWriteModel(stream, &hmm->model[i]);
	}

	if (hmm->cons_models > 0) {
		(void)fprintf(stream, "\n<CONSEQ_MODELS>\n");
		for (i = 0; i < hmm->cons_models; i++) {
			zWriteModel(stream, &hmm->model[i]);
		}
	}

	if (hmm->gtf_conv != NULL) {
		(void)fprintf(stream, "\n<GTF_CONVERSION\n");
		zWriteGTFConvInfo(stream, hmm->gtf_conv);
	}
        fputs("\n<STATE_INCREMENTS>\n", stream);
	for (i = 0; i < hmm->states; i++) {
		zStrIdx name = hmm->state[i].name;
		fprintf(stream, "%32s\t%d\t%d\n", zStrIdx2Char(name), zGetGenomicIncrement(hmm, i), zGetCDnaIncrement(hmm, i));
	}

}

void zFreeHMM (zHMM *hmm) {
	int i,j;

	zFree(hmm->name);
	for (i = 0; i < hmm->durations; i++) zFreeDurationGroup(&hmm->duration[i]);
	for (i = 0; i < hmm->models; i++)      zFreeModel(&hmm->model[i]);
	for (i = 0; i < hmm->cons_models; i++) zFreeModel(&hmm->cons_model[i]);
	zFree(hmm->cons_model);

	for (i = 0; i < hmm->est_models; i++) zFreeModel(&hmm->est_model[i]);
	zFree(hmm->est_model);
	
	zFree(hmm->orig_state);
	for (i = 0; i < hmm->states; i++)    zFree(hmm->state[i].init);
	zFree(hmm->state);
	for (i = 0; i < hmm->transitions; i++) {
		zFree(hmm->transition[i].prob);
		zFree(hmm->transition[i].score);
	}
	zFree(hmm->transition);
	zFree(hmm->duration);
	zFree(hmm->model);
	
	if (hmm->gtf_conv != NULL) {
		zFreeGTFConvInfo(hmm->gtf_conv);
		zFree(hmm->gtf_conv);    hmm->gtf_conv = NULL;
	}

	zFree(hmm->iso_state);
	zFree(hmm->iso_transition);
	
	for (i = 0; i < hmm->orig_states; i++) zFreeIVec(&hmm->simap[i]);
	zFree(hmm->simap);
	zFree(hmm->somap);
	zFree(hmm->reverse_somap);
	zFree(hmm->mmap);
	zFree(hmm->dmap);
	
	for (i = 0; i < hmm->states; i++) {
		if (hmm->jmap[i] != NULL) {
			zFreeIVec(hmm->jmap[i]);
			zFree(hmm->jmap[i]);
		}
		if (hmm->fmap[i] != NULL) {
			zFreeIVec(hmm->fmap[i]);
			zFree(hmm->fmap[i]);
		}
		for (j = 0; j < hmm->states; j++) {
			zFree(hmm->tmap[i][j]);
		}
		zFree(hmm->tmap[i]);
		if (hmm->increments != NULL) zFree(hmm->increments[i]);
	}
	if (hmm->increments != NULL) zFree(hmm->increments);
	zFree(hmm->tmap);
	zFree(hmm->jmap);
	zFree(hmm->fmap);
	zFree(hmm->cmmap);
	zFree(hmm->ammap);
	zFree(hmm->inter_continue);
}

int zGetStateIdxFromStrIdx(const zHMM* hmm, zStrIdx stridx) {
	zIVec *vec;
	if (stridx < 0 || stridx >= hmm->feature_count) return -1;

	vec = &hmm->simap[hmm->somap[stridx]];
	if (NULL == vec || vec->size == 0) return -1;

	return vec->elem[0];
}

bool zUsedInExplicitState(const zHMM *hmm, zStrIdx stridx)
{
  /* Figures out whether feature identified by stridx */
  /* corresponds to an EXPLICIT state                 */

  int i;

  for (i = 0; i < hmm->states; i++)
     {
       if ((hmm->state[i].type == EXPLICIT)
          && (hmm->state[i].model == stridx))
         {
           return(1);
         }
     }

  return(0);
}

int zFillInStateIdxFromFeature(const zHMM* hmm, zDNA* dna, 
							   zDNA* rdna, zSfeature *f) {
	int i;
	zPhase_t feature_phase;

	if (f->strand == '-') {
		feature_phase =  zSFBeginPhase(rdna, f);
	} else {
		feature_phase =  zSFEndPhase(dna, f);
	}

	for (i = 0; i < hmm->states; ++i) {
		if (hmm->state[i].strand == f->strand &&
			hmm->state[i].model  == f->name &&
			hmm->state[i].phase  == feature_phase) {
			f->state = i;
			f->name  = hmm->state[i].name;
			return 1;
		}
	}
	return 0;
}

static void zIncludeStopInSfeature(zSfeature *f) {
	if ('-' == f->strand) {
		f->start -=3; 
	} else {
		f->end +=3;
	}
}

void zGTFVec2SFVec(const zHMM *hmm, zDNA *dna, zDNA *rdna,
				   zGTFVec *gtfvec, zSFVec *sfv) {
	int           i;
	zGTFfeature  *gf, *carry;
	const char   *cur_tid;  /* current transcript id */
	zSfeature     sfeature;

	/* This sort is necessary beacuse the conv. algorithm assumes that the 
	 * start and stop codon features bookend a transcript, if they are
	 * annotated.
	 */
	zSortGTFVec(gtfvec);

	zClearSfeature(&sfeature);
	
	carry = NULL;
	cur_tid = "";
	for (i = 0; i < gtfvec->size; i++) {
		gf = &gtfvec->elem[i];
		
		switch (gf->type) {
		case zGTF_START:
		case zGTF_STOP:
			if (strcmp(cur_tid, gf->transcript_id)) {
				/* starting a new transcript*/
				carry = gf;
			} else { /* go back and correct the prev feature */
				if (i < 1) zWarn("zGTF2Zoe: unexpected start_codon");
				else {
					if (sfv->last->name == hmm->gtf_conv->gtft2StrIdx[zGTF_CDS]) {
						sfv->last->name = hmm->gtf_conv->gtft2StrIdx[gf->type];
					} else {
						sfv->last->name = hmm->gtf_conv->gtft2StrIdx[zGTF_ESNGL];
					}
					if (gf->type == zGTF_STOP) 
						zIncludeStopInSfeature(sfv->last);
				}
			}
			break;
		case zGTF_CDS:
			sfeature.start  = gf->start;
			sfeature.end    = gf->end;
			sfeature.strand = gf->strand;
			sfeature.score  = gf->score;
			sfeature.lfrag  = gf->frame;
			sfeature.rfrag  = (gf->end - gf->start + 1 - gf->frame) % 3;
			if (sfeature.strand == '+') {
				sfeature.frame = (sfeature.start + sfeature.lfrag + 2) % 3;
			} else {
				sfeature.frame = (dna->length - sfeature.end 
								  + 1 + sfeature.lfrag) % 3;
			}
			sfeature.group  = NULL;
			
			sfeature.name   = hmm->gtf_conv->gtft2StrIdx[gf->type];
			if (NULL != carry) { /* this really should be Einit or Eterm */
				if (carry->type == zGTF_STOP) 
					zIncludeStopInSfeature(&sfeature);
				sfeature.name = hmm->gtf_conv->gtft2StrIdx[carry->type];
				carry = NULL;
			}
			zPushSFVec(sfv, &sfeature);
			break;
		default:
			break;
		}
		cur_tid = gf->transcript_id;
	}

	for (i = 0; i < sfv->size; i++) {
		if (!zFillInStateIdxFromFeature(hmm, dna,rdna,&sfv->elem[i])) {
		zWarn("Couldn't find appropriate state for %s",
			  zStrIdx2Char(sfv->elem[i].name));
		}
	}
}

static void zFillInGTFStart(zGTFfeature *start, zSfeature *f) {
	start->type   = zGTF_START;
	if (f->strand == '+') {
		start->start = f->start;
		start->end   = start->start +2;
	} else {
		start->end   = f->end;
		start->start = start->end -2;
	}
	start->strand = f->strand;
	start->score  = MIN_SCORE;
	start->frame  = 0;
}

static void zFillInGTFStop(zGTFfeature *stop, zSfeature *f) {
	stop->type   = zGTF_STOP;
	if (f->strand == '+') {
		f->end -= 3;
		stop->start = f->end +1;
		stop->end   = stop->start +2;
	} else {
		f->start +=3;
		stop->end   = f->start -1;
		stop->start = stop->end -2;
	}
	stop->strand = f->strand;
	stop->score  = MIN_SCORE;
	stop->frame  = 0;
}

static void zFillInGTFCDS(zGTFfeature *cds, zSfeature *f) {
	cds->type   = zGTF_CDS;
	cds->start  = f->start;
	cds->end    = f->end;
	cds->strand = f->strand;
	cds->score  = f->score;
	cds->frame  = f->lfrag;
}

static void zFillInGTF5UTR(zGTFfeature *utr, zSfeature *f) {
	utr->type   = zGTF_5UTR;
	utr->start  = f->start;
	utr->end    = f->end;
	utr->strand = f->strand;
	utr->score  = f->score;
	utr->frame  = 0;
}

void zSFVec2GTFVec(const zHMM* hmm, const zSFVec* sfvec, zGTFVec* gtfvec, 
				   char* seqname) {
	int          i, j;
	zStrIdx      model, lookahead_model;
	zSfeature    sf, lookahead_sf;
	zGTFConvInfo *conv = hmm->gtf_conv;
	gtf_t        cur_gtf, lookahead_gtf;
	int          gene_count, transcript_count;
	char        *src;
	const char  *gene_id_fmt;
	const char  *transcript_id_fmt;
	char         gid[64], tid[64];
	zGTFfeature  gfeature;
	int          in_gene = 0;
        int          found_feature = 0;

        lookahead_gtf = 0;
	cur_gtf = zGTF_START;

	gene_count = 0;
	transcript_count = 1;

	src = zGetProgramName();
	gene_id_fmt       = "%1$s.%3$03d";
	transcript_id_fmt = "%1$s.%3$03d.%4$d";

	snprintf(gid, sizeof(gid), gene_id_fmt,
			 seqname, src, gene_count, transcript_count);
	snprintf(tid, sizeof(tid), transcript_id_fmt,
			 seqname, src, gene_count, transcript_count);

	zClearGTFfeature(&gfeature);
	gfeature.seqname = seqname;
	gfeature.source  = src;
	gfeature.gene_id = gid;
	gfeature.transcript_id = tid;

	for (i = 0; i < sfvec->size; ++i) {
		zCopySfeature(&sfvec->elem[i], &sf);
		sf.state = zGetStateIdxFromStrIdx(hmm, sf.name);
		model = hmm->state[sf.state].model;
		cur_gtf = conv->strIdx2gtft[model];

		if (cur_gtf < zGTF_UNDEF) {  /* action needs to happen */
			/* update the gene count */
			if (in_gene == 0) {
				++gene_count;
				snprintf(gid, sizeof(gid), gene_id_fmt,
						 seqname, src, gene_count, transcript_count);
				snprintf(tid, sizeof(tid), transcript_id_fmt,
						 seqname, src, gene_count, transcript_count);
			}

			if (zGTF_STOP == cur_gtf) {
				/* add the stop_codon */
				zFillInGTFStop(&gfeature, &sf);
				zPushGTFVec(gtfvec, &gfeature);

				/* add the CDS */
				zFillInGTFCDS(&gfeature, &sf);
				zPushGTFVec(gtfvec, &gfeature);

				in_gene = (sf.strand == '-');
			} else if (zGTF_CDS == cur_gtf) {
				/* add the CDS */
				zFillInGTFCDS(&gfeature, &sf);
				zPushGTFVec(gtfvec, &gfeature);

				in_gene = 1;
                        } else if ((zGTF_START == cur_gtf)
                                  || (zGTF_ESNGL == cur_gtf)
                                  || (zGTF_5UTR == cur_gtf)) {

                                if (zGTF_5UTR == cur_gtf) {
                                   /* add the 5UTR */
                                   zFillInGTF5UTR(&gfeature, &sf);
                                   zPushGTFVec(gtfvec, &gfeature);
                                }
                                else {
				   /* add the start_codon */
				   zFillInGTFStart(&gfeature, &sf);
				   zPushGTFVec(gtfvec, &gfeature);

                                   if (zGTF_ESNGL == cur_gtf) {
     				       /* add the stop_codon */
				       zFillInGTFStop(&gfeature, &sf);
				       zPushGTFVec(gtfvec, &gfeature);
                                   }

				   /* add the CDS */
				   zFillInGTFCDS(&gfeature, &sf);
				   zPushGTFVec(gtfvec, &gfeature);
                                }

                                if (sf.strand == '+') {
                                        in_gene = (zGTF_ESNGL != cur_gtf);
                                }
                                else {
                                        /* Look ahead--if then next feature is a minus-strand UTR,
                                           we're still in the same gene */
                                        found_feature = 0;
                                        j = i+1;
                                        while ((j < sfvec->size) && (!found_feature)) {
                                                zCopySfeature(&sfvec->elem[j], &lookahead_sf);
                                                lookahead_sf.state = zGetStateIdxFromStrIdx(hmm, lookahead_sf.name);
                                                lookahead_model = hmm->state[lookahead_sf.state].model;
                                                lookahead_gtf = conv->strIdx2gtft[lookahead_model];
                                                if (lookahead_gtf < zGTF_UNDEF) {
                                                        found_feature = 1;
                                                }
                                                j++;
                                        }
                                        if ((lookahead_sf.strand == '-') && (zGTF_5UTR == lookahead_gtf)) {
                                                in_gene = 1;
                                        }
                                        else {
                                                in_gene = 0;
                                        }
                                }

                        }
		}
	}
}

score_t zGetInitProb(zHMM* hmm, int state_index, int iso_group){
	if (iso_group < 0 || iso_group >= hmm->iso_states) {
		zDie("zGetInitProb: iso_group %d not found\n", iso_group);
	}
	return  hmm->state[state_index].init[iso_group];
}

/*********************************************************************************
 * The exit score for the GEOMETRIC state is included in its initial probability *
 * to save tim comparing the length of the GEOMETRIC state later. Saves a lot of *
 * time in PairHMM at least. Used in PairHMM ONLY, as of now.                    *
 *********************************************************************************/

score_t zGetFixedInitProb(zHMM* hmm, int state_index, int iso_group, float gc){
	score_t score;

	if (hmm->mode != GPAIRHMM) {
		zDie("Cannot Nullify initial probabilities for non-GPAIRHMM mode");
	}

	score = zGetInitProb(hmm, state_index, iso_group);
	if (hmm->state[state_index].type == INTERNAL || hmm->state[state_index].type == GINTERNAL) {
		zDurationGroup* group = hmm->dmap[hmm->state[state_index].duration];
		int iso_group = zGetDurationIsochoreGroup(group, gc);
		score += zScoreDuration(&group->duration[iso_group], 1);
	}
	return score;
}

/*********************************************************************************
 * This now includes the exit probability of the to state, which we used to add  *
 * in the first from->to transition. May be I should check for geometric states? *
 * Right now everything else has a 0 score for length 1. I check for internal    *
 * states - inherently assuming they are geometric                               * 
 *********************************************************************************/

int zFixInternalTransitions(zHMM* hmm) {
	int i, j;

	if (hmm->mode != GPAIRHMM) {
		zDie("Cannot fix internal transitions for non-GPAIRHMM mode");
	}		/* get transistion score for each isochore */

	for (i = 0; i < hmm->transitions; i++) { /* process each transistion given in parameter file */
		int from  = hmm->somap[hmm->transition[i].from];
		int to    = hmm->somap[hmm->transition[i].to];
		for (j = 0; j < hmm->iso_transitions; j++) {
			/* Get the mean gc content of the transition isochore and get the duration score corresponding to that */
			float gc = hmm->iso_transition[j]; 
			if (j > 0) { 
				gc += hmm->iso_transition[j-1];
			}
			gc /= 2;
			if (hmm->state[to].type == INTERNAL || hmm->state[to].type == GINTERNAL) { 
				hmm->tmap[from][to][j] += zScoreDurationGroup(hmm->dmap[hmm->state[to].duration], 1, gc);
			}
		}
	}
	
	return 1;
}

/*
  Score in the NULL models
*/

static int zNullifyTransitions (zHMM *hmm, const char *genomic_null_duration_name, const char *cdna_null_duration_name) {
	int i;
	zTransition *t;
	float gc = 50.0; /* It's just a number, so NULL model assumes 50% GC for parameters irrespective of sequence GC percent*/
	zDurationGroup *dur;
	zStrIdx random_genomic = zChar2StrIdx(genomic_null_duration_name);
	zStrIdx random_cdna = zChar2StrIdx(cdna_null_duration_name);
	score_t genomic_null_trans, cdna_null_trans;

	if (hmm->mode != GPAIRHMM) {
		zDie("Cannot Nullify transitions for non-GPAIRHMM mode");
	}

	if (random_genomic >= hmm->feature_count) {
		return 0;
	}
	if (random_cdna >= hmm->feature_count) {
		return 0;
	}

	for (i = 0; i < hmm->transitions; i++) {
		int from_state, state, j;
		t          = &hmm->transition[i];
		from_state = hmm->somap[t->from];
		state      = hmm->somap[t->to];
		dur        = hmm->dmap[random_genomic];
		genomic_null_trans = (zScoreDurationGroup(dur, dur->duration[0].min+2+zGetGenomicIncrement(hmm, state), gc) - zScoreDurationGroup(dur, dur->duration[0].min+2, gc));
		dur        = hmm->dmap[random_cdna];
		cdna_null_trans    = (zScoreDurationGroup(dur, dur->duration[0].min+2+zGetCDnaIncrement(hmm, state), gc)    - zScoreDurationGroup(dur, dur->duration[0].min+2, gc));
		if (zIsGenomicOnly(t->to)) {
			for (j = 0; j < hmm->iso_transitions; j++) 
				hmm->tmap[from_state][state][j] -= genomic_null_trans;
		} else if (zIsCDnaOnly(t->to)) {
			for (j = 0; j < hmm->iso_transitions; j++) 
				hmm->tmap[from_state][state][j] -= cdna_null_trans;
		} else if (zIsMatch(t->to)) {
			for (j = 0; j < hmm->iso_transitions; j++) 
				hmm->tmap[from_state][state][j] = hmm->tmap[from_state][state][j] - cdna_null_trans - genomic_null_trans; 
		}
	}

	for (i = 0; i < hmm->states; i++) {
		int j;
		dur = hmm->dmap[random_genomic];
		genomic_null_trans = (zScoreDurationGroup(dur, dur->duration[0].min+2+zGetGenomicIncrement(hmm, i), gc) - zScoreDurationGroup(dur, dur->duration[0].min+2, gc));
		dur = hmm->dmap[random_cdna];
		cdna_null_trans    = (zScoreDurationGroup(dur, dur->duration[0].min+2+zGetCDnaIncrement(hmm, i), gc) - zScoreDurationGroup(dur, dur->duration[0].min+2, gc));
		for (j = 0; j < hmm->iso_transitions; j++) {
			if (hmm->tmap[i][i][j] != MIN_SCORE) {
				if (zGetGenomicIncrement(hmm, i) > 0) {
					hmm->tmap[i][i][j] -= genomic_null_trans;
				} 
				if (zGetCDnaIncrement(hmm, i) > 0) {
					dur = hmm->dmap[random_cdna];
					hmm->tmap[i][i][j] -= cdna_null_trans;
				} 
			} else if (hmm->state[i].type == EXPLICIT) {
				zDistribution *dist = &hmm->dmap[hmm->state[i].duration]->duration[0].distribution[0];
				score_t step;
				coor_t k;

				if (zGetGenomicIncrement(hmm, i) > 0) {
					step = genomic_null_trans;
					for (k = dist->start; k <= dist->end; k++) {
						dist->param[k - dist->start] -= k*step;
					}
				}
				if (zGetCDnaIncrement(hmm, i) > 0) {
					step = cdna_null_trans;
					for (k = dist->start; k <= dist->end; k++) {
						dist->param[k - dist->start] -= k*step;
					}
				}
			}
		}
	}

	return 1;
}

static int zNullifyInitProbs (zHMM *hmm, const char* genomic_null_duration_name, const char* cdna_null_duration_name) {
	zStrIdx random_genomic = zChar2StrIdx(genomic_null_duration_name);
	zStrIdx random_cdna = zChar2StrIdx(cdna_null_duration_name);
	zDurationGroup *dur;
	score_t genomic_null_trans, cdna_null_trans, init_score;
	float gc = 50.0;
	int state_index, i, iso_group = 0;

	if (hmm->mode != GPAIRHMM) {
		zDie("Cannot Nullify initial probabilities for non-GPAIRHMM mode");
	}

	if (random_genomic >= hmm->feature_count) {
		return 0;
	}
	if (random_cdna >= hmm->feature_count) {
		return 0;
	}

	for (i = 0; i < hmm->iso_states; i++) {
		if (gc < hmm->iso_state[i]) {
			iso_group = i;
		}
	}

	dur        = hmm->dmap[random_genomic];
	genomic_null_trans = (zScoreDurationGroup(dur, 1, gc));
	dur        = hmm->dmap[random_cdna];
	cdna_null_trans    = (zScoreDurationGroup(dur, 1, gc));

	for (state_index = 0; state_index < hmm->states; state_index++) {
		init_score =  hmm->state[state_index].init[iso_group];
		init_score -= genomic_null_trans;
		init_score -= cdna_null_trans;
		hmm->state[state_index].init[iso_group] = init_score;
	}
	return 1;
}

static int zNullifyEmissions (zHMM *hmm, const char* suffix) {
	int i;

	if (hmm->mode != GPAIRHMM) {
		zDie("Cannot Nullify emissions for non-GPAIRHMM mode");
	}

	for (i = 0; i < hmm->models; i++) {
		int j;
		zStrIdx null_model_name;
		zModel *model, *null_model;
		int data_size=0;
		char buffer[64];

		model = &hmm->model[i];
		if (strstr(model->name, suffix) != NULL) continue;
		sprintf(buffer, "%s%s", model->name, suffix);
		null_model_name = zChar2StrIdx(buffer);
		if (null_model_name >= hmm->feature_count) {
			return 0;
		}
		null_model = hmm->mmap[null_model_name];
		if (model->type == LUT) {
			data_size = zPOWER[model->symbols][model->length];
		} else if (model->type == WMM) {
			data_size = model->symbols*model->length;
		} else {
			zDie("Cannot nullify non-LUT/WMM models yet");
		}
		for (j = 0; j < data_size; j++) {
			model->data[j] -= null_model->data[j];
		}
	}

	return 1;
}

int zNullifyHMM(zHMM* hmm) {
	if (hmm->mode != GPAIRHMM) {
		zDie("Cannot Nullify HMM for non-GPAIRHMM mode");
	}

	return (zNullifyEmissions(hmm, NULL_MODEL_SUFFIX) &&
		zNullifyInitProbs(hmm, GENOMIC_NULL_DURATION_NAME, CDNA_NULL_DURATION_NAME) &&
		zNullifyTransitions(hmm, GENOMIC_NULL_DURATION_NAME, CDNA_NULL_DURATION_NAME));
}

#ifdef DEBUG

int zGetGenomicIncrement(const zHMM* hmm, int state) {
	return hmm->increments[state][0];
}

int zGetCDnaIncrement(const zHMM* hmm, int state) {
	return hmm->increments[state][1];
}

score_t zGetTransitionScore(zHMM* hmm, int from_state, int to_state, int iso_group) {
	if (iso_group < 0 || iso_group >= hmm->iso_transitions ) {
		zDie("zGetTransitionScore: iso_group %d not found from %d to %d\n", 
			 iso_group, from_state, to_state );
	}
	return hmm->tmap[from_state][to_state][iso_group];
}

#endif /* DEBUG */

#endif /* ZOE_HMM_C */
