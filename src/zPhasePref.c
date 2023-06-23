/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zPhasePref.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_PHASEPREF_C
#define ZOE_PHASEPREF_C

#include "zPhasePref.h"

int zReadPhasePref (FILE *stream, zPhasePref *pp) {
	int i;
	
	for (i = 0; i < PHASECOUNT; i++) {
		if (fscanf(stream, "%f", &pp->prob[i]) != 1) {
			zWarn("zReadPhasePref fscanf error");
			return 0;
		}
		/* convert to score too */
		pp->score[i] = zFloat2Score(pp->prob[i]);
	}
	
	pp->Esngl = (zStrIdxExists("Esngl")) ? zChar2StrIdx("Esngl") : -1;
	pp->Einit = (zStrIdxExists("Einit")) ? zChar2StrIdx("Einit") : -1;
	pp->Exon  = (zStrIdxExists("Exon"))  ? zChar2StrIdx("Exon")  : -1;
	pp->Eterm = (zStrIdxExists("Eterm")) ? zChar2StrIdx("Eterm") : -1;
	
	return 1;
}

void zWritePhasePref (FILE *stream, const zPhasePref *pp) {
	int i;
	
	for (i = 0; i < PHASECOUNT; i++) {
		(void)fprintf(stream, "%f\n", pp->prob[i]);
	}
}

/*
 * Score the phase change when transitioning from an exon
 */
score_t zScorePhaseFromE (zPhasePref *pp, zStrIdx from, zPhase_t to, int lfrag) {
	if (from == pp->Einit) {
		switch (to) {
			case Phase0:                               return pp->score[Ei_I0];
			case Phase1: case Phase1T:                 return pp->score[Ei_I1];
		    case Phase2: case Phase2TA: case Phase2TG: return pp->score[Ei_I2];
		    default: return 0.0; /*zDie("zScorePhase doesFromE not allow to %s", zStrIdx2Char(to)); */
		}
	} else if (from == pp->Exon) {
		switch (lfrag) {
			case 0: /* corresponds to E0 */
				switch (to) {
					case Phase0:                  return pp->score[E0_I0];
					case Phase1: case Phase1T:    return pp->score[E0_I1];
					case Phase2: case Phase2TA: 
					case Phase2TG:                return pp->score[E0_I2];	
					default: zDie("zScorePhase impossible A %d", to);
				}
			case 1: /* corresponds to E2 */
				switch (to) {
					case Phase0:                  return pp->score[E2_I0];
					case Phase1: case Phase1T:    return pp->score[E2_I1];
					case Phase2: case Phase2TA: 
					case Phase2TG:                return pp->score[E2_I2];
					default: zDie("zScorePhase impossible B %d", to);
				}
			case 2: /* corresponds to E1 */
				switch (to) {
					case Phase0:                  return pp->score[E1_I0];
					case Phase1: case Phase1T:    return pp->score[E1_I1];
					case Phase2: case Phase2TA: 
					case Phase2TG:                return pp->score[E1_I2];
					default: zDie("zScorePhase impossible C %d", to);
				}
			default: zDie("zScorePhase impossible D");
		}
	} 
	return 0;
}

/*
 * Score the phase change when transitioning to an exon
 */
score_t zScorePhaseToE (zPhasePref *pp, zPhase_t from, zStrIdx to) {
	switch (from) {
		case Phase0:
			     if (to == pp->Exon)  return pp->score[I0_E0];
			else if (to == pp->Eterm) return pp->score[I0_Et];
			else if (to == pp->Einit) return 0;
            else return 0.0; /*zDie("zScorePhaseToE does not allow to %s", zStrIdx2Char(to));*/
		case Phase1: case Phase1T:
			     if (to == pp->Exon)  return pp->score[I1_E1];
			else if (to == pp->Eterm) return pp->score[I1_Et];
			else if (to == pp->Einit) return 0;
            else return 0.0; /*zDie("zScorePhaseToE does not allow to %s", zStrIdx2Char(to)); */
		case Phase2: case Phase2TA: case Phase2TG:
			     if (to == pp->Exon)  return pp->score[I2_E2];
			else if (to == pp->Eterm) return pp->score[I2_Et];
			else if (to == pp->Einit) return 0;
            else return 0.0; /*zDie("zScorePhaseToE does not allow to %d", zStrIdx2Char(to)); */
		default: return 0;
	}
	
	/*return 0;*/

}

#endif
