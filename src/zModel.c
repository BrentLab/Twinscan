/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zModel.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_MODEL_C
#define ZOE_MODEL_C

#include "zModel.h"

int zReadModel (FILE *stream, zModel *model) {
	char type[33];
	char name[33];
	char data[33];
	int length, focus, symbols, submodels, mem, i;
	float number;
	char seq_type[33];
	
	/* parse the model header */
	if (fscanf(stream, "%32s %32s %32s %d %d %d %d",
			name, type, seq_type, &length, &focus, &symbols, &submodels) != 7) {
		zWarn("zReadModel header failed (%s)", name);
		return 0;
	}


	/* set model attributes */
	model->length    = length;
	model->focus     = focus;
	model->order     = 0;
	model->symbols   = symbols;
	model->submodels = submodels;
	model->window_size = 0;
	model->name = zMalloc(strlen(name) + 1, "zModel name");
	(void)strcpy(model->name, name);
		
	     if (strcmp(type, "WMM")  == 0) model->type = WMM;
	else if (strcmp(type, "LUT")  == 0) model->type = LUT;
	else if (strcmp(type, "WWAM") == 0) model->type = WWAM;
	else if (strcmp(type, "SAM")  == 0) model->type = SAM;
	else if (strcmp(type, "SDT")  == 0) model->type = SDT;
	else if (strcmp(type, "CDS")  == 0) model->type = CDS;
	else if (strcmp(type, "MIX")  == 0) model->type = MIX;
	else if (strcmp(type, "ISO")  == 0) model->type = ISO;
 	else if (strcmp(type, "SIG")  == 0) model->type = SIG; 
	else {
		zWarn("zReadModel sequence model unknown (%s)", type);
		return 0;
	}

		 if      (strcmp(seq_type, "DNA")    == 0) model->seq_type = DNA;
		 else if (strcmp(seq_type, "CONSEQ") == 0) model->seq_type = CONSEQ;
		 else if (strcmp(seq_type, "ESTSEQ") == 0) model->seq_type = ESTSEQ;
		 else if (strcmp(seq_type, "GENOMIC")== 0) model->seq_type = GENOMIC;
		 else if (strcmp(seq_type, "PAIR") == 0)   model->seq_type = PAIR;
		 else {
			 zWarn("zReadModel sequence model type unknown (%s)", seq_type);
			 return 0;
		 }	
		 
	/* parse the model data */
	if (model->type == WMM) {
		model->submodel = NULL;
		mem = length * symbols;
		model->data = zMalloc(mem * sizeof(score_t), "zModel WMM data");
		for (i = 0; i < mem; i++) {
			if (fscanf(stream, "%32s", data) != 1) {
				zWarn("zReadModel wmm fscanf error");
				return 0;
			}
			model->data[i] = zText2Score(data);
		}

	} else if (model->type == SIG) {
		model->submodel = NULL;
		mem = 2*zPOWER[model->symbols][model->length];
		model->data = zMalloc(mem * sizeof(score_t), "zModel SIG data");
		for (i = 0; i < mem; i++) {
			if (fscanf(stream, "%32s", data) != 1) {
				zWarn("zReadModel sig fscanf error");
				return 0;
			}
			model->data[i] = zText2Score(data);	
		}
	} else if (model->type == LUT) {
		model->submodel = NULL;
		mem = zPOWER[model->symbols][model->length];
		model->data = zMalloc(mem * sizeof(score_t), "zModel LUT data");
		for (i = 0; i < mem; i++) {
			if (fscanf(stream, "%32s", data) != 1) {
				zWarn("zReadModel lut fscanf error");
				return 0;
			}
			model->data[i] = zText2Score(data);
	
		}
	} 

	else if (model->type == WWAM) {
		model->submodel = NULL;
	    /* read order */
		if (fscanf(stream, "%d", &model->order) != 1) {
			
			zWarn("zReadModel WWAM couldn't read order of model");
			return 0;
		}
		else {
			model->order += 1;
		}
		/* read data */
		mem = zPOWER[model->symbols][model->order] * model->length;
		model->data = zMalloc(mem * sizeof(score_t), "zModel WWAM data");
		for (i = 0; i < mem; i++) {
			if (fscanf(stream, "%32s", data) != 1) {
				zWarn("zReadModel wwam fscanf error");
				return 0;
			}
			
			/* 			printf("zReadModel:  reading WWAM model data index %i is %f\n", i, zText2Score(data)); */
			
			model->data[i] = zText2Score(data);	
			
		}
	}
	
	else { /* composite models */
		if (MIX == model->type) {
			mem = model->submodels;
			model->data = zMalloc (mem * sizeof(score_t), "zModel MIX data");
			for (i = 0; i < mem; i++ ) {
				if (fscanf(stream, "%f", &number) != 1) {
					zWarn("zReadModel mix fscanf error");
					return 0;
				}
				model->data[i] = zFloat2Score(number);
			}
		} else if (ISO == model->type) {
			/* check for GC data */
			(void)fscanf(stream, " GC_WINDOW=%u", &(model->window_size));
			mem = model->submodels;
			if (mem <= 0) {
				zWarn("zReadModel iso without any submodels");
				return 0;
			}
			model->data = zMalloc (mem * sizeof(score_t), "zModel ISO data");
/* 			if (fscanf(stream, "levels: %f", &model->data[0]) != 1) { */
/* 				zWarn("zReadModel iso couldn't find tag (levels)"); */
/* 				return 0; */
/* 			} */
			for (i = 0; i < mem; i++ ) {
				if (fscanf(stream, "%lf", &model->data[i]) != 1) {
					zWarn("zReadModel iso fscanf error");
					return 0;
				}
				model->data[i] /= 100.0;
			}
		}		
		model->submodel = zMalloc(sizeof(zModel) * model->submodels, "zModel");
		if ( SDT == model->type || CDS == model->type || SAM == model->type) {
			model->data = NULL;  /* I wonder what this was doing here unconditionally!!! Modified so that MIX and ISO will have data */
		}
		for (i = 0; i < model->submodels; i++) {
			if (!zReadModel(stream, &model->submodel[i])) {
				zWarn("zReadModel submodel");
				return 0;
			}
		}
	}
	
	return 1;
}

int zReadModelHeader (FILE *stream, zModel *model, int pseudocount) {
	char type[33];
	char name[33];
	int length, focus, symbols, submodels, mem, i;	
	
	/* parse the model header */
	if (fscanf(stream, "%32s %32s %d %d %d %d",
			name, type, &length, &focus, &symbols, &submodels) != 6) {

		zWarn("zReadModel header failed");

		return 0;
	}
	
	/* set model attributes */
	model->length    = length;
	model->focus     = focus;
	model->symbols   = symbols;
	model->submodels = submodels;
	model->name = zMalloc(strlen(name) + 1, "zModel name");
	(void)strcpy(model->name, name);
		
	     if (strcmp(type, "WMM") == 0) model->type = WMM;
	else if (strcmp(type, "LUT") == 0) model->type = LUT;
	else if (strcmp(type, "SAM") == 0) model->type = SAM;
	else if (strcmp(type, "SDT") == 0) model->type = SDT;
	else if (strcmp(type, "CDS") == 0) model->type = CDS;
	else if (strcmp(type, "MIX") == 0) model->type = MIX;
	else if (strcmp(type, "SIG") == 0) model->type = SIG;
	else {
		zWarn("zReadModel sequence model unknown (%s)", type);
		return 0;
	}
		
	/* parse the model data */
	if (model->type == WMM) {
		model->submodel = NULL;
		mem = length * symbols;
		model->data = zMalloc(mem * sizeof(score_t), "zModel WMM data");
		for (i = 0; i < mem; i++) model->data[i] = pseudocount;
	} else if (model->type == SIG) {
		model->submodel = NULL;
		mem = length * symbols;
		model->data = zMalloc(mem * sizeof(score_t), "zModel SIG data");
		for (i = 0; i < mem; i++) model->data[i] = pseudocount;
	} else if (model->type == LUT) {
		model->submodel = NULL;
		mem = zPOWER[model->symbols][model->length];
		model->data = zMalloc(mem * sizeof(score_t), "zModel LUT data");
		for (i = 0; i < mem; i++) model->data[i] = pseudocount;
	} else {
		model->submodel = zMalloc(sizeof(zModel) * model->submodels, "zModel");
		model->data = NULL;
		for (i = 0; i < model->submodels; i++) {
			if (!zReadModelHeader(stream, &model->submodel[i], pseudocount)) {
				zWarn("zReadModel submodel");
				return 0;
			}
		}
	}
	
	return 1;
}

static void zWriteModelPrivate(FILE *stream, const zModel *model, int indent) {
	char typestring[16], seqtypestring[16], data[16];
	int i, j, p;
	coor_t pos;
		
	switch (model->type) {
	case WMM:  strcpy(typestring, "WMM"); break;
	case LUT:  strcpy(typestring, "LUT"); break;
	case WWAM: strcpy(typestring, "WWAM"); break;
	case SAM:  strcpy(typestring, "SAM"); break;
	case SDT:  strcpy(typestring, "SDT"); break;
	case CDS:  strcpy(typestring, "CDS"); break;
	case MIX:  strcpy(typestring, "MIX"); break;
	case ISO:  strcpy(typestring, "ISO"); break;
	case SIG:  strcpy(typestring, "SIG"); break;
	default:   strcpy(typestring, "???"); break;
	}
	
	switch (model->seq_type) {
	case DNA:     strcpy(seqtypestring, "DNA");    break;
	case CONSEQ:  strcpy(seqtypestring, "CONSEQ"); break;
	case GENOMIC: strcpy(seqtypestring, "GENOMIC"); break;
	case PAIR:    strcpy(seqtypestring, "PAIR"); break;
	default:      strcpy(seqtypestring, "???");    break;
	}

	/* print header */
	for (i = 0; i < indent; i++) (void)fputc('\t', stream);
	(void)fprintf(stream, "%s %s %s %d %d %d %d\n",
				  model->name, typestring, seqtypestring, (int)model->length, 
				  (int)model->focus, model->symbols,	model->submodels);
	
	/* print body */
	if (model->type == WMM) {
		for (pos = 0; pos < model->length; pos++) {
			for (j = 0; j <= indent; j++) (void)fputc('\t', stream); /* padding */
			for (j = 0; j < model->symbols; j++) {
				zScore2Text(model->data[ (pos * model->symbols) + j], data);
				(void)fprintf(stream, "%s\t", data);
			}
			(void)fprintf(stream, "\n");
		}
	} else if (model->type == LUT) {
		p = zPOWER[model->symbols][model->length];
		for (i = 0; i < p; i++) {
			if (i % model->symbols == 0) { /* padding */
				for (j = 0; j <= indent; j++) (void)fputc('\t', stream);
			}
			zScore2Text(model->data[i], data);
			(void)fprintf(stream, "%s\t", data);
			if ( (i+1) % model->symbols == 0) (void)fprintf(stream, "\n");
		}
	} else if (model->type == WWAM) {
		/* print order */
		for (j = 0; j <=indent ; j++) fputc('\t', stream);
		(void)fprintf(stream, "%d\n", model->order - 1);
		/* print data */
		p = zPOWER[model->symbols][model->order] * model->length;
		for (i = 0; i < p; i++) {
			if (i % model->symbols == 0) { /* pad */
				for (j = 0; j <= indent; j++) (void)fputc('\t', stream);
			}
			zScore2Text(model->data[i], data);
			fprintf(stream, "%s\t", data);
			if ( (i+1) % model->symbols == 0) fprintf(stream, "\n");
		}
	} else if (model->type == SIG) {
		p = 2*zPOWER[model->symbols][model->length];
		for( i = 0; i < p; i++) {
			if (i % model->symbols == 0) { /* padding */
				for (j = 0; j <= indent; j++) (void)fputc('\t', stream);
			}
			zScore2Text(model->data[i], data);
			(void)fprintf(stream, "%s\t", data);
			if ( (i+1) % model->symbols == 0) (void)fprintf(stream, "\n");
		}
	} else {
		for (i = 0; i < model->submodels; i++) {
			(void)zWriteModelPrivate(stream, &model->submodel[i], indent +1);
		}
	}
}

void zWriteModel(FILE *stream, const zModel *model) {
	zWriteModelPrivate(stream, model, 0); /* just sets indent to zero */
}

void zFreeModel (zModel *model) {
	int i;

	for (i = 0; i < model->submodels; i++) zFreeModel(&model->submodel[i]);
        if (model->submodel != NULL) zFree(model->submodel);
	zFree(model->name);
	zFree(model->data);
}


void zAmbiguateModel (zModel *model, score_t score) {

	/*
	  When scoring real DNA sequences you have to deal with Ns
	  but when defining a model you only want to define it on
	  the 4 DNA bases.  This function takes a model which is 
	  defined for 4 symbols A, C, G, and T and transforms it 
	  into a model defined on 5 symbols A, C, G, T, and N.  The
	  score input to the function is used to fill in scores of
	  sequences containing Ns. (actually different scores are
	  given to Ns on different models to correspond to twinscan
	  check the model you are interested in to make sure what
	  gets used)
	*/
	
	int i, j, index, offset_4, offset_5;
	size_t mem;
	score_t *d5;
	char string[50];
	
	if (model->symbols == 5) {
		return; /* Nothing to do */
	} else if (model->symbols !=4) {
		zDie("zAmbiguateModel requires 4-symbol models");
	}
	model->symbols = 5;
	
	if (model->type == WMM) {
		mem = model->length * 5 * sizeof(score_t);
		d5 = zMalloc(mem, "zAmbiguateModel WMM data");
		
		for (i = 0; i < (int)model->length; i++) {
			/* record the 4-symbol scores in the 5-symbol array */
			for (j = 0; j < 4; j++) {
				d5[(i * 5) + j] = model->data[(i * 4) + j];
			}
			/* set all N's to score unless the model length is 1, in
			which case the impossible feature score is used */
			if (model->length == 1) d5[(i * 5) + 4] = MIN_SCORE;
			else                    d5[(i * 5) + 4] = score;
		}
		zFree(model->data);
		model->data = d5;
		
	} else if (model->type == SIG) {

		mem = 2*zPOWER[5][model->length] * sizeof(score_t);
		d5 = zMalloc(mem, "zAmbiguateModel SIG data");
		
		/* set all values to score to begin */
		for (i = 0; i < 2*zPOWER[5][model->length]; i++) d5[i] = 0.0;
		
		/* this reflects what genscan does to the signal peptide if gives */
		/* sequences with an n a score of 0  */
		/* see CState.cpp line 467 */

		/* convert all 4-symbol values into the 5-symbol array */
		for (i = 0; i < 2*zPOWER[4][model->length]; i++) {
			zDecToBase(i, 4, string);      /* convert to base4 string */
			index = zBaseToDec(5, string); /* convert to base5 string */
			d5[index] = model->data[i];
		}
		zFree(model->data);
		model->data = d5;	

	} else if (model->type == LUT) {
		mem = zPOWER[5][model->length] * sizeof(score_t);
		d5 = zMalloc(mem, "zAmbiguateModel LUT data");
		
		/* set all values to score to begin */

/* 		for (i = 0; i < zPOWER[5][model->length]; i++) d5[i] = score; */
		for (i = 0; i < zPOWER[5][model->length]; i++) d5[i] = -10.0;
		
		/* convert all 4-symbol values into the 5-symbol array */
		for (i = 0; i < zPOWER[4][model->length]; i++) {
			zDecToBase(i, 4, string);      /* convert to base4 string */
			index = zBaseToDec(5, string); /* convert to base5 string */
			d5[index] = model->data[i];
		}
		zFree(model->data);
		model->data = d5;
		
	} 

	else if (model->type == WWAM) {
		
		mem = zPOWER[5][model->order] * model->length * sizeof(score_t);
		d5 = zMalloc(mem, "zAmbiguateModel WWAM data");

		/* and then a miracle occurs... */

/* 		for(i=0;i < zPOWER[5][model->order] * model->length;i++) d5[i] = score; */

/* this reflects what happens in twinscan, it assigns -1 to ns */
/* see CState.cpp line 252 */

		for(i=0;i < zPOWER[5][model->order] * (int)model->length;i++) {
			d5[i] = -1.0;
		}

		for(i=0; i < (int)model->length; i++) {
			offset_5 = i*zPOWER[5][model->order];
			offset_4 = i*zPOWER[4][model->order];
			for(j=0;j<zPOWER[4][model->order];j++) {
				zDecToBase(j, 4, string);      /* convert to base4 string */
				index = zBaseToDec(5, string); /* convert to base5 string */
				d5[index+offset_5] = model->data[j+offset_4];

/* 				printf("zAmbiguateModel:  offset_5 = %i offset_4 = %i i = %i moving %f from index %i in 4-symbol rep to index %i in 5-symbol rep\n", offset_5, offset_4, i, model->data[j+offset_4], j+offset_4, index+offset_5); */

			}

		}

		zFree(model->data);
		model->data = d5;
	}

	else {
		for (i = 0; i < model->submodels; i++) {
			zAmbiguateModel(&model->submodel[i], score);
		}
	}
}

void zDeambiguateModel (zModel *model) {
	int i, index;
	coor_t j;
	size_t mem;
	score_t *d4;
	char string[50];
	int  containsN;
			
	if (model->symbols != 5) {
		zDie("zDeambiguateModel requires 5-symbol models");
	}
	model->symbols = 4;
		
	if (model->type == WMM) {
		mem = model->length * 4 * sizeof(score_t);
		d4 = zMalloc(mem, "zDeambiguateModel WMM data");
		
		for (i = 0; i < (int)model->length; i++) {
			for (j = 0; j < 4; j++) {
				d4[(i * 4) + j] = model->data[(i * 5) + j];
			}
		}
		zFree(model->data);
		model->data = d4;
		
	} else if (model->type == LUT) {
		mem = zPOWER[4][model->length] * sizeof(score_t);
		d4 = zMalloc(mem, "zDeambiguateModel LUT data");
				
		/* convert all 5-symbol values into the 4-symbol array */
		for (i = 0; i < zPOWER[5][model->length]; i++) {
			zDecToBase(i, 5, string); /* base5 string */
			containsN = 0;
			for (j = 0; j < strlen(string); j++) {
				if (string[j] == '4') {
					containsN = 1;
					break;
				}
			}
			if (containsN) continue; /* skip N's */
			index = zBaseToDec(4, string);
			d4[index] = model->data[i];
		}
		zFree(model->data);
		model->data = d4;

	} else if (model->type == WWAM) {
		mem = zPOWER[4][model->order] * model->length;
		d4 = zMalloc(mem * sizeof(score_t), "zDeambiguate WWAM data");
		
		for (i = 0; i < zPOWER[5][model->order]* (int)model->length; i++) {
			zDecToBase(i % zPOWER[5][model->order], 5, string); /* base5 str */
			containsN = 0;
			for (j = 0; j < strlen(string); j++) containsN |= (string[j]=='4');
			if (containsN) continue; /* skip N's */
			index = zBaseToDec(4, string)
				+ zPOWER[4][model->order]*(i/zPOWER[5][model->order]);
			d4[index] = model->data[i];
		}
		zFree(model->data);
		model->data = d4;

	} else if (model->type == SIG) {
		mem = 2 * zPOWER[4][model->length] * sizeof(score_t);
		d4 = zMalloc(mem, "zDeambiguate SIG data");

		/* convert all 5-symbol values into the 4-symbol array */
		for (i = 0; i < 2 * zPOWER[5][model->length]; i++) {
			zDecToBase(i, 5, string); /* base5 string */
			containsN = 0;
			for (j = 0; j < strlen(string); j++) {
				if (string[j] == '4') {
					containsN = 1;
					break;
				}
			}
			if (containsN) continue; /* skip N's */
			index = zBaseToDec(4, string);
			d4[index] = model->data[i];
		}
		zFree(model->data);
		model->data = d4;


	} else {
		for (i = 0; i < model->submodels; i++) {
			zDeambiguateModel(&model->submodel[i]);
		}
	}
}

int zAntiModel (zModel* model) {
	score_t *buffer;
	char string[50];
	int length, i, to;
	int base = model->symbols;
	int stringlength;
	int j;
	if (model->type != LUT && model->type != WMM) {
		zDie("zAntiModel can only be done on LUT and WMM models");
	}
	if (model->symbols != 4 && model->symbols != 5) {
		zDie("zAntiModel can only be done on a 4 or 5 symbol model");
	}
	/* To reverse complement the S5 int of 4 or 5 symbol model -> (8-x)%5 */
	if (model->type == LUT) {
		length = zPOWER[model->symbols][model->length];
		buffer = (score_t*) zMalloc(sizeof(score_t)*length, "zAntiModel: buffer");

		/* Complement and reverse - reverse only for length > 1 */
		for (i = 0; i < length; i++) {
			/* Get the string */
			zDecToBase(i, base, string);
			stringlength = (int)strlen(string);

			/* Pad it with 0s in the left and make it model->length long */
			string[model->length]='\0';
			for (j = stringlength - 1; j >= 0; j--) {
				string[j+model->length-stringlength] = string[j];
			}
			for (j = 0; j < (int)model->length-stringlength; j++) {
				string[j] = '0';
			}
			stringlength = (int)strlen(string);

			/* Complement */
			for (j = 0; j < stringlength; j++) {
				string[j] = '0'+((8-(string[j]-'0'))%5);
			}
			if (model->seq_type != PAIR) {
				/* Reverse */
				for (j = 0; j < stringlength/2; j++) {
					char c = string[j];
					string[j] = string[stringlength-1-j];
					string[stringlength-1-j] = c;
				}
			}
			to = zBaseToDec(base, string);
			buffer[to] = model->data[i];
		}
		for (i = 0; i < length; i++) {
			model->data[i] = buffer[i];
		}
		zFree(buffer);
	} else if (model->type == WMM) {
		length = model->symbols*model->length;
		buffer = (score_t*) zMalloc(sizeof(score_t)*length, "zAntiModel: buffer");
		for (i = 0; i < (int)model->length; i++) {
			for (j = 0; j < model->symbols; j++) {
				int to   = i*model->symbols + j;
				int from = (model->length - 1 - i)*model->symbols + (8 - j)%5;
				buffer[to] = model->data[from];
			}
		}
		for (i = 0; i < length; i++) {
			model->data[i] = buffer[i];
		}
		zFree(buffer);
	}
	return 1;
}

void zMarginalizeModel (zModel *model, int base) {
	int max = base - 1;
	if (model->symbols != base) {
		zDie("Ambiguate the model to %d before calling zMarginalizeModel", base);
	}
	if (model->type == LUT) {
		int i;
		char string[50];
		for (i = 0; i < zPOWER[base][model->length]; i++) {
			float sum = 0;
			int position = -1;
			int index = -1;
			coor_t j = 0;
			zDecToBase(i, base, string);      /* convert to base_base string */
			for (j = 0; j < strlen(string); j++) {
				if (string[j] == '0' + max) {
					position = j;
					break;
				}
			}
			if (position != -1) {
				for (j = 0; j < (coor_t)max; j++) {
					string[position] = '0'+j;
					index = zBaseToDec(base, string); /* convert to base5 string */
					sum += zScore2Float(model->data[index]);
				}
				string[position] = '0' + max;
				index = zBaseToDec(base, string); /* convert to base5 string */
				model->data[index] = zFloat2Score(sum/(max));
			}
		}
	} else if (model->type == WMM) {
		coor_t i;
		for (i = 0; i < model->length; i++) {
			float sum = 0;
			coor_t j;
			for (j = 0; j < (coor_t) max; j++) {
				sum += zScore2Float(model->data[i*model->symbols+j]);
			}
			model->data[i*model->symbols+max] = zFloat2Score(sum/max);
		}
	} else {
		zDie("Can only marginalize LUT and WMM models!");
	}
}

#endif
