/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
zGTF.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2003 Charles J. Vaske

\******************************************************************************/

#include "zGTF.h"

/******************************************************************************\
  gtf_t
\******************************************************************************/
#define GTF_START_STRING  "start_codon"
#define GTF_STOP_STRING   "stop_codon"
#define GTF_CDS_STRING    "CDS"
#define GTF_EXON_STRING   "exon"
#define GTF_ESNGL_STRING  "zEsngl"
#define GTF_5UTR_STRING   "5UTR"
#define GTF_UNDEF_STRING  "undefined"

gtf_t  zText2GTF(const char* s) {
    if (!strcmp(s, GTF_START_STRING)) {
		return zGTF_START;
    } else if (strcmp(s, GTF_STOP_STRING) == 0) {
		return zGTF_STOP;
    } else if (strcmp(s, GTF_CDS_STRING) == 0) {
		return zGTF_CDS;
    } else if (strcmp(s, GTF_EXON_STRING) == 0) {
		return zGTF_EXON;
    } else if (strcmp(s, GTF_EXON_STRING) == 0) {
		return zGTF_EXON;
    } else if (strcmp(s, GTF_ESNGL_STRING) == 0) {
		return zGTF_ESNGL;
    } else if (strcmp(s, GTF_5UTR_STRING) == 0) {
		return zGTF_5UTR;
    } else {
		zWarn("Unidentified GTF feature: %s", s);
		return zGTF_UNDEF;
    }
}

void zGTF2Text(const gtf_t ftype, char* s) {
    switch (ftype) {
    case zGTF_START: strcpy(s, GTF_START_STRING); break;
    case zGTF_STOP:  strcpy(s, GTF_STOP_STRING);  break;
    case zGTF_CDS:   strcpy(s, GTF_CDS_STRING);   break;
    case zGTF_EXON:  strcpy(s, GTF_EXON_STRING);  break;
	case zGTF_ESNGL: strcpy(s, GTF_ESNGL_STRING); break;
	case zGTF_5UTR:  strcpy(s, GTF_5UTR_STRING);  break;
    case zGTF_UNDEF:
    default:         strcpy(s, GTF_UNDEF_STRING);
    }
}

/******************************************************************************\
  GTF Features
\******************************************************************************/
#define GTF_READ_FORMAT "%127s\t%127s\t%31s\t%9s\t%9s\t%9s\t%1s\t%1s\tgene_id \"%127[^\"]\"; transcript_id \"%127[^\"]\";"
#define GTF_WRITE_FORMAT "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\";\n"

int zReadGTFfeature (FILE *stream, zGTFfeature *f) {
	coor_t  i;
	int  number_read;
	char name[128], source[128], type[32], start[10], end[10],
	    score[10], strand[2], frame[2], gene_id[128], transcript_id[128];
	char line[1024];

	do { /* loop over blank lines */
		if (fgets(line, sizeof(line)-1, stream) == NULL) return 0;
		line[sizeof(line)-1] = '\0'; /* being extra careful... */
	
		/* Remove any comments on the line */
		for (i = 0; i < sizeof(line); i++) {
			if ('\0' == line[i] || '#'  == line[i] || '\n' == line[i]) {
				line[i] = '\0'; break; 
			}
		}
		
		/* read all the attributes */
		number_read = sscanf(line, GTF_READ_FORMAT,
							 name, source, type, start, end, score, 
							 strand, frame, gene_id, transcript_id);
	} while (0 == number_read || EOF == number_read); /* blank line */
	if (number_read != 10) {
		return 0;
	}
		
	/* alloc memory for the struct */
	f->seqname       = zMalloc(strlen(name)+1,   "zReadGTFFeature name");
	f->source        = zMalloc(strlen(source)+1, "zReadGTFFeature source");
	f->gene_id       = zMalloc(strlen(gene_id)+1,"zReadGTFFeature gene_id");
	f->transcript_id = zMalloc(strlen(transcript_id)+1,
							   "zReadGTFFeature transcript_id");
	
	/* fill in the struct */
	strcpy(f->seqname, name);
	strcpy(f->source,  source);
	f->type          = zText2GTF(type);
	f->start         = zText2Coor(start)-1; /* we use 0-based coords */
	f->end           = zText2Coor(end)-1;   /* we use 0-based coords */
	f->score         = zText2Score(score);
	f->strand        = zText2Strand(strand);
	f->frame         = zText2Frame(frame);
	strcpy(f->gene_id, gene_id);
	strcpy(f->transcript_id, transcript_id);
	return 1;
}

void zClearGTFfeature (zGTFfeature *f) {
	f->seqname       = NULL;
	f->source        = NULL;
	f->type          = zGTF_UNDEF;
	f->start         = UNDEFINED_COOR;
	f->end           = UNDEFINED_COOR;
	f->score         = MIN_SCORE;
	f->strand        = UNDEFINED_STRAND;
	f->frame         = UNDEFINED_FRAME;
	f->gene_id       = NULL;
	f->transcript_id = NULL;
}

void zFreeGTFfeature (zGTFfeature *f) {
	zFree(f->seqname);
	zFree(f->source);
	zFree(f->gene_id);
	zFree(f->transcript_id);
	zClearGTFfeature(f);
}

void zWriteGTFfeature(FILE *stream, zGTFfeature *f) {
	char type[32], start[10], end[10], score[10], strand[2], frame[2];
	
	zGTF2Text(f->type, type);
	zCoor2Text(f->start+1, start); /* convert to 1-based coords */
	zCoor2Text(f->end+1, end);     /* convert to 1-based coords */
	zScore2Text(f->score, score);
	zStrand2Text(f->strand, strand);
	zFrame2Text(f->frame, frame);

	fprintf(stream, GTF_WRITE_FORMAT, f->seqname, f->source, type,
			start, end, score, strand, frame, f->gene_id, f->transcript_id);
}

void zCopyGTFfeature (const zGTFfeature* src, zGTFfeature* f) {
	f->seqname = NULL;
	if (NULL != src->seqname) {
		f->seqname = zMalloc(strlen(src->seqname)+1,"zCopyGTFfeature seqname");
		strcpy(f->seqname, src->seqname);
	}
	f->source = NULL;
	if (NULL != src->source) {
		f->source = zMalloc(strlen(src->source)+1, "zCopyGTFfeature source");
		strcpy(f->source, src->source);
	}
	f->type   = src->type;
	f->start  = src->start;
	f->end    = src->end;
	f->score  = src->score;
	f->strand = src->strand;
	f->frame  = src->frame;
	if (NULL != src->gene_id) {
		f->gene_id = zMalloc(strlen(src->gene_id)+1,"zCopyGTFfeature gene_id");
		strcpy(f->gene_id, src->gene_id);
	}
	if (NULL != src->transcript_id) {
		f->transcript_id = zMalloc(strlen(src->transcript_id)+1,
								   "zCopyGTFfeature transcript_id");
		strcpy(f->transcript_id, src->transcript_id);
	}

}

int zGTFfeatureCmp (const zGTFfeature* a, const zGTFfeature* b) {
	int r;
	/* compare first based on the transcricpt, then gene id's */
	if ((r = strcmp(a->transcript_id, b->transcript_id))) return r;
	if ((r = strcmp(a->gene_id, b->gene_id)))             return r;

	/* order by coordinates */
	if      (a->end < b->end) return -1;
	else if (a->end > b->end) return  1;
	if      (a->start < b->start) return -1;
	else if (a->start > b->start) return  1;

	/* place the "smaller" features first  (start, stop, CDS, then exon)*/
	if      (a->type < b->type) return -1;
	else if (a->type > b->type) return  1;

	return 0;
}

int zGTFpfeatureCmp (const void* a, const void* b) {
	return zGTFfeatureCmp((zGTFfeature*)a, (zGTFfeature*)b);
}

/******************************************************************************\
  zGTFVec
\******************************************************************************/

void zInitGTFVec (zGTFVec *vec, size_t limit) {
	vec->size  = 0;
	vec->limit = limit;
	vec->elem  = (vec->limit > 0) ? zMalloc(limit*sizeof(zGTFfeature), "GTFVec")
		: NULL;
	vec->last  = NULL;
}

void zPushGTFVec (zGTFVec *vec, const zGTFfeature *feature) {
	if (vec->limit == vec->size) {
		vec->limit = (vec->limit == 0) ? 1 : vec->limit*2;
		vec->elem = zRealloc(vec->elem, vec->limit*sizeof(zGTFfeature),
							 "zPushGTFVec");
	}
	vec->last = &vec->elem[vec->size];
	zCopyGTFfeature(feature, vec->last);
	vec->size++;
}

void zFreeGTFVec (zGTFVec *vec) {
	int i;
	for (i = 0; i < vec->size; ++i) {
		zFreeGTFfeature(&vec->elem[i]);
	}
	zFree(vec->elem); 
	vec->elem  = NULL;
	vec->size  = 0;
	vec->limit = 0;
	vec->last  = NULL;
}

zGTFVec* zReadGTFVec (FILE* stream) {
	zGTFVec     *vec;
	zGTFfeature  f;

	vec = zMalloc(sizeof(zGTFVec), "zReadGTFVec");
	zInitGTFVec(vec, 1);
	while (zReadGTFfeature(stream, &f)) {
		zPushGTFVec(vec, &f);
		zFreeGTFfeature(&f);
	}
	return vec;
}

void zWriteGTFVec (FILE* stream, const zGTFVec *vec) {
	int i;
	for (i = 0; i < vec->size; ++i) {
		zWriteGTFfeature(stream, &vec->elem[i]);
	}
}

void zSortGTFVec (zGTFVec *vec) {
        if (vec->size > 0)
	    qsort(vec->elem, vec->size, sizeof(zGTFfeature), zGTFpfeatureCmp);
}

/******************************************************************************\
  zGTFConversion
\******************************************************************************/

void zInitGTFConvInfo(zGTFConvInfo *conv) {
	memset(conv->gtft2StrIdx, -1, sizeof(conv->gtft2StrIdx));
	conv->maxStrIdx = zStringPoolCount();
	conv->strIdx2gtft = zMalloc(conv->maxStrIdx * sizeof(gtf_t),
								"zReadGTFConvInfo");
	memset(conv->strIdx2gtft, zGTF_UNDEF, conv->maxStrIdx * sizeof(gtf_t));
}

void zFreeGTFConvInfo(zGTFConvInfo *conv) {
	zFree(conv->strIdx2gtft);
	conv->strIdx2gtft = NULL;
}

void zReadGTFConvInfo(FILE* stream, zGTFConvInfo *conv) {
	char    line[2048], model[128], gtff[128];
	gtf_t   gtf_type;
	zStrIdx stridx;
	int     number_read;

	zInitGTFConvInfo(conv);
	while (fgets(line, sizeof(line)-1, stream) != NULL) {
		line[sizeof(line)-1] = '\0';

		number_read = sscanf(line, "%127s => %127s", model, gtff);
		model[sizeof(model)-1] = '\0';
		gtff [sizeof(gtff) -1] = '\0';

		if (EOF == number_read) continue;
		if (2 != number_read) return;
		stridx   = zChar2StrIdx(model);
		gtf_type = zText2GTF(gtff);
		if (stridx < 0 || stridx >= conv->maxStrIdx) {
			zWarn("The model type %s is unrecognized", model);
		} else if (gtf_type >= zGTF_UNDEF) {
			zWarn("GTF feature %s is not allowed", gtff);
		} else {
			/* enter the conversion */
			conv->gtft2StrIdx[gtf_type] = stridx;
			conv->strIdx2gtft[stridx]   = gtf_type;
		}
	}
}

void zWriteGTFConvInfo(FILE* stream, const zGTFConvInfo *conv) {
	char   gtf_type[32];
	int    i;

	for (i = 0; i < conv->maxStrIdx; ++i) {
		if (conv->strIdx2gtft[i] != zGTF_UNDEF) {
			zGTF2Text(conv->strIdx2gtft[i], gtf_type);
			fprintf(stream, "%s => %s\n", zStrIdx2Char(i), gtf_type);
		}
	}
}


int zCDSMatchesStop(zGTFfeature *cds, zGTFfeature *stop) {
	if (cds->type != zGTF_CDS || stop->type != zGTF_STOP) return 0;
	
	if (cds->strand != stop->strand) return 0;
	
	if (strcmp(cds->gene_id, stop->gene_id)) return 0;

	if (cds->strand == '+') {
		return (cds->end+1 == stop->start);
	} else if (cds->strand == '-') {
		return (cds->start-1 == stop->end);
	} else {
		return 0;
	}
}

int zCDSMatchesStart(zGTFfeature *cds, zGTFfeature *start) {
	if (cds->type != zGTF_CDS || start->type != zGTF_START) return 0;

	if (cds->strand != start->strand) return 0;

	if (strcmp(cds->gene_id, start->gene_id)) return 0;
	
	if (cds->strand == '+') {
		return (cds->start == start->start);
	} else if (cds->strand == '-') {
		return (cds->end == start->end);
	} else {
		return 0;
	}
}
