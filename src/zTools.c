/* -*- Mode: C; tab-width: 2; c-basic-offset: 2; indent-tabs-mode: t -*- */  
/******************************************************************************\
 zTools.c - part of the ZOE library for genomic analysis
 
 Copyright (C) 2001-2002 Ian F. Korf

\******************************************************************************/

#ifndef ZOE_TOOLS_C
#define ZOE_TOOLS_C

#include "zTools.h"
#include <libgen.h>

const coor_t   UNDEFINED_COOR   = UINT_MAX;
const frame_t  UNDEFINED_FRAME  = -1;
/*EVAN score_t is now a double not a float */
const score_t  MIN_SCORE        = -FLT_MAX;
const score_t  MAX_SCORE        = FLT_MAX;
const strand_t UNDEFINED_STRAND = -1;

/******************************************************************************\
 Library Information
\******************************************************************************/

static char zVersionNumber[] = "2002-03-25";
void zLibInfo (void) {
	(void)fprintf(stderr, "ZOE library version %s\n", zVersionNumber);
}


/******************************************************************************\
 Program Name
\******************************************************************************/

#define ZPROG_NAME_LN 256
static char PROGRAM_NAME[ZPROG_NAME_LN] = "unnamed program";

void zSetProgramName (char *string) {
	char* prog_name = basename(string);
	if (NULL == prog_name) prog_name = string;
	(void)strncpy(PROGRAM_NAME, prog_name, ZPROG_NAME_LN);
	PROGRAM_NAME[ZPROG_NAME_LN-1] = '\0';
}

char* zGetProgramName (void) {
	return PROGRAM_NAME;
}

/******************************************************************************\
 Commandline Processing
\******************************************************************************/

static zTVec zARGUMENTS;
static zHash zOPTIONS;
static char zTRUE_OPTION[5] = "true";

void zParseOptions (int *argc, char **argv) {
	int    i;
	size_t j;
	char  *token, *attr, *value;
	
	zInitHash(&zOPTIONS);
	zInitTVec(&zARGUMENTS, 1);
	
	/* parse */
	for (i = 0; i < *argc; i++) {
		token = argv[i];
		attr  = token+1;
		value = zTRUE_OPTION;
		if (token[0] == '-' && strlen(token) > 1) {
			for (j = 0; j < strlen(token); j++) {
				if (token[j] == '=') {
					token[j] = '\0';
					value = token + j + 1;
					break;
				}
			}
			zSetHash(&zOPTIONS, attr, value);
		} else {
			zPushTVec(&zARGUMENTS, argv[i]);
		}
	}
	
	/* reset the command line */
	*argc = zARGUMENTS.size;
	for (i = 0; i < zARGUMENTS.size; i++) {
		argv[i] = zARGUMENTS.elem[i];
	}
}

void zFreeOptions() {
	zFreeTVec(&zARGUMENTS);
	zFreeHash(&zOPTIONS);
}

char* zOption (const char *attr) {
	return (char*)zGetHash(&zOPTIONS, attr);
}


/******************************************************************************\
 Library Verbosity

Specifies verbosity levels for library

\******************************************************************************/

static zVerbosityObject* LIBVERBOSITY = NULL;

void zSetVerbosityLevel(int verbosity) {

  if(LIBVERBOSITY == NULL) {
	LIBVERBOSITY = zMalloc(sizeof(zVerbosityObject), "zSetVerbosityLevel");
	LIBVERBOSITY->verbosity_level = verbosity;
	zInitTVec(&LIBVERBOSITY->features, 10);

  }
  else {
	LIBVERBOSITY->verbosity_level = verbosity;
  }
}

void zFreeVerbosityGlobalVariable() {
	zFreeTVec(&LIBVERBOSITY->features);
	zFree(LIBVERBOSITY);
}

void zInitVerbosityObject(zVerbosityObject* verbosity_object, int verbosity) {
	verbosity_object = zMalloc(sizeof(zVerbosityObject), "zInitVerbosityObject");
	verbosity_object->verbosity_level = verbosity;
	zInitTVec(&verbosity_object->features, 10);
}

void zFreeVerbosityObject(zVerbosityObject* verbosity) {
	zFreeTVec(&verbosity->features);
}

int zGetVerbosityLevel() {
  return (LIBVERBOSITY->verbosity_level);
}

void zAddVerbosityFeature(char* s) {

  zPushTVec(&LIBVERBOSITY->features, s);

}

zTVec* zGetVerbosityFeatures() {

  return &LIBVERBOSITY->features;

}

/******************************************************************************\
 Important typedef input/output
\******************************************************************************/

static void reverse_string (char *s) {
	int c, i, j;
	
	for (i = 0, j = strlen(s) -1; i < j; i++, j--) {
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}
}

static void itoa (int n, char *s) {
	int i, sign;
	
	if ((sign = n) < 0) n = -n;
	i = 0;
	do {
		s[i++] = n % 10 + '0';
	} while ((n /= 10) > 0);
	if (sign < 0)
		s[i++] = '-';
		s[i] = '\0';
		reverse_string(s);
}

/* coor_t */
void zCoor2Text (coor_t val, char *s) {
	if (val == UNDEFINED_COOR) strcpy(s, "...");
	else                       itoa(val, s);
}

coor_t zText2Coor (const char *s) {
	int val;
	
	if (strcmp(s, "...") == 0) return UNDEFINED_COOR;
	if ( (sscanf(s, "%d", &val)) != 1) {
		zWarn("zText2Coor sscanf failure (%s)", s);
		return UNDEFINED_COOR;
	}
	return (coor_t)val;
}

/* frame_t */
void zFrame2Text (frame_t val, char *s) {
	if (val == UNDEFINED_FRAME) strcpy(s, ".");
	else                        itoa(val, s);
}

frame_t zText2Frame (const char *s) {
	int val;
	
	if (strcmp(s, ".") == 0) return UNDEFINED_FRAME;
	if ( (sscanf(s, "%d", &val)) != 1) {
		zWarn("zText2Frame sscanf failure (%s)", s);
		return UNDEFINED_FRAME;
	}
	return (frame_t)val;
}

/* score_t */
void zScore2Text (score_t val, char *s) {
	     if (val == MIN_SCORE) strcpy(s, ".");
	else if (val == MAX_SCORE)  strcpy(s, "*");
	else                              sprintf(s, "%.0f", val);
}

void zScore2FloatText (score_t val, char *s) {
	     if (val == MIN_SCORE) strcpy(s, ".");
	else if (val == MAX_SCORE)  strcpy(s, "*");
	else                              sprintf(s, "%f", val);
}

score_t zText2Score (const char *s) {
	int val;
	
	     if (strcmp(s, ".") == 0) return MIN_SCORE;
	else if (strcmp(s, "*") == 0) return MAX_SCORE;
	if ( (sscanf(s, "%d", &val)) != 1) {
		zWarn("zText2Score sscanf failure (%s)", s);
		return MIN_SCORE;
	}
	return (score_t)val;
}

/* strand_t */
void zStrand2Text (strand_t val, char *s) {
	     if (val == '+') strcpy(s, "+");
	else if (val == '-') strcpy(s, "-");
	else                 strcpy(s, ".");
}

strand_t zText2Strand (const char *s) {	
	if (strcmp(s, ".") == 0) return UNDEFINED_STRAND;
	if (strcmp(s, "+") == 0) return '+';
	if (strcmp(s, "-") == 0) return '-';
	
	zWarn("zText2Strand sscanf failure (%s)", s);
	return UNDEFINED_STRAND;
}

#define PHASE0   "0"
#define PHASE1   "1"
#define PHASE1T  "1T"
#define PHASE2   "2"
#define PHASE2TA "2TA"
#define PHASE2TG "2TG"

zPhase_t zText2Phase (const char *s) {
	if (!strcmp(s, PHASE0))   return Phase0;
	if (!strcmp(s, PHASE1))   return Phase1;
	if (!strcmp(s, PHASE1T))  return Phase1T;
	if (!strcmp(s, PHASE2))   return Phase2;
	if (!strcmp(s, PHASE2TA)) return Phase2TA;
	if (!strcmp(s, PHASE2TG)) return Phase2TG;

	zWarn("zText2Phase unrecognized phase");
	return Phase0;
}

void zPhase2Text (const zPhase_t phase, char *s) {
	switch (phase) {
	case Phase0:   strcpy(s, PHASE0);   break;
	case Phase1:   strcpy(s, PHASE1);   break;
	case Phase1T:  strcpy(s, PHASE1T);  break;
	case Phase2:   strcpy(s, PHASE2);   break;
	case Phase2TA: strcpy(s, PHASE2TA); break;
	case Phase2TG: strcpy(s, PHASE2TG); break;
	default:       strcpy(s, "U");
	}
}

/******************************************************************************\
 Warnings and Errors
\******************************************************************************/

void zWarn (const char *fmt, ...) {
	va_list args;
	
	fprintf(stderr, "ZOE WARNING (from %s): ", zGetProgramName());
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	
	if (errno) {
		fprintf(stderr, "ERRNO %d: %s\n", errno, strerror(errno));
		errno = 0;
	}
}

void zDie (const char *fmt, ...) {
	va_list args;
	
	fflush(stdout);
	
	fprintf(stderr, "ZOE ERROR (from %s): ", zGetProgramName());
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	
	if (errno) fprintf(stderr, "ERRNO %d: %s\n", errno, strerror(errno));

	exit(2);
}


/******************************************************************************\
 Tracing
\******************************************************************************/

#if TRACE == 1
	void zTrace1 (const char *fmt, ...) {
		va_list args;
	
		if (errno) fprintf(stderr, "ERRNO %d: %s\n", errno, strerror(errno));
		fprintf(stderr, "TRACE1: ");
		va_start(args, fmt);
		vfprintf(stderr, fmt, args);
		va_end(args);
		fprintf(stderr, "\n");
	}
void zTrace2 (const char *fmt, ...) {const char*a=fmt;a=0;/*compiler hush*/}
void zTrace3 (const char *fmt, ...) {const char*a=fmt;a=0;/*compiler hush*/}

#elif TRACE == 2
	void zTrace1 (const char *fmt, ...) {
		va_list args;
	
		if (errno) fprintf(stderr, "ERRNO %d: %s\n", errno, strerror(errno));
		fprintf(stderr, "TRACE1: ");
		va_start(args, fmt);
		vfprintf(stderr, fmt, args);
		va_end(args);
		fprintf(stderr, "\n");
	}
	void zTrace2 (const char *fmt, ...) {
		va_list args;
	
		if (errno) fprintf(stderr, "ERRNO %d: %s\n", errno, strerror(errno));
		fprintf(stderr, "\tTRACE2: ");
		va_start(args, fmt);
		vfprintf(stderr, fmt, args);
		va_end(args);
		fprintf(stderr, "\n");
	}
void zTrace3 (const char *fmt, ...) {const char*a=fmt;a=0;/*compiler hush*/}

#elif TRACE >= 3
	void zTrace1 (const char *fmt, ...) {
		va_list args;
	
		if (errno) fprintf(stderr, "ERRNO %d: %s\n", errno, strerror(errno));
		fprintf(stderr, "TRACE1: ");
		va_start(args, fmt);
		vfprintf(stderr, fmt, args);
		va_end(args);
		fprintf(stderr, "\n");
	}
	void zTrace2 (const char *fmt, ...) {
		va_list args;
	
		if (errno) fprintf(stderr, "ERRNO %d: %s\n", errno, strerror(errno));
		fprintf(stderr, "\tTRACE2: ");
		va_start(args, fmt);
		vfprintf(stderr, fmt, args);
		va_end(args);
		fprintf(stderr, "\n");
	}
	void zTrace3 (const char *fmt, ...) {
		va_list args;
	
		if (errno) fprintf(stderr, "ERRNO %d: %s\n", errno, strerror(errno));
		fprintf(stderr, "\t\tTRACE3: ");
		va_start(args, fmt);
		vfprintf(stderr, fmt, args);
		va_end(args);
		fprintf(stderr, "\n");
	}

#else
void zTrace1 (const char *fmt, ...) { const char*a=fmt;a=0;/*compiler hush*/}
void zTrace2 (const char *fmt, ...) { const char*a=fmt;a=0;/*compiler hush*/}
void zTrace3 (const char *fmt, ...) { const char*a=fmt;a=0;/*compiler hush*/}

#endif


/******************************************************************************\
 Memory Tools
\******************************************************************************/

void* zMalloc (size_t size, const char *str) {
	void *buffer;
	
	zTrace3("zMalloc (%d) from %s", size, str);
	if ((buffer = malloc(size)) == NULL) zDie("malloc(%d) %s", size, str);
	return buffer;
}

void* zCalloc (size_t nobj, size_t size, const char *str) {
	void *buffer;

	zTrace3("zCalloc (%d, %d) from %s", nobj, size, str);
	if ((buffer = calloc(nobj, size)) == NULL) zDie("calloc %s", str);
	return buffer;  
}

void* zRealloc (void *p, size_t size, const char *str) {
	void *buffer;
	
	zTrace3("zRealloc (%d) from %s", size, str);
	if (p == NULL) {
		buffer = zMalloc(size, "zRealloc redirection to zMalloc");
	} else {
		buffer = realloc(p, size);
		if (buffer == NULL) zDie("realloc %s", str);
	}
	return buffer;  
}

void zFree (void *p) {
	if (p != NULL) {
		free(p);
		p = NULL;
	}
}

#ifdef __G_LIB_H__

/* GCC slice malloc */
gpointer zSliceAlloc (gsize size, const char *str) {
	gpointer buffer;

	zTrace3("zSliceAlloc (%d) from %s", size, str);
	if ((buffer = g_slice_alloc(size)) == NULL) zDie("g_slice_alloc(%d) %s", size, str);
	return buffer;
}

void zSliceFree (gsize size, gpointer p) {
	if (p != NULL) {
		g_slice_free1(size, p);
		p = NULL;
	}
}

#else

void* zSliceAlloc (size_t size, const char *str) {
	void *buffer;
	
	zTrace3("zSliceAlloc (%d) from %s", size, str);
	if ((buffer = malloc(size)) == NULL) zDie("malloc(%d) %s", size, str);
	return buffer;
}

void zSliceFree (size_t size, void* p) {
	if (p != NULL) {
		if (size > 0) free(p);
		p = NULL;
	}
}

#endif /* __G_LIB_H__ */

/******************************************************************************\
 Comparison Functions
\******************************************************************************/

int zTcmp (const void *a, const void *b) {
	return strcmp( *(char**)a, *(char**)b );
}
int zIcmp (const void *a, const void *b) {
	return *(int*)a - *(int*)b;
}
int zFcmp (const void *a, const void *b) {
	float f = *(float*)a - *(float*)b;
	if (f > 0) return  1;
	else if (f < 0) return -1;
	else            return  0;
}

/******************************************************************
 * -1 can occur only at the end of traceback, where it all starts.* 
 * Therefore I assume that -1 is less than everything else. But   *
 * (coor_t) -1 = MAX_UNSIGNED_INT, and it will break comparisons. *
 * Hence (coor_t) -1 is treated as a special case.                * 
 ******************************************************************/
int zCoortCmp (const void *a, const void *b) {
	coor_t minus_one = (coor_t) -1;
	coor_t c1 = *(coor_t*)a;
	coor_t c2 = *(coor_t*)b;

	if (c1 == minus_one && c2 == minus_one) {
		return 0;
	}
	if (c1 == minus_one) {
		return -1;
	} else if (c2 == minus_one) {
		return 1;
	}
	if (c1 < c2) {
		return -1;
	} else if (c1 > c2) {
		return 1;
	} else {
		return 0;
	}
}

/******************************************************************************\
 zIVec
\******************************************************************************/
void zFreeIVec (zIVec *vec) {
	zFree(vec->elem);
}
void zInitIVec (zIVec *vec, int limit) {
	vec->size  = 0;
	vec->limit = limit;
	if (vec->limit > 0) vec->elem = zMalloc(limit * sizeof(int), "zInitIVec");
	else                vec->elem = NULL;
}
void zPushIVec (zIVec *vec, int val) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zRealloc(vec->elem, vec->limit * sizeof(int), "zPushIVec");
	}
	vec->elem[vec->size] = val;
	vec->last = val;
	vec->size++;
}


/******************************************************************************\
 zFVec
\******************************************************************************/
void zFreeFVec (zFVec *vec) {
	zFree(vec->elem);
}
void zInitFVec (zFVec *vec, int limit) {
	vec->size  = 0;
	vec->limit = limit;
	if (vec->limit > 0) vec->elem = zMalloc(limit * sizeof(float), "zInitFVec");
	else                vec->elem = NULL;
}
void zPushFVec (zFVec *vec, float val) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zRealloc(vec->elem, vec->limit * sizeof(float), "zPushFVec");
	}
	vec->elem[vec->size] = val;
	vec->last = val;
	vec->size++;
}


/******************************************************************************\
 zTVec
\******************************************************************************/
void zFreeTVec (zTVec *vec) {
	int i;
	for (i = 0; i < vec->size; i++) {
		zFree(vec->elem[i]);
	}
	zFree(vec->elem);
}
void zInitTVec (zTVec* vec, int limit) {
	vec->size  = 0;
	vec->limit = limit;
	vec->last  = NULL;
	if (vec->limit > 0) vec->elem = zMalloc(limit * sizeof(char*), "zInitTVec");
	else                vec->elem = NULL;
}
void zPushTVec (zTVec* vec, const char* text) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zRealloc(vec->elem, vec->limit * sizeof(char*), "zPushTVec");
	}
	vec->elem[vec->size] = zMalloc(strlen(text) +1, "zPushTVec");
	strcpy(vec->elem[vec->size], text);
	vec->last = vec->elem[vec->size];
	vec->size++;
}


/******************************************************************************\
 zVec
\******************************************************************************/
void zFreeVec (zVec *vec) {
	zFree(vec->elem);
}
void zInitVec (zVec *vec, int limit) {
	vec->size  = 0;
	vec->limit = limit;
	if (vec->limit > 0) vec->elem = zMalloc(limit * sizeof(void*), "zInitVec");
	else                vec->elem = NULL;
	vec->last  = NULL;
}
void zPushVec (zVec *vec, void *thing) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zRealloc(vec->elem, vec->limit * sizeof(void*), "zPushVec");
	}
	vec->elem[vec->size] = thing;
	vec->last = vec->elem[vec->size];
	vec->size++;
}



/******************************************************************************\
 zHash

The hash function is my own creature. I used the multiplication method as
described  in Cormen et al. but added a 7 periodic component so that for
example, ATG and TGA don't hash to the same index. Because DNA might be hashed,
I thought it would be a good idea to avoid some multipe of 3 and because number
systems are 2 or 10 based, I stayed away from those too. The 7 values chosen
were kind of arbitrary. I have tested this on some rather large text files, and
the hash function separation is quite good even after several levels of
re-hashing. Performance is good too, about 6x faster than the C++ map and 3
times faster than a Perl hash.

\******************************************************************************/

static double zHASH_MULTIPLIER[7] = {
	3.1415926536, /* PI */
	2.7182818285, /* e */
	1.6180339887, /* golden mean */
	1.7320508076, /* square root of 3 */
	2.2360679775, /* square root of 5 */
	2.6457513111, /* square root of 7 */
	3.3166247904, /* square root of 11 */
};

int zHashFunc (const zHash *hash, const char *key) {
	size_t i;
	double sum;

	sum = 0;
	for (i = 0; i < strlen(key); i++) {
		sum += key[i] * zHASH_MULTIPLIER[i % 7];
	}
	
	return (int) (hash->slots * (sum - floor(sum)));
}

static float zMAX_HASH_DEPTH = 2.0;         /* The hash will remain between */
static int zHashLevelToSlots (int level) {  /* half full and twice-filled */
	return pow(4, level);                   /* with these values */
}

static void zExpandHash (zHash *hash) {
	int   i, j;
	char *key;
	void *val;
	int   oldslots = hash->slots;
	zVec *oldkey   = hash->key;
	zVec *oldval   = hash->val;
	zVec *kvec, *vvec;
		
	/* create the new hash */
	hash->level = hash->level +1;
	hash->slots = zHashLevelToSlots(hash->level);
	hash->key   = zMalloc(hash->slots * sizeof(zVec), "zExpandHash key");
	hash->val   = zMalloc(hash->slots * sizeof(zVec), "zExpandHash val");
	for (i = 0; i < hash->slots; i++) {
		zInitVec(&hash->key[i], 0);
		zInitVec(&hash->val[i], 0);
	}
	
	/* brand new hash? */
	if (hash->keys == 0) return;
	else                 hash->keys = 0; /* will be set below */
	
	/* transfer old stuff to new hash and free old key */
	for (i = 0; i < oldslots; i++) {
		kvec = &oldkey[i];
		vvec = &oldval[i];
		for (j = 0; j < kvec->size; j++) {
			key = kvec->elem[j];
			val = vvec->elem[j];
			zSetHash(hash, key, val);
			zFree(key);
		}
		zFreeVec(kvec);
		zFreeVec(vvec);
	}
	
	zFree(oldkey);
	zFree(oldval);
}

void zInitHash (zHash *hash) {
	hash->level = 0;
	hash->slots = 0;
	hash->keys  = 0;
	hash->key   = NULL;
	hash->val   = NULL;
	zExpandHash(hash);
}

void zFreeHash (zHash *hash) {
	int i, j;

	for (i = 0; i < hash->slots; i++) {
		for (j = 0; j < hash->key[i].size; j++) {
			zFree(hash->key[i].elem[j]);
		}
		zFreeVec(&hash->key[i]);
		zFreeVec(&hash->val[i]);
	}
	zFree(hash->key);
	zFree(hash->val);
}

void* zGetHash (const zHash *hash, const char *key) {
	int i, index;
	char *string;

	index = zHashFunc(hash, key);
	/* resolve collisions */
	for (i = 0; i < hash->key[index].size; i++) {
		string = hash->key[index].elem[i];
		if (strcmp(key, string) == 0) {
			return hash->val[index].elem[i];
		}
	}
	return NULL; /* return is NULL if not found */
}

void zSetHash (zHash *hash, const char *key, void *val) {
	int i, index;
	char *string;
	char *actual_key; /* MA: Added this so that on the fly keys can  *
				*     be used. Otherwise, freeing the    *
				*     hash won't free the value corres-  *
				*     ponding to that key. Querying may  *
				*     return empty set too!!             *
				*     Also fixed zFreeHash() and         *
				*     zExpandHash() to recognize this    */
	int new_key = 1;

	index = zHashFunc(hash, key);

	/* reassign unless new key */
	for (i = 0; i < hash->key[index].size; i++) {
		string = hash->key[index].elem[i];
		if (strcmp(key, string) == 0) {
			hash->val[index].elem[i] = val;
			new_key = 0;
			break;
		}
	}

	if (new_key) {
		actual_key = (char*) zMalloc(sizeof(char)*(strlen(key) + 1), "zSetHash for key copy");
		strcpy(actual_key, key);
		zPushVec(&hash->val[index], val);
		zPushVec(&hash->key[index], (void*)actual_key);
		hash->keys++;
	}

	/* check if we have to expand the hash */
	if ((float)hash->keys / (float)hash->slots >= zMAX_HASH_DEPTH) {
		zExpandHash(hash);
	}
}

zVec* zKeysOfHash (const zHash *hash) {
	int i, j;
	zVec *vec;
	vec = zMalloc(sizeof(zVec), "zKeysOfHash");

	zInitVec(vec, hash->keys);
	for (i = 0; i < hash->slots; i++) {
		for (j = 0; j < hash->key[i].size; j++) {
			zPushVec(vec, hash->key[i].elem[j]);
		}
	}
	
	return vec;
}

zVec* zValsOfHash (const zHash *hash) {
	int i, j;
	zVec *vec;
	vec = zMalloc(sizeof(zVec), "zValsOfHash");

	zInitVec(vec, hash->keys);
	for (i = 0; i < hash->slots; i++) {
		for (j = 0; j < hash->val[i].size; j++) {
			zPushVec(vec, hash->val[i].elem[j]);
		}
	}
	return vec;
}

void zStatHash (const zHash *hash) {
	int i, max, min, total, count;

	max = 0;
	min = 1000000;
	total = 0;
	for (i = 0; i < hash->slots; i++) {
		count = hash->val[i].size;
		total += count;
		if (count > max) max = count;
		if (count < min) min = count;
	}
	(void)printf("HashStats: level=%d slots=%d keys=%d min=%d max=%d ave=%f\n",
		hash->level, hash->slots, hash->keys, min, max,
		(float)total / (float)hash->slots);
}

void zEStatHash (const zHash *hash) {
	int i, max, min, total, count;

	max = 0;
	min = 1000000;
	total = 0;
	for (i = 0; i < hash->slots; i++) {
		count = hash->val[i].size;
		total += count;
		if (count > max) max = count;
		if (count < min) min = count;
	}
	(void)fprintf(stderr, "HashStats: level=%d slots=%d keys=%d min=%d max=%d ave=%f",
		hash->level, hash->slots, hash->keys, min, max,
		(float)total / (float)hash->slots);
	
	for (i = 0; i < hash->slots; i++) {
		(void)fprintf(stderr, " %d", hash->val[i].size);
	}
	fprintf(stderr, "\n");
}

/******************************************************************************\
  String Interning
\******************************************************************************/

static zHash *zINTERNED = NULL;
static zHash *zStringPool   = NULL;
static zVec  *zStringLookup = NULL;

char* zStrIntern (char* str, int freeStr) {
	char* val;
	
	/* No interns yet */

	if (zINTERNED == NULL) {
		zINTERNED = zMalloc(sizeof(zHash), "Init StrIntern Hash");
		zInitHash(zINTERNED);
	}
	
	val = (char*)zGetHash(zINTERNED, str);
	
	/* The string has already been interned */

	if (val != NULL) {
		return val;
	}

	/* Create new copy for interning */
	/* This will be returned to zStrIdx2Char where it
           will be put in the zStringLookup vector. That 
           should be freed by zFree() in zStringPoolFree(). */

	val = zMalloc(strlen(str)+1, "StrIntern key");
	strcpy(val, str);
	
	/* Set that in the intern hash and return it */

	zSetHash(zINTERNED, val, val);
	
	if (freeStr) zFree(str);
	
	return val;
}

static void zStringPoolInit(void) {
	zStringPool   = zMalloc(sizeof(zHash), "String Pool init");
	zStringLookup = zMalloc(sizeof(zVec), "String Pool init");
	
	zInitHash (zStringPool);
	zInitVec  (zStringLookup, 10); /* size is arbitrary */
}

/* The values in zStringPool are the idx's for the strings. They are not freed elsewhere */

void zStringPoolFree(void) {
	int i;
	zVec *values;

	/* Internalized copies of strings (See zStrIntern) */
	for (i = 0; i < zStringLookup->size; i++) {
		zFree(zStringLookup->elem[i]);
	}
	zFreeVec(zStringLookup); 
	zFree(zStringLookup);

	values = zValsOfHash(zStringPool);
	for (i = 0; i < values->size; i++) {
		zFree(values->elem[i]);
	}
	zFreeVec(values);
	zFree(values);
	zFreeHash (zStringPool);
	zFree(zStringPool);

	zFreeHash(zINTERNED);
	zFree(zINTERNED);
}

int zStringPoolCount(void) {
	return (NULL == zStringLookup) ? 0 : zStringLookup->size;
}

zStrIdx zChar2StrIdx(const char* str) {
	char* key;
	zStrIdx *  idx;

	if (NULL == zStringPool) zStringPoolInit();

	idx = (zStrIdx*)zGetHash(zStringPool, str);

	if (NULL == idx) { /* enter new string */
		key = zStrIntern((char*)str, 0);

		zPushVec(zStringLookup, key);
		idx  = zMalloc(sizeof(zStrIdx), "String Pool index");
		*idx = zStringLookup->size - 1;

		zSetHash(zStringPool, key, idx);
	}
	
	return *idx;
}

int zStrIdxExists(const char* str) {
	return (NULL != zStringPool && NULL != zGetHash(zStringPool, str));
}

const_string_t zStrIdx2Char(const zStrIdx idx) {
	return (NULL == zStringLookup
			|| idx < 0
			|| idx >= zStringLookup->size) ? NULL : zStringLookup->elem[idx];
}

/**********************************************************************\
    zList
\**********************************************************************/

zListNode* zListGetNode(zList* l){
	zListNode* ln;
	if((l->dead == NULL) || (l->dead->next == NULL)){
		ln = zMalloc(sizeof(zListNode),"zGetNode ln");		
		ln->data = l->init_func();
	}
	else{
		ln = l->dead->next;
		l->dead->next = ln->next;
		if(ln->next != NULL){
			ln->next->prev = l->dead;
		}
		l->reset_func(ln->data);
	}
	ln->next = NULL;
	ln->prev = NULL;
	return ln;
}

void zListReleaseNode(zList* l,zListNode* ln){
	ln->next = l->dead->next;
	ln->prev = l->dead;
	if(l->dead->next != NULL){
		l->dead->next->prev = ln;
	}
	l->dead->next = ln;
}

void  zInitList(zList* l,zListInitFunc init_func,zListFreeFunc free_func,
								zListResetFunc reset_func){
  l->init_func = init_func;
	l->free_func = free_func;
	l->reset_func = reset_func;
	l->head = zMalloc(sizeof(zListNode),"zInitList ln");		
	l->tail = zMalloc(sizeof(zListNode),"zInitList ln");		
	l->dead = zMalloc(sizeof(zListNode),"zInitList ln");		
	l->current = l->head;
	l->head->prev = NULL;
	l->head->next = l->tail;
	l->tail->prev = l->head;
	l->tail->next = NULL;
	l->dead->prev = NULL;
	l->dead->next = NULL;
	l->head->data = NULL;
	l->tail->data = NULL;
	l->dead->data = NULL;
	l->size = 0;
}

void  zFreeList(zList* l){
	l->current = l->head->next;
	while(l->current != l->tail){
		l->current = l->current->next;
		l->free_func(l->current->prev->data);
		zFree(l->current->prev);
	}
	zFree(l->head);
	zFree(l->tail);
	l->head = NULL;
	l->tail = NULL;
	l->size = 0;
	l->current = l->dead->next;
	if(l->current != NULL){
		while(l->current->next != NULL){
			l->current = l->current->next;
			l->free_func(l->current->prev->data);
			zFree(l->current->prev);
		}
		l->free_func(l->current->data);
		zFree(l->current);
	}
	zFree(l->dead);
	l->dead = NULL;
	l->current = NULL;
}

void  zResetList(zList* l){
	l->current = l->head->next;
	while(l->current != l->tail){
		l->current = l->current->next;
		zListReleaseNode(l,l->current->prev);
	}
	l->head->next = l->tail;
	l->tail->prev = l->head;
	l->current = l->head;
	l->size = 0;
}

int zListGetSize(zList* l){
	return l->size;
}

void* zListMoveFirst(zList* l){
	l->current = l->head->next;
	return l->current->data;
}

void* zListMoveLast(zList* l){
	l->current = l->tail->prev;
	return l->current->data;
}

void* zListMoveNext(zList* l){
	if(l->current != l->tail){
		l->current = l->current->next;
	}
	return l->current->data;
}

void* zListMovePrev(zList* l){
	if(l->current != l->head){
		l->current = l->current->prev;
	}
	return l->current->data;
}

void* zListGetCurrent(zList* l){
	return l->current->data;
}

bool  zListHasNext(zList* l){
	if((l->current != l->tail) &&
		 (l->current->next != l->tail)){
		return true;
	}
	else{
		return false;
	}
}
															
bool  zListHasPrev(zList* l){
	if((l->current != l->head) &&
		 (l->current->prev != l->head)){
		return true;
	}
	else{
		return false;
	}
}

void* zListAddNext(zList* l){
	zListNode* ln;
	if(l->current == l->tail){
		return zListAddPrev(l);
	}
	ln = zListGetNode(l);
	ln->prev = l->current;
	ln->next = l->current->next;
	l->current->next = ln;
	ln->next->prev = ln;
	l->size++;
	return ln->data;
}

void* zListAddPrev(zList* l){
	zListNode* ln;
	if(l->current == l->head){
		return zListAddNext(l);
	}
	ln = zListGetNode(l);
	ln->next = l->current;
	ln->prev = l->current->prev;
	l->current->prev = ln;
	ln->prev->next = ln;
	l->size++;
	return ln->data;
}

void* zListAddFirst(zList* l){
	zListNode* ln = zListGetNode(l);
	ln->prev = l->head;
	ln->next = l->head->next;
	l->head->next = ln;
	ln->next->prev = ln;
	l->size++;
	return ln->data;
}

void* zListAddLast(zList* l){
	zListNode* ln = zListGetNode(l);
	ln->next = l->tail;
	ln->prev = l->tail->prev;
	l->tail->prev = ln;
	ln->prev->next = ln;
	l->size++;
	return ln->data;
}

void zListRemoveFirst(zList* l){
	zListNode* ln = l->head->next;
	if(ln == l->tail){
		return;
	}
	if(ln == l->current){
		l->current = ln->next;
	}
	l->head->next = ln->next;
	ln->next->prev = l->head;
	l->size--;
	zListReleaseNode(l,ln);
}

void zListRemoveLast(zList* l){
	zListNode* ln = l->tail->prev;
	if(ln == l->head){
		return;
	}
	if(ln == l->current){
		l->current = ln->prev;
	}
	l->tail->prev = ln->prev;
	ln->prev->next = l->tail;
	l->size--;
	zListReleaseNode(l,ln);
}

void zListRemoveCurrent(zList* l){
	zListNode* n = l->current;
	if(n == l->head || n == l->tail){
		return;
	}
	n->next->prev = n->prev;
	n->prev->next = n->next;
	l->current = n->prev;
	l->size--;
	zListReleaseNode(l,n);
}

zListPos zListGetPos(zList* l){
	return l->current;
}

void zListSetPos(zList* l,zListPos p){
	l->current = p;
}

void zListMovePosToNext(zList* l,zListPos p){
	zListNode* n = p;
	zListNode* c = l->current;

	if(c == NULL || c == n || n == NULL || n == l->head || n == l->tail){
		return;
	}
	n->prev->next = n->next;
	n->next->prev = n->prev;
	if(c == l->tail){
		c = c->prev;
	}
	n->next = c->next;
	n->prev = c;
	c->next->prev = n;
	c->next = n;
}

void zListMovePosToPrev(zList* l,zListPos p){
	zListNode* n = p;
	zListNode* c = l->current;
	
	if(c == NULL || c == n || n == NULL || n == l->head || n == l->tail){
		return;
	}
	n->prev->next = n->next;
	n->next->prev = n->prev;
	if(c == l->head){
		c = c->next;
	}
	n->next = c;
	n->prev = c->prev;
	c->prev->next = n;
	c->prev = n;
}

void  zListMoveCurrentToFirst(zList* l){
	if(l->current == NULL || l->current == l->head || l->current == l->tail || 
		 l->current == l->head->next){
		return;
	}
	l->current->prev->next = l->current->next;
	l->current->next->prev = l->current->prev;
	l->head->next->prev = l->current;
	l->current->next = l->head->next;
	l->head->next = l->current;
	l->current->prev = l->head;
}

void  zListMoveCurrentToLast(zList* l){
	if(l->current == NULL || l->current == l->head || l->current == l->tail || 
		 l->current == l->tail->prev){
		return;
	}
	l->current->prev->next = l->current->next;
	l->current->next->prev = l->current->prev;
	l->tail->prev->next = l->current;
	l->current->prev = l->tail->prev;
	l->tail->prev = l->current;
	l->current->next = l->tail;
}

void* zListGetFirst(zList* l){
	if(l->head->next == l->tail){
		return NULL;
	}
	return l->head->next->data;
}

void* zListGetLast(zList* l){
	if(l->tail->prev == l->head){
		return NULL;
	}
	return l->tail->prev->data;
}

/**********************************************************************\
    zPtrList
\**********************************************************************/

/* list of dead nodes ready for reuse by any zPtrList linked only with
   next pointers */
zPtrListNode* dead_plnodes = NULL;
/* count of active lists. when all lists are freed we need to free all 
   dead nodes */
int plcount = 0;

#ifdef DEBUG
int ptrlistnode_id = 0;
#endif

zPtrListNode* zPtrListGetNode(){
	zPtrListNode* ln;
	if(dead_plnodes == NULL){
		ln = zMalloc(sizeof(zPtrListNode),"zPtrListGetNode ln");		
#ifdef DEBUG
		ln->id = ptrlistnode_id++;
#endif
	}
	else{
		ln = dead_plnodes;
		dead_plnodes = dead_plnodes->next;
	}
	ln->next = NULL;
	ln->prev = NULL;
	ln->data = NULL;
	return ln;
}

void zPtrListReleaseNode(zPtrListNode* ln){
	ln->next = dead_plnodes;
	dead_plnodes = ln;
}

#ifdef DEBUG
int ptrlist_id = 0;
#endif

void  zInitPtrList(zPtrList* l){
	l->head = zPtrListGetNode();
	l->tail = zPtrListGetNode();
	l->current = l->head;
	l->head->next = l->tail;
	l->tail->prev = l->head;
	l->size = 0;
#ifdef DEBUG
	l->id = ptrlist_id++;
#endif
	plcount++;
}

void  zFreePtrList(zPtrList* l){
	l->current = l->head->next;
	while(l->current != l->tail){
		l->current = l->current->next;
		zFree(l->current->prev);
	}
	zFree(l->head);
	zFree(l->tail);
	plcount--;
	if(plcount == 0){
		l->current = dead_plnodes;
		while(l->current != NULL){
			l->head = l->current;
			l->current = l->current->next;
			zFree(l->head);
		}
	}
	l->head = NULL;
	l->tail = NULL;
	l->current = NULL;
	l->size = 0;
}

void  zResetPtrList(zPtrList* l){
	l->current = l->head->next;
	while(l->current != l->tail){
		l->current = l->current->next;
		zPtrListReleaseNode(l->current->prev);
	}
	l->head->next = l->tail;
	l->tail->prev = l->head;
	l->current = l->head;
	l->size = 0;
}

int zPtrListGetSize(zPtrList* l){
	return l->size;
}

void* zPtrListMoveFirst(zPtrList* l){
	l->current = l->head->next;
	return l->current->data;
}

void* zPtrListMoveLast(zPtrList* l){
	l->current = l->tail->prev;
	return l->current->data;
}

void* zPtrListMoveNext(zPtrList* l){
	if(l->current != l->tail){
		l->current = l->current->next;
	}
	return l->current->data;
}

void* zPtrListMovePrev(zPtrList* l){
	if(l->current != l->head){
		l->current = l->current->prev;
	}
	return l->current->data;
}

void* zPtrListGetCurrent(zPtrList* l){
	return l->current->data;
}

bool  zPtrListHasNext(zPtrList* l){
	if((l->current != l->tail) &&
		 (l->current->next != l->tail)){
		return true;
	}
	else{
		return false;
	}
}
															
bool  zPtrListHasPrev(zPtrList* l){
	if((l->current != l->head) &&
		 (l->current->prev != l->head)){
		return true;
	}
	else{
		return false;
	}
}

void zPtrListAddNext(zPtrList* l,void* ptr){
	zPtrListNode* ln;
	if(l->current == l->tail){
		zPtrListAddPrev(l,ptr);
		return;
	}
	ln = zPtrListGetNode();
	ln->data = ptr;
	ln->prev = l->current;
	ln->next = l->current->next;
	l->current->next = ln;
	ln->next->prev = ln;
	l->size++;
}

void zPtrListAddPrev(zPtrList* l,void* ptr){
	zPtrListNode* ln;
	if(l->current == l->head){
		zPtrListAddNext(l,ptr);
		return;
	}
	ln = zPtrListGetNode();
	ln->data = ptr;
	ln->next = l->current;
	ln->prev = l->current->prev;
	l->current->prev = ln;
	ln->prev->next = ln;
	l->size++;
}

void zPtrListAddFirst(zPtrList* l,void* ptr){
	zPtrListNode* ln = zPtrListGetNode();
	ln->data = ptr;
	ln->prev = l->head;
	ln->next = l->head->next;
	l->head->next = ln;
	ln->next->prev = ln;
	l->size++;
}

void zPtrListAddLast(zPtrList* l,void* ptr){
	zPtrListNode* ln = zPtrListGetNode();
	ln->data = ptr;
	ln->next = l->tail;
	ln->prev = l->tail->prev;
	l->tail->prev = ln;
	ln->prev->next = ln;
	l->size++;
}

void zPtrListRemoveFirst(zPtrList* l){
	zPtrListNode* ln = l->head->next;
	if(ln == l->tail){
		return;
	}
	if(ln == l->current){
		l->current = ln->next;
	}
	l->head->next = ln->next;
	ln->next->prev = l->head;
	l->size--;
	zPtrListReleaseNode(ln);
}

void zPtrListRemoveLast(zPtrList* l){
	zPtrListNode* ln = l->tail->prev;
	if(ln == l->head){
		return;
	}
	if(ln == l->current){
		l->current = ln->prev;
	}
	l->tail->prev = ln->prev;
	ln->prev->next = l->tail;
	l->size--;
	zPtrListReleaseNode(ln);
}

void zPtrListRemoveCurrent(zPtrList* l){
	zPtrListNode* n = l->current;
	if(n == l->head || n == l->tail){
		return;
	}
	n->next->prev = n->prev;
	n->prev->next = n->next;
	l->current = n->prev;
	l->size--;
	zPtrListReleaseNode(n);
}


void zPtrListRemoveNext(zPtrList* l){
	zPtrListNode* n = l->current;
	if(l->current == l->tail){
		return;
	}
	l->current = l->current->next;
	zPtrListRemoveCurrent(l);
	l->current = n;
}


void zPtrListRemovePrev(zPtrList* l){
	zPtrListNode* n = l->current;
	if(l->current == l->head){
		return;
	}
	l->current = l->current->prev;
	zPtrListRemoveCurrent(l);
	l->current = n;
}

zPtrListPos zPtrListGetPos(zPtrList* l){
	return l->current;
}

void zPtrListSetPos(zPtrList* l,zPtrListPos p){
	l->current = p;
}

void zPtrListMovePosToNext(zPtrList* l,zPtrListPos p){
	zPtrListNode* n = p;
	zPtrListNode* c = l->current;

	if(c == NULL || c == n || n == NULL || n == l->head || n == l->tail){
		return;
	}
	n->prev->next = n->next;
	n->next->prev = n->prev;
	if(c == l->tail){
		c = c->prev;
	}
	n->next = c->next;
	n->prev = c;
	c->next->prev = n;
	c->next = n;
}

void zPtrListMovePosToPrev(zPtrList* l,zPtrListPos p){
	zPtrListNode* n = p;
	zPtrListNode* c = l->current;
	
	if(c == NULL || c == n || n == NULL || n == l->head || n == l->tail){
		return;
	}
	n->prev->next = n->next;
	n->next->prev = n->prev;
	if(c == l->head){
		c = c->next;
	}
	n->next = c;
	n->prev = c->prev;
	c->prev->next = n;
	c->prev = n;
}

void  zPtrListMoveCurrentToFirst(zPtrList* l){
	if(l->current == NULL || l->current == l->head || l->current == l->tail || 
		 l->current == l->head->next){
		return;
	}
	l->current->prev->next = l->current->next;
	l->current->next->prev = l->current->prev;
	l->head->next->prev = l->current;
	l->current->next = l->head->next;
	l->head->next = l->current;
	l->current->prev = l->head;
}

void  zPtrListMoveCurrentToLast(zPtrList* l){
	if(l->current == NULL || l->current == l->head || l->current == l->tail || 
		 l->current == l->tail->prev){
		return;
	}
	l->current->prev->next = l->current->next;
	l->current->next->prev = l->current->prev;
	l->tail->prev->next = l->current;
	l->current->prev = l->tail->prev;
	l->tail->prev = l->current;
	l->current->next = l->tail;
}

void* zPtrListGetFirst(zPtrList* l){
	if(l->head->next == l->tail){
		return NULL;
	}
	return l->head->next->data;
}

void* zPtrListGetLast(zPtrList* l){
	if(l->tail->prev == l->head){
		return NULL;
	}
	return l->tail->prev->data;
}

#endif

