/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
  zTools.h - part of the ZOE library for genomic analysis
  
  Copyright (C) 2001-2002 Ian F. Korf
  
  zTools provides some generally useful definitions, structures, and functions.
  All components of the ZOE library include zTools.
   
\******************************************************************************/

#ifndef ZOE_TOOLS_H
#define ZOE_TOOLS_H

#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAS_GLIB
#include <glib.h>
#endif

#define PADDING 50

typedef enum bool {false, true} bool;

/******************************************************************************\
 Library Information

Library build date.

\******************************************************************************/

void zLibInfo (void);

/******************************************************************************\
 Program Name

Every program that uses this library should register its program name. When a
program dies and reports an error message the program name is also reported.
This is especially important in sequence analysis pipelines where many programs
might be called.

\******************************************************************************/

void  zSetProgramName (char *string);
char *zGetProgramName (void);

/******************************************************************************\
 Commandline Processing

Commandline arguments are given as -flag or -attribute=value

\******************************************************************************/

void  zParseOptions (int*, char**);
void  zFreeOptions(); 
char* zOption (const char*);

/******************************************************************************\
 Important typedefs

These typedefs are all some form of int or float, but you should use these
rather than int for clarity. Each one has its own conversion function because
certain constants have string-based input/output formats but numeric internal
formats.

Description
===========

coor_t    for sequence coordinates, non-negative
frame_t   for DNA frame {0, 1, 2}, no negative values
score_t   for scores (log-odds ratios or log prob)
strand_t  for DNA strand {+, -, .}


Type      Special Value     String Representation
========  ================  =====================
coor_t    UNDEFINED_COOR     ...
frame_t   UNDEFINED_FRAME    .
score_t   MIN_SCORE          .
          MAX_SCORE          *
strand_t  UNDEFINED_STRAND   .

\******************************************************************************/

typedef unsigned int   coor_t;
typedef          char  frame_t;
typedef          double score_t;
typedef          char  strand_t;

extern const coor_t   UNDEFINED_COOR;
extern const frame_t  UNDEFINED_FRAME;
extern const score_t  MIN_SCORE;
extern const score_t  MAX_SCORE;
extern const strand_t UNDEFINED_STRAND;
void     zCoor2Text   (coor_t, char*);
coor_t   zText2Coor   (const char*);
void     zFrame2Text  (frame_t, char*);
frame_t  zText2Frame  (const char*);
void     zScore2Text  (score_t, char*);
void     zScore2FloatText  (score_t, char*);
score_t  zText2Score  (const char*);
void     zStrand2Text (strand_t, char*);
strand_t zText2Strand (const char*);

typedef enum {
	Phase0,
	Phase1,
	Phase2,
	Phase1T,  /* overhanging T */
	Phase2TA, /* overhanging TA */
	Phase2TG /* overhanging TG */
} zPhase_t;

zPhase_t zText2Phase(const char*);
void     zPhase2Text(const zPhase_t, char*);

/******************************************************************************\
 Error Messages

There are two functions for reporting errors. zWarn and zDie throw a message to
stderr and zDie terminates the program. The syntax is identical to printf.

\******************************************************************************/
void zWarn (const char*, ...);
void zDie  (const char*, ...);

/******************************************************************************\
 Tracing
 
The trace functions are for reporting the progress of the program as a
debugging tool. These are conditionally compiled by setting the TRACE level
(see below). zTrace1 is for high level messages, zTrace2 is for medium level
messages, and zTrace3 is for low level messages.

 #define TRACE 1
 gcc ... -DTRACE=1

\******************************************************************************/
void zTrace1 (const char*, ...);
void zTrace2 (const char*, ...);
void zTrace3 (const char*, ...);

/******************************************************************************\
 Memory Tools

The standard memory allocation tools have wrappers that take an additional
string argument that should be the name of the calling function and perhaps a
bit more information. If malloc, for example fails, it uses the string as an
error message and then terminates the program. This is one of the few places
where there are fatal error messages.

\******************************************************************************/
void* zMalloc (size_t, const char*);
void* zCalloc (size_t, size_t, const char*);
void* zRealloc (void*, size_t, const char*);
void  zFree (void*);

#ifdef __G_LIB_H__

gpointer zSliceAlloc (size_t, const char*); 
void  zSliceFree (size_t, gpointer); 

#else

void* zSliceAlloc(size_t size, const char *str);
void  zSliceFree (size_t size, void* p); 

#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

#endif /* __G_LIB_H__ */

/******************************************************************************\
 Comparison Functions

There are three main variable types: int, float, and text. Comparison operators
for each are provided with prototypes that will work with qsort and bsearch.

\******************************************************************************/
int zIcmp(const void*, const void*);
int zFcmp(const void*, const void*);
int zTcmp(const void*, const void*);
int zCoortCmp(const void*, const void*);

/******************************************************************************\
 Integer Vector

zIVec is a dynamic, insert-only array of integers. Before you can use zIVec,
you must initialize it with a default size. If you exceed this initial size,
the array will double in size (and continue doubling as needed).

	zIVec vector;
	zInitIVec(&vector, some_non_negative_number);

You addd elements with zPushIVec.

	zPushIVec(&vector, 20);
	zPushIVec(&vector, -1);

Arrays are often sorted. Use the stdlib qsort function to sort zIVec with the
zIcmp function.

	qsort(vector.elem, vector.size, sizeof(int), zIcmp);

You can iterate using standard array[] syntax using the attributes "size" and
"elem" as shown below.

	int i;
	for (i = 0; i < vector.size; i++) {
		printf("%d %d\n", i, vector.elem[i]);
	}

As a convenience, you can get the last element of the array with "last".

	printf("%d\n", vector.last);

When you're done, you should free the memory.

	zFreeIVec(&vector);

\******************************************************************************/
struct zIVec {
	int   *elem;
	int    size;
	int    limit;
	int    last;
};
typedef struct zIVec zIVec;
void zFreeIVec (zIVec*);
void zInitIVec (zIVec*, int);
void zPushIVec (zIVec*, int);

/******************************************************************************\
 Float Vector

zFVec is identical to zIVec except that it holds floats. See the zIVec
documentation above (note: use zFcmp for qsort).

\******************************************************************************/
struct zFVec {
	float *elem;
	int    size;
	int    limit;
	float  last;
};
typedef struct zFVec zFVec;
void zFreeFVec (zFVec*);
void zInitFVec (zFVec*, int);
void zPushFVec (zFVec*, float);

/******************************************************************************\
 Text Vector

zTVec is like zIVec except that it holds text (char*). zTVec copies each
string. If you do not want to copy the string, use the generic (void*) zVec.
For qsort, use the zTvec function.

\******************************************************************************/
struct zTVec {
	char **elem;
	int    size;
	int    limit;
	char  *last;
};
typedef struct zTVec zTVec;
void zFreeTVec (zTVec*);
void zInitTVec (zTVec*, int);
void zPushTVec (zTVec*, const char*);

/******************************************************************************\
 Generic Vector

zVec is a generic vector holding void* values. It works like the other vectors.

\******************************************************************************/
struct zVec {
	void **elem;
	int    size;
	int    limit;
	void  *last;
};
typedef struct zVec zVec;
void zInitVec (zVec*, int);
void zFreeVec (zVec*);
void zPushVec (zVec*, void*);
void zSortVec (zVec *vec, size_t size, int(*compar)(const void*, const void*)); 

/******************************************************************************\
 Generic Hash

zHash is an insert-only hash that provides lookup by string for void* values.
You must initialze the hash before using it.

	zHash hash;
	zInitHash(&hash);

To add values to the hash, you use "zSetHash" with a key-value pair. 
THE VALUE IS NOT COPIED, YOU MUST ALLOCATE STORAGE YOURSELF.

	zSetHash(&hash, key, value);
	zSetHash(&hash, "foo", value);

You should NOT do something like this:

	zSetHash(&hash, key, "bar");
	zSetHash(&hash, "foo", "bar");

You retrieve a value with "zGetHash". If there is no value associated with the
key, zGetHash returns NULL. In this way you can use zGetHash to test for the
existence of a hash key.

	val = zGetHash(&hash, key);
	if (val == NULL) do_something();
	else             do_otherthing();

If you want to retrieve all the keys or values of the hash, you use zKeysOfHash
or zValsOfHash. These return zVec pointers (that you must free after use) to
the keys and values. The keys and values will be unsorted.

	int i;
	void *val;
	zVec *keys;
	
	keys = zKeysOfHash(&hash);
	qsort(keys->elem, keys->size, sizeof(void*), zTcmp);
	for (i = 0; i < keys->size; i++) {
		val = zGetHash(&hash, keys.elem[i]);
	}

When you're done, you should ONLY free the memory used by the hash itself. 
The memory used by the keys in hash is freed automatically by the following function. 
The memory used by the values in the hash MUST BE FREED ELSEWHERE, WHEREVER THE POINTERS
WERE ALLOCATED BEFORE SENDING THEM INTO THE HASH.

	zFreeHash(&hash);

	(or)

	hash_ptr = zMalloc(sizeof(zHash),"something");
		...
	zFreeHash(hash_ptr);
	zFree(hash_ptr);

\******************************************************************************/
struct zHash {
	int   level;
	int   slots;
	int   keys;
	zVec *key;
	zVec *val;
};
typedef struct zHash zHash;

void  zInitHash (zHash*);
void  zFreeHash (zHash*);
int   zHashFunc (const zHash*, const char*);
void  zSetHash (zHash*, const char*, void*);
void* zGetHash (const zHash*, const char*);
zVec* zKeysOfHash (const zHash*);
zVec* zValsOfHash (const zHash*);
void  zStatHash (const zHash*);  /* to be deprecated */
void  zEStatHash (const zHash*); /* to be deprecated */

/******************************************************************************\
  String Pooling, Interning (Interenalization)

  String interning allows efficient comparison of strings using
  the == operator.  That it is, it converts strings to a canonical
  form, so that all char* point to the same object.
  
	char* s1 = "hello";
	char* s2;
	s2 = malloc(strlen(s1)+1);
	strcpy(s2, s1);
	
	assert(strcmp(s1,s2));
	assert(s1 != s2);
	
	s1 = zStrIntern(s1, 0);
	s2 = zStrIntern(s2, 1);
	
	assert(s1 == s2);
	assert(strcmp(s1,s2));
	
  zStrIntern()'s second parameter tells whether to free() the char*
  passed in.  A zero value indicates that the memory should not be
  freed.
\******************************************************************************/
typedef int zStrIdx; /* this has to be an int so that -1 can be used as an *
					  * impossible string index                            */

char*     zStrIntern (char* str, int freeStr);

/* Converts the string to an index into the pool, adding the string as needed */
zStrIdx   zChar2StrIdx(const char*);

/* returns 1 if the string is in the pool, 0 if not */
int       zStrIdxExists(const char*);

/* Gets the string for an index into the pool */
typedef const char* const_string_t;
const_string_t zStrIdx2Char(const zStrIdx);

/* Current number of strings in the pool */
int zStringPoolCount();

/* Free the global pool */
void zStringPoolFree(); 


/* Verbosity object */

void zSetVerbosityLevel(int);
int zGetVerbosityLevel();

void zAddVerbosityFeature(char*);
zTVec* zGetVerbosityFeatures();

typedef struct zVerbosityObject {

	int verbosity_level;
	zTVec features;

} zVerbosityObject;

void zInitVerbosityObject(zVerbosityObject* verbosity, int verbose);
void zFreeVerbosityObject(zVerbosityObject* verbosity); 
void zFreeVerbosityGlobalVariable(); 

/**********************************************************************\
  Linked List of objects
   this is a linked list object which allows all standard list operations
   the object takes functions on initialization which allows it to create
   (allocate) new objects (structs) stored in the list, delete those 
   obejcts and reset the objects (for reuse).  If you do not want the list
   to manage the creation/deletion of the objects (for example if you need 
   the object to persist after removed from the list) use the zPtrList 
   below. 
\**********************************************************************/

struct zListNode;
typedef struct zListNode zListNode;
struct zListNode {
	zListNode* next;
	zListNode* prev;
	void*      data;
};

typedef void* (*zListInitFunc)();
typedef void (*zListFreeFunc)(void*);
typedef void (*zListResetFunc)(void*);

typedef struct {
	zListNode* head;
	zListNode* tail;
	zListNode* current;
	zListNode* dead;
	int        size;
	zListInitFunc  init_func;
	zListFreeFunc  free_func;
	zListResetFunc reset_func;
} zList;

typedef zListNode* zListPos;

void  zInitList(zList*,zListInitFunc,zListFreeFunc,zListResetFunc);
void  zFreeList(zList*);
void  zResetList(zList*);
int   zListGetSize(zList*);
void* zListMoveFirst(zList*);
void* zListMoveLast(zList*);
void* zListMoveNext(zList*);
void* zListMovePrev(zList*);
void* zListGetCurrent(zList*);
bool  zListHasNext(zList*);
bool  zListHasPrev(zList*);
void* zListAddNext(zList*);
void* zListAddPrev(zList*);
void* zListAddFirst(zList*);
void* zListAddLast(zList*);
void  zListRemoveFirst(zList*);
void  zListRemoveLast(zList*);
void  zListRemoveCurrent(zList*);
zListPos zListGetPos(zList*);
void  zListSetPos(zList*,zListPos);
void  zListMovePosToNext(zList*,zListPos);
void  zListMovePosToPrev(zList*,zListPos);
void  zListMoveCurrentToFirst(zList*);
void  zListMoveCurrentToLast(zList*);
void* zListGetFirst(zList*);
void* zListGetLast(zList*);

/**********************************************************************\
  Linked List of pointers
    Just as the linked list above except that this list requires that 
    you allocate/free the memory of the objects stored in the list.
\**********************************************************************/

struct zPtrListNode;
typedef struct zPtrListNode zPtrListNode;
struct zPtrListNode {
	zPtrListNode* next;
	zPtrListNode* prev;
	void*         data;
#ifdef DEBUG
	int id;
#endif
};

typedef struct {
	zPtrListNode* head;
	zPtrListNode* tail;
	zPtrListNode* current;
	int        size;
#ifdef DEBUG
	int id;
#endif
} zPtrList;

typedef zPtrListNode* zPtrListPos;

void  zInitPtrList(zPtrList*);
void  zFreePtrList(zPtrList*);
void  zResetPtrList(zPtrList*);
int   zPtrListGetSize(zPtrList*);
void* zPtrListMoveFirst(zPtrList*);
void* zPtrListMoveLast(zPtrList*);
void* zPtrListMoveNext(zPtrList*);
void* zPtrListMovePrev(zPtrList*);
void* zPtrListGetCurrent(zPtrList*);
bool  zPtrListHasNext(zPtrList*);
bool  zPtrListHasPrev(zPtrList*);
void  zPtrListAddNext(zPtrList*,void*);
void  zPtrListAddPrev(zPtrList*,void*);
void  zPtrListAddFirst(zPtrList*,void*);
void  zPtrListAddLast(zPtrList*,void*);
void  zPtrListRemoveFirst(zPtrList*);
void  zPtrListRemoveLast(zPtrList*);
void  zPtrListRemoveCurrent(zPtrList*);
void  zPtrListRemovePrev(zPtrList*);
void  zPtrListRemoveNext(zPtrList*);
zPtrListPos zPtrListGetPos(zPtrList*);
void  zPtrListSetPos(zPtrList*,zPtrListPos);
void  zPtrListMovePosToNext(zPtrList*,zPtrListPos);
void  zPtrListMovePosToPrev(zPtrList*,zPtrListPos);
void  zPtrListMoveCurrentToFirst(zPtrList*);
void  zPtrListMoveCurrentToLast(zPtrList*);
void* zPtrListGetFirst(zPtrList*);
void* zPtrListGetLast(zPtrList*);

/**********************************************************************\
 snprintf - ANSI 
\**********************************************************************/

#ifdef __STRICT_ANSI__
extern int snprintf (char *__restrict __s, size_t __maxlen,
                     __const char *__restrict __format, ...)
     __THROW __attribute__ ((__format__ (__printf__, 3, 4)));
#endif /* __STRICT_ANSI__ */
#endif




