#ifndef NODE_H
#define NODE_H

#include <string.h>
#include <stdlib.h>
#include "cache.h"

#define STRLEN 1024

typedef struct node Node;

struct node {
  Node *parent;
  Node **children;
  int num_children;  
  int *descendants;
  int num_descendants;
  int num_states;
  int id;
  char seq[STRLEN];
  char name[STRLEN];

  float  **P; /* substitution probability matrix from parent to this node */
  double **E; /* expected number of substitutions from parent to this node */
  double **W; /* weight coefficients */

  float *pL_exact;   /* first pass */
  float *pL_partial; /* second pass */
  Cache *pL_exact_cache; /* caches first pass calculations */
  Cache *pL_partial_cache; /* caches second pass calculations */
};

Node* NewNode(int id);

#endif
