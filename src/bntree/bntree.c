#ifndef BNTREE_C
#define BNTREE_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "bntree.h"
#include "node.h"
#include "assert.h"

BNTree* NewBNTree(int cache_size) {
  int i;
  Node *root;

  BNTree *bntree = (BNTree *)malloc(sizeof (BNTree));

  bntree->cache_size = cache_size;

  /* set up the symbol map */
  /* predefined for now */
  for (i=0; i<NUM_CHARS; i++) {
    if (i == 'A')
      bntree->symbol_map[i] = 0;
    else if (i == 'C')
      bntree->symbol_map[i] = 1;
    else if (i == 'G')
      bntree->symbol_map[i] = 2;
    else if (i == 'T')
      bntree->symbol_map[i] = 3;
    else if (i == '_')
      bntree->symbol_map[i] = 4;
    else if (i == '.')
      bntree->symbol_map[i] = 5;
    else 
      bntree->symbol_map[i] = -1;
  }

  /* default name */
  strcpy (bntree->name, "Tree");

  /* set up the root */
  root = NewNode(0);
  strcpy (root->name, "root");
  strcpy (root->seq, "1");
  root->parent = NULL;
  bntree->root = root;
  bntree->num_nodes = 1;

  return bntree;
}

BNTree* NewBNTreeFromFile (char *filename, int cache_size) {
  return NewBNTreeFromFP (fopen(filename, "r"), NULL, -1, cache_size);
}

BNTree* NewBNTreeFromFP (FILE *fp, char *name, int order, int cache_size) {

  BNTree *bntree;
  char model_type[STRLEN];
  char treestr[STRLEN];
  char namestr[STRLEN];
  char parentstr[STRLEN];
  Node *node = NULL;
  int i, j; 
  int matrix_size;
  double val;
  int cpd_num;

  bntree = NewBNTree(cache_size);

  if (name == NULL) {  /* read header */
    fscanf (fp, "%s %s %d", bntree->name, model_type, &(bntree->order));
    if (strcmp(model_type, "BNTREE") != 0) {
      printf ("Error, model type %s, should be BNTREE\n", model_type);
      exit(-1);
    }
  }
  else { /* we're already past the header */
    strcpy (bntree->name, name);
    bntree->order = order;
  }
  
  bntree->num_states = intpow(ALPHABET_SIZE, bntree->order+1);
  bntree->root->num_states = bntree->num_states;

  /* read in tree topology */
  fscanf (fp, "%s", treestr);
  parse_tree_string (bntree, treestr);

  /* initialize the local probability distributions for each node */
  matrix_size = bntree->num_states;
  for (i=0; i<bntree->num_nodes; i++) {
    node = bntree->postorder[i];
    
    if (node->parent == NULL) continue;  /* no P for root */

    node->P = (float **)malloc(matrix_size * sizeof(float *));
    for (j=0; j<matrix_size; j++)
      node->P[j] = (float *)malloc(matrix_size * sizeof(float));
  }

  /* now read in all local probability distributions */
  for (cpd_num = 0; cpd_num < bntree->num_nodes-1; cpd_num++) {
    if (fscanf (fp, "%s --> %s", parentstr, namestr) != 2) {
      fprintf (stderr, "Error reading CPD header\n");
      exit(-1);
    }

    /* find the correct node */
    for (i=0; i<bntree->num_nodes; i++) {
      if (strcmp (namestr, bntree->postorder[i]->name) == 0) {
	node = bntree->postorder[i];
      }
    }
    
    /* read in P */
    for (i=0; i<matrix_size; i++) {
      for (j=0; j<matrix_size; j++) {
	fscanf (fp, "%lf", &val);
	node->P[i][j] = val;
      }
    }
  }    

  return bntree;
}

void parse_tree_string (BNTree *bntree, char *treestr) {

  char *current_seq = NULL;
  char *current_name = NULL;
  Node *node, *new_node;
  int current_id = 1;
  char c;
  int in_seq = 0;
  unsigned int j;
  int i;
  int index, ancestor_num;
  int key_length, array_length, entry_size, num_slots;

  node = bntree->root;

  for (j=0; j<strlen(treestr); j++) {
    c = treestr[j];

    if (in_seq) {
      if (isdigit(c)) {
	strncat (current_seq, &c, 1);
	continue;
      }
      else {
	in_seq = 0;
      }
    }

    if (c == '(' || c == '[') {
      new_node = NewNode(current_id++);
      new_node->num_states = bntree->num_states;
      node->children[node->num_children] = new_node;
      node->num_children++;
      new_node->parent = node;
      node = new_node;
      current_name = node->name;
      current_seq = node->seq;
      bntree->num_nodes++;
    }
    else if (c == ',') {
      new_node = NewNode(current_id++);
      new_node->num_states = bntree->num_states;
      node->parent->children[node->parent->num_children] = 
	new_node;
      node->parent->num_children++;
      new_node->parent = node->parent;
      node = new_node;
      current_name = node->name;
      current_seq = node->seq;
      bntree->num_nodes++; 
    }
    else if (c == ')' || c ==']') {
      node = node->parent;
      /* check for null right child */
      if (strcmp("none", node->children[1]->name) == 0) {
	node->num_children--;
	bntree->num_nodes--;
	current_id--;
	free (node->children[1]->children);
	free (node->children[1]);
	node->children[1] = NULL;
      }
      current_name = NULL;
    }
    else if (c == ':') {
      in_seq = 1;
    }
    else if (c == ';') {
      /* ignore terminating character */
    }
    else if (current_name != NULL) {
      if (!isspace(c) || strlen(current_name) > 0) {
	strncat(current_name, &c, 1);
      }
    }
  }

  /* create the postorder traversal list */
  /* also figure out the number of sequences here */
  bntree->num_seqs = 0;
  bntree->postorder = malloc(bntree->num_nodes * sizeof (Node *));
  index = 0;
  create_postorder (bntree, bntree->root, &index);

  /* now create the preorder traversal list */
  /* also set the ancestor node names here */
  bntree->preorder = malloc(bntree->num_nodes * sizeof (Node *));
  index = 0;
  ancestor_num = 1;
  create_preorder (bntree, bntree->root, &index, &ancestor_num);  

  /* for each node, figure out which sequences descend from it */
  for (i=0; i<bntree->num_nodes; i++) {
    node = bntree->preorder[i];
    node->descendants = (int *) calloc (bntree->num_seqs, sizeof (int));
    node->num_descendants = 0;
    get_descendants (node, node);
  }

  /* set up the caches for all non-leaf nodes */
  for (i=0; i<bntree->num_nodes; i++) {
    if (bntree->cache_size != 0) {
      node = bntree->preorder[i];
      if (node->num_children > 0) {
	/* first pass cache */
	/* figure out the key length and array length */
	key_length = ceil (1.0/8 * log (longpow (longpow (ALPHABET_SIZE, bntree->order+1), 
						 node->num_descendants)) / log (2));

	if (node->parent == NULL) {
	  /* root node, separate entry for each pL_exact */
	  array_length = 1;
	}
	else {
	  array_length = bntree->num_states;
	}
	entry_size = sizeof (CacheEntry) + key_length + array_length * sizeof (float);
	num_slots = bntree->cache_size / entry_size;
	node->pL_exact_cache = NewCache (num_slots, bntree->cache_size, key_length, array_length);
	
	/* second pass cache */
	if (bntree->order > 0) {
	  /* figure out the key length and array length */
	  key_length = ceil (1.0/8 * log (longpow (longpow (ALPHABET_SIZE, bntree->order+1), 
						node->num_descendants)) / log (2));
	  if (node->parent == NULL) {
	    /* root node, separate entry for each pL_exact */
	    array_length = 1;
	  }
	  else {
	    array_length = bntree->num_states;
	  }
	  entry_size = sizeof (CacheEntry) + key_length + array_length * sizeof (float);
	  num_slots = bntree->cache_size / entry_size;
	  node->pL_partial_cache = NewCache (num_slots, bntree->cache_size, key_length, array_length);
	}			       
      }
      else {
	/* no cache for this particular node */
	node->pL_exact_cache = NULL;
	node->pL_partial_cache = NULL;
      }
    }
    else {
      /* no cacheing */
      node->pL_exact_cache = NULL;
      node->pL_partial_cache = NULL;
    }
  }
}

void get_descendants (Node *ancestor_node, Node *current_node) {
  int seq_num;
  int i;

  if (current_node->num_children == 0) {
    seq_num = atoi (current_node->seq) - 1;
    ancestor_node->descendants[seq_num] = 1;
    ancestor_node->num_descendants++;
  }
  else if (current_node->parent == NULL) {
    /* For convenience, consider the root node to have descended from itself */
    seq_num = atoi (current_node->seq) - 1;
    ancestor_node->descendants[seq_num] = 1;
    ancestor_node->num_descendants++;
    for (i=0; i<current_node->num_children; i++) {
      get_descendants (ancestor_node, current_node->children[i]);
    }
  }
  else { /* internal node */
    for (i=0; i<current_node->num_children; i++) {
      get_descendants (ancestor_node, current_node->children[i]);
    }
  }
}

void print_tree (FILE *fp, BNTree *bntree) {
  int u;
  int i, j;
  Node *node;

  /* header */
  fprintf (fp, "%s\tBNTREE\t%d\n", bntree->name, bntree->order);

  /* tree topology */
  print_subtree (fp, bntree->root);
  fprintf (fp, "\n");

  /* CPDs */
  for (u=0; u<bntree->num_nodes; u++) {
    node = bntree->preorder[u];

    if (node->parent == NULL) continue;

    fprintf (fp, "%s --> %s\n", node->parent->name, node->name);
    for (i=0; i<bntree->num_states; i++) {
      for (j=0; j<bntree->num_states; j++) {
	fprintf (fp, "%g ", node->P[i][j]);
      }
      fprintf (fp, "\n");
    }
  }
}

void print_subtree (FILE *fp, Node *node) {
  
  if (node->num_children == 0)
    fprintf (fp, "%s:%s", node->name, node->seq);
  else {
    fprintf (fp, "(");
    print_subtree (fp, node->children[0]);
    fprintf (fp, ",");
    if (node->children[1] != NULL)
      print_subtree (fp, node->children[1]);
    else
      fprintf (fp, "none");
    fprintf (fp, ")");
  }
}

/* creates postorder and also figures out total number of sequences */
void create_postorder (BNTree *bntree, Node *node, int *index) {
  int i;

  if (node->parent == NULL || node->num_children == 0)
    bntree->num_seqs++;

  for (i=0; i<node->num_children; i++) {
    create_postorder (bntree, node->children[i], index);
  }
  
  bntree->postorder[*index] = node;
  (*index)++;
}

void create_preorder (BNTree *bntree, Node *node, int *index, int *ancestor_num) {
  int i;
  char ancestorstr[STRLEN];

  if (node->parent != NULL && node->num_children != 0) {
    /* ancestor node */
    sprintf (ancestorstr, "a%d", *ancestor_num);
    strcpy (node->name, ancestorstr); 
    (*ancestor_num)++;
  }

  bntree->preorder[*index] = node;
  (*index)++;

  for (i=0; i<node->num_children; i++) {
    create_preorder (bntree, node->children[i], index, ancestor_num);
  }  
}

/* don't call with negative arguments */
int intpow (int x, int y) {
  int i;
  int result = 1;
  
  for (i=0; i<y; i++)
    result *= x;

  return result;
}

unsigned long longpow (unsigned long x, unsigned long y) {
  unsigned long i;
  unsigned long result = 1;

  for (i=0; i<y; i++)
    result *= x;
  
  return result;
}

double calc_tuple_str_prob (BNTree *bntree, char *tuple, int conditional) {
  int *assignment = malloc (bntree->num_seqs * sizeof(int));
  double result;
  int tuple_idx;
  int i, j;

  if ((int)strlen(tuple) != (bntree->order + 1) * (bntree->num_seqs + 1) - 1) {
    fprintf (stderr, "Wrong length for tuple %s\n", tuple);
    exit (-1);
  }

  /* convert the tuple string into an assignment array */
  for (i=0; i<bntree->num_seqs; i++)
    assignment[i] = 0;

  tuple_idx = 0;
  for (i=bntree->order; i>=0; i--) {
    for (j=0; j<bntree->num_seqs; j++) {
      assignment[j] += intpow (ALPHABET_SIZE, i) * bntree->symbol_map[(int)tuple[tuple_idx]];
      tuple_idx++;
    }
    tuple_idx++;
  }
  
  result = calc_tuple_prob (bntree, assignment, conditional);
  free (assignment);
  return result;
}

double calc_tuple_prob (BNTree *bntree, int *assignment, int conditional) {
  int num_passes;

  if (bntree->order > 0 && conditional)
    num_passes = 2;
  else
    num_passes = 1;

  /* new method */
  if (num_passes == 1) {
    return calc_pL (bntree->root, assignment, assignment[0], 0);
  }
  else {
    return calc_pL (bntree->root, assignment, assignment[0], 0) /
      calc_pL (bntree->root, assignment, assignment[0], 1);
  }
}

/* This is a bad way to do it, and won't work for key length > 8 */
/* Should be rewritten */
char* assignment_to_key (Node *node, int *assignment, int pass) {
  unsigned long base10;
  int i;
  int j;
  int key_length;
  char *key;

  /* first convert the assignment of this node's descendants to
     a base-10 number */
  base10 = 0;
  i=0;
  j=0;
  if (pass == 1 && node->parent == NULL) {
    /* special case, partial match for everything but root */
    base10 += assignment[0];
    for (i=1; i<node->num_descendants; i++)
      base10 += longpow (node->num_states / ALPHABET_SIZE, i+1) *
	(assignment[i] - assignment[i] % ALPHABET_SIZE);
  }
  else {
    while (i < node->num_descendants) {
      if (node->descendants[j]) {
	if (pass == 0)
	    base10 += longpow (node->num_states, i) * assignment[j]; 
	else
	  base10 += longpow (node->num_states / ALPHABET_SIZE, i) *
	    (assignment[j] - assignment[j] % ALPHABET_SIZE);
	i++;
      }
      j++;
    }
  }

  /* printf ("key is %lld\n", base10); */

  /* now convert the base-10 number to a byte array (base 256) */
  if (pass == 0)
    key_length = node->pL_exact_cache->key_length;
  else
    key_length = node->pL_partial_cache->key_length;

  key = (char *) malloc (key_length);

  for (i=key_length-1; i>=0; i--) {
    key[i] = base10 / longpow (256, i);
    base10 -= base10 / longpow (256, i);
  }

  assert (base10 == 0);

  return key;
}

void calc_pLs (Node *node, int *assignment, int pass) {
  int i;
  int j;
  int k;
  char *key = NULL;
  float *vals = NULL;
  int seq_num;
  int leaf_assignment;
  float sum;

  if (node->num_children == 0) {
    /* this node is a leaf, base case */
    seq_num = atoi (node->seq) - 1;

    if (pass == 0) {
      leaf_assignment = assignment[seq_num];
      node->pL_exact = (float *) malloc (node->num_states * sizeof (float));
      for (i=0; i<node->num_states; i++) {
	if (i == leaf_assignment)
	  node->pL_exact[i] = 1;
	else
	  node->pL_exact[i] = 0;
      }
    }
    else if (pass == 1) {
      /* we have a partial assignment only */
      leaf_assignment = assignment[seq_num] - assignment[seq_num] % ALPHABET_SIZE;    
      node->pL_partial = (float *) malloc (node->num_states * sizeof (float));
      for (i=0; i<node->num_states; i++) {
	if (i - i % ALPHABET_SIZE == leaf_assignment)
	  node->pL_partial[i] = 1;
	else
	  node->pL_partial[i] = 0;
      }
    }

    return;
  }

  /* non-leaf node, recursive case */
  /* first look in cache to see if we can remember the answer */
  if (pass == 0) {
    if (node->pL_exact_cache != NULL) {
      key = assignment_to_key (node, assignment, pass);
      vals = cache_get (node->pL_exact_cache, key);
    }
    if (vals != NULL) {
      /* cache hit, no need to calculate anything */
      node->pL_exact = vals;
      free (key);
      return;
    }
  }
  else if (pass == 1) {
    if (node->pL_partial_cache != NULL) {
      key = assignment_to_key (node, assignment, pass);
      vals = cache_get (node->pL_partial_cache, key);
    }
    if (vals != NULL) {
      /* cache hit, no need to calculate anything */
      node->pL_partial = vals;
      free (key);
      return;
    }
  }

  /* cache miss, we have to do some calculations */

  /* first calculate the appropriate pLs for this node's children */
  for (i=0; i<node->num_children; i++) {
    calc_pLs (node->children[i], assignment, pass);
  }

  /* now use that information to calculate this node's pLs */
  if (pass == 0) {
    node->pL_exact = (float *) malloc (node->num_states * sizeof (float));
    for (i=0; i<node->num_states; i++) { /* this node's value */
      node->pL_exact[i] = 1;
      for (j=0; j<node->num_children; j++) { /* child number */
	sum = 0;
	for (k=0; k<node->num_states; k++) { /* child's value */
	  sum += node->children[j]->P[i][k] * node->children[j]->pL_exact[k];
	}
	node->pL_exact[i] *= sum;
      }
    }

    /* store results in cache */
    if (node->pL_exact_cache != NULL) {
      cache_insert (node->pL_exact_cache, key, node->pL_exact);
      free (key);
    }
    
    /* free children's pLs */
    for (i=0; i<node->num_children; i++) {
      if (node->children[i]->pL_exact_cache == NULL ||
	  node->children[i]->num_children == 0) {
	/* if we're not cacheing or child is a leaf, free up memory */
	free (node->children[i]->pL_exact);
      }
    }
  }
  else {
    assert (pass == 1);
    node->pL_partial = (float *) malloc (node->num_states * sizeof (float));
    for (i=0; i<node->num_states; i++) { /* this node's value */
      node->pL_partial[i] = 1;
      for (j=0; j<node->num_children; j++) { /* child number */
	sum = 0;
	for (k=0; k<node->children[j]->num_states; k++) { /* child's value */
	  sum += node->children[j]->P[i][k] * node->children[j]->pL_partial[k];
	}
	node->pL_partial[i] *= sum;
      }
    }

    /* store results in cache */
    if (node->pL_partial_cache != NULL) {
      cache_insert (node->pL_partial_cache, key, node->pL_partial);
      free (key);
    }

    /* free children's pLs */
    for (i=0; i<node->num_children; i++) {
      if (node->children[i]->pL_partial_cache == NULL ||
	  node->children[i]->num_children == 0) {
	/* if we're not cacheing or child is a leaf, free up memory */
	free (node->children[i]->pL_partial);
      }
    }
  }
}

float calc_pL (Node *node, int *assignment, int node_value, int pass) {
  int i;
  int j;
  int k;
  char *key = NULL;
  float *vals = NULL;
  float pL;
  float sum;
  
  /* first look in cache to see if we can remember the answer */
  if (pass == 0) {
    if (node->pL_exact_cache != NULL) {
      key = assignment_to_key (node, assignment, pass);
      vals = cache_get (node->pL_exact_cache, key);
    }
    if (vals != NULL) {
      /* cache hit, no need to calculate anything */
      free (key);
      return *vals;
    }
  }
  else if (pass == 1) {
    if (node->pL_partial_cache != NULL) {
      key = assignment_to_key (node, assignment, pass);
      vals = cache_get (node->pL_partial_cache, key);
    }
    if (vals != NULL) {
      /* cache hit, no need to calculate anything */
      free (key);
      return *vals;
    }
  }
  
  /* cache miss, time to calculate */
  /* first calculate the appropriate pLs for this node's children */
  for (i=0; i<node->num_children; i++) {
    calc_pLs (node->children[i], assignment, pass);
  }

  /* now use that information to calculate this node's specific pL */
  if (pass == 0) {
    pL = 1;
    for (j=0; j<node->num_children; j++) { /* child number */
      sum = 0;
      for (k=0; k<node->children[j]->num_states; k++) { /* child's value */
	sum += node->children[j]->P[node_value][k] * node->children[j]->pL_exact[k];
      }
      pL *= sum;
    }

    /* store results in cache */
    if (node->pL_exact_cache != NULL) {
      vals = (float *) malloc (sizeof (float));
      vals[0] = pL;
      cache_insert (node->pL_exact_cache, key, vals);
      free (key);
    }

    /* free children's pLs */
    for (i=0; i<node->num_children; i++) {
      if (node->children[i]->pL_exact_cache == NULL ||
	  node->children[i]->num_children == 0) {
	/* if we're not cacheing or child is a leaf, free up memory */
	free (node->children[i]->pL_exact);
      }
    }

    return pL;
  }
  else {
    assert (pass == 1);
    pL = 1;
    for (j=0; j<node->num_children; j++) { /* child number */
      sum = 0;
      for (k=0; k<node->children[j]->num_states; k++) { /* child's value */
	sum += node->children[j]->P[node_value][k] * node->children[j]->pL_partial[k];
      }
      pL *= sum;
    }

    /* store results in cache */
    if (node->pL_partial_cache != NULL) {
      vals = (float *) malloc (sizeof (float));
      vals[0] = pL;
      cache_insert (node->pL_partial_cache, key, vals);
      free (key);
    }

    /* free children's pLs */
    for (i=0; i<node->num_children; i++) {
      if (node->children[i]->pL_partial_cache == NULL ||
	  node->children[i]->num_children == 0) {
	/* if we're not cacheing or child is a leaf, free up memory */
	free (node->children[i]->pL_partial);
      }
    }

    return pL;
  }
}

void free_all_caches (BNTree *bntree) {
  int i;
  Node *node;

  for (i=0; i<bntree->num_nodes; i++) {
    node = bntree->preorder[i];
    if (node->pL_exact_cache != NULL)
      cache_free (node->pL_exact_cache);
    if (node->pL_partial_cache != NULL)
      cache_free (node->pL_partial_cache);
  }
}

void calc_expected_num_subst (BNTree *bntree, char **tuple_arr,
			      double *count_arr, int num_tuples) {
  int i, j, k;
  int u, v;
  int t;
  int tuple_idx;
  Node *node;
  double **pL;
  double **pU;
  int *assignment;

  double tuple_prob;
  double denom;
  double subst_prob;

  pL = (double **) malloc (bntree->num_nodes * sizeof(double *));
  for (i=0; i<bntree->num_nodes; i++)
    pL[i] = (double *) malloc (bntree->num_states * sizeof (double));
  
  pU = (double **) malloc (bntree->num_nodes * sizeof(double *));
  for (i=0; i<bntree->num_nodes; i++)
    pU[i] = (double *) malloc (bntree->num_states * sizeof (double));

  assignment = (int *) malloc (bntree->num_seqs * sizeof (int));

  for (u = 0; u < bntree->num_nodes; u++) {
    node = bntree->preorder[u];
    if (node->parent == NULL) continue;

    for (i = 0; i < bntree->num_states; i++) {
      for (j = 0; j < bntree->num_states; j++) {
	node->E[i][j] = 0;
      }
    }
  }

  for (t = 0; t < num_tuples; t++) {

    char *tuple = tuple_arr[t];

    /* convert the tuple string into an assignment array */
    for (i=0; i<bntree->num_seqs; i++)
      assignment[i] = 0;

    tuple_idx = 0;
    for (i=bntree->order; i>=0; i--) {
      for (j=0; j<bntree->num_seqs; j++) {
	assignment[j] += intpow (ALPHABET_SIZE, i) * bntree->symbol_map[(int)tuple[tuple_idx]];
	tuple_idx++;
      }
      tuple_idx++;
    }

    /* first traverse in postorder and calculate pL */
    for (u=0; u<bntree->num_nodes; u++) {
      Node *node = bntree->postorder[u];
      
      if (node->num_children == 0) { /* leaf node, base case */
	int seq_num = atoi(node->seq) - 1;
	for (i=0; i<bntree->num_states; i++) {
	  if (assignment[seq_num] == i)
	    pL[node->id][i] = 1;
	  else
	    pL[node->id][i] = 0;
	}
      }

      else if (node->parent == NULL) { /* root node, observed */
	for (i=0; i<bntree->num_states; i++) {
	  if (assignment[0] == i) {
	    pL[node->id][i] = 1;
	    for (v=0; v<node->num_children; v++) {
	      Node *child = node->children[v];
	      double child_total = 0;
	      for (j=0; j<bntree->num_states; j++) {
		child_total += pL[child->id][j] * child->P[i][j];
	      }
	      pL[node->id][i] *= child_total;
	    }
	  }
	  else {
	    pL[node->id][i] = 0;
	  }
	}
      }

      else { /* internal, non-root node */
	for (i=0; i<bntree->num_states; i++) {
	  pL[node->id][i] = 1;
	  for (v=0; v<node->num_children; v++) {
	    Node *child = node->children[v];
	    double child_total = 0;
	    for (j=0; j<bntree->num_states; j++) {
	      child_total += pL[child->id][j] * child->P[i][j];
	    }
	    pL[node->id][i] *= child_total;
	  }
	}
      }
    }

    /* now traverse in preorder and calculate pU */
    for (u = 0; u < bntree->num_nodes; u++) {
      node = bntree->preorder[u];
      
      if (node->parent == NULL) { /* root, base case */
	for (i=0; i<bntree->num_states; i++) {
	  if (assignment[0] == i)
	    pU[node->id][i] = 1;
	  else
	    pU[node->id][i] = 0;
	}
      }
      else {
	/* For now we are assuming max 1 sibling...otherwise
	   this does not work! */
	assert (node->parent->num_children <= 2);
	
	if (node->parent->num_children == 1) { /* no sibling */
	  for (i=0; i<bntree->num_states; i++) {
	    pU[node->id][i] = 0;
	    for (j=0; j<bntree->num_states; j++) {
	      pU[node->id][i] += pU[node->parent->id][j] * node->P[j][i];
	    }
	  }
	}
 
	else if (node->parent->num_children == 2) {
	  Node *sibling = node->parent->children[0] == node ?
	    node->parent->children[1] : node->parent->children[0];
	  double *tmp = (double *) malloc (bntree->num_states * sizeof (double));
	  for (j=0; j<bntree->num_states; j++) {
	    tmp[j] = 0;
	    for (k=0; k<bntree->num_states; k++) {
	      tmp[j] += pU[node->parent->id][j] * pL[sibling->id][k] *
		sibling->P[j][k];
	    }
	  }
	  
	  for (i=0; i<bntree->num_states; i++) {
	    pU[node->id][i] = 0;
	    for (j=0; j<bntree->num_states; j++) {
	      pU[node->id][i] += tmp[j] * node->P[j][i];
	    }
	  }
	  free (tmp);
	}
      }
    }

    for (u=0; u<bntree->num_nodes; u++) {
      node = bntree->preorder[u];

      if (node->parent == NULL) continue;

      tuple_prob = 0;
      for (i=0; i<bntree->num_states; i++)
	tuple_prob += pL[node->id][i] * pU[node->id][i];
      if (tuple_prob == 0) continue; /* is this right? */

      for (i=0; i<bntree->num_states; i++) {
	denom = 0;
	for (j=0; j<bntree->num_states; j++)
	  denom += node->P[i][j] * pL[node->id][j];
	if (denom == 0) continue; /* is this right? */

	for (j=0; j<bntree->num_states; j++) {
	
	  /* probability of parent being i given the observed tuple */	  
	  subst_prob = pL[node->parent->id][i] * pU[node->parent->id][i] /
	    tuple_prob;

	  /* and this node being j given the observed tuple */
	  subst_prob *= node->P[i][j] * pL[node->id][j] /
	    denom;

	  node->E[i][j] += count_arr[t] * subst_prob;
	}
      }
    }
  }

  for (i=0; i<bntree->num_nodes; i++) {
    free (pL[i]);
    free (pU[i]);
  }
  free (pL);
  free (pU);
  free (assignment);
}

#endif
