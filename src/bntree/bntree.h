#ifndef BNTREE_H
#define BNTREE_H

#include <stdio.h>
#include "node.h"
#include "celllistitem.h"

#define ALPHABET_SIZE 6
#define NUM_CHARS 256

typedef struct bntree BNTree;

struct bntree {
  Node *root;
  char name[STRLEN];
  int num_states;
  int order;
  int num_nodes;
  int num_seqs;
  int symbol_map[NUM_CHARS];
  int cache_size;

  Node *current_node; /* for training */
  int num_params;
  int **param_map; /* map of cell row and col to param index */
  CellListItem **inv_param_map; /* map of param numbers to cells */
  int **weight_idx; /* where to look for the weight coefficients given cell */

  Node **preorder;  /* array of node pointers in preorder traversal order */
  Node **postorder; /* array of node pointers in postorder traversal order */
};

BNTree* NewBNTree(int cache_size);
BNTree* NewBNTreeFromFile (char *fname, int cache_size);
BNTree* NewBNTreeFromFP (FILE *fp, char *name, int order, int cache_size);
void parse_tree_string (BNTree *bntree, char *treestr);

double calc_tuple_int_prob (BNTree *bntree, int tuple, int conditional);
double calc_tuple_str_prob (BNTree *bntree, char *tuple, int conditional);
double calc_tuple_prob (BNTree *bntree, int *assignment, int conditional);

/* helper functions for calc_tuple_prob */
void calc_pLs (Node *node, int *assignment, int pass);
float calc_pL (Node *node, int *assignment, int node_value, int pass);

/* for training */
void calc_expected_num_subst (BNTree *bntree, char **tuple_arr,
			      double *count_arr, int num_tuples);

void print_tree (FILE *fp, BNTree *bntree);
void free_all_caches (BNTree *bntree);

/* misc. helper functions */
void print_subtree (FILE *fp, Node *node);
void create_postorder (BNTree *bntree, Node *node, int *index);
void create_preorder (BNTree *bntree, Node *node, int *index, int *ancestor_num);
void get_descendants (Node *ancestor_node, Node *current_node);
int intpow (int x, int y);
unsigned long longpow (unsigned long x, unsigned long y);

#endif
