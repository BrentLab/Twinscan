#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "training.h"
#include "bntree.h"
#include "math.h"


#define NEG_INF -1e20
#define PSEUDOCOUNT 0.01  /* just to prevent 0 counts */
#define EM_TOL 1e-5
#define MAX_EM_ITERS 1000

BNTree* train_bntree (char *treestr, char *modelstr, char *ss_filename, char *model_name, char *output_filename) {
  FILE *ss_fp, *out_fp;
  int model_type;
  double *background_probs;
  char **tuple_arr;
  double *count_arr;
  int num_tuples = 0;
  double ll;
  double last_ll;
  int iteration;
  int converged;
  char tag_str[STRLEN];
  char val_str[STRLEN];
  int u;
  int i, j, k;
  Node *node;
  int param_num;
  double row_sum;
  double param_val;
  CellListItem *cli;

  BNTree *bntree = NewBNTree(0); /* no cache for training (for now) */

  strcpy (bntree->name, model_name);

  if (strcmp (modelstr, "TT0") == 0) {
    model_type = TT0;
    bntree->order = 0;
    bntree->num_states = 6;
  }
  else if (strcmp (modelstr, "TT1") == 0) {
    model_type = TT1;
    bntree->order = 1;
    bntree->num_states = 36;
  }
  else if (strcmp (modelstr, "TT2") == 0) {
    model_type = TT2;
    bntree->order = 2;
    bntree->num_states = 216;
  }
  else if (strcmp (modelstr, "R0") == 0) {
    model_type = R0;
    bntree->order = 0;
    bntree->num_states = 6;
    bntree->num_params = 18;
  }
  else if (strcmp (modelstr, "R1") == 0) {
    model_type = R1;
    bntree->order = 1;
    bntree->num_states = 36;
    bntree->num_params = 236;
  }
  else if (strcmp (modelstr, "R2") == 0) {
    model_type = R2;
    bntree->order = 2;
    bntree->num_states = 216;
  }
  else if (strcmp (modelstr, "G0") == 0) {
    model_type = G0;
    bntree->order = 0;
    bntree->num_states = 6;
    bntree->num_params = 36;
  }
  else if (strcmp (modelstr, "G1") == 0) {
    model_type = G1;
    bntree->order = 1;
    bntree->num_states = 36;
    bntree->num_params = 1296;
  }
  else if (strcmp (modelstr, "G2") == 0) {
    model_type = G2;
    bntree->order = 2;
    bntree->num_states = 216;
    bntree->num_params = 46656;
  }
  else {
    fprintf (stderr, "Unknown probability model type: %s\n", modelstr);
    exit(-1);
  }
    
  /* initialize tree structure */
  parse_tree_string (bntree, treestr);

  /* parse the sufficient statistics file */
  ss_fp = fopen (ss_filename, "r");
  while (fscanf (ss_fp, "%s = %s\n", tag_str, val_str) == 2) {
    if (strcmp (tag_str, "NTUPLES") == 0)
	num_tuples = atoi(val_str);
  }
  fclose (ss_fp);

  tuple_arr = (char **) malloc (num_tuples * sizeof (char *));
  for (i=0; i<num_tuples; i++)
    tuple_arr[i] = (char *) malloc (STRLEN * sizeof (char));

  count_arr = (double *) malloc (num_tuples * sizeof (double));
  background_probs = (double *) malloc (bntree->num_states * sizeof (double));

  parse_ss (bntree, ss_filename, tuple_arr, count_arr, background_probs);


  /* initialize parameters */
  bntree->param_map = (int **)malloc(bntree->num_states * sizeof(int *));
  for (i=0; i<bntree->num_states; i++)
    bntree->param_map[i] = (int *)malloc(bntree->num_states * sizeof(int));

  bntree->weight_idx = (int **)malloc(bntree->num_states * sizeof(int *));
  for (i=0; i<bntree->num_states; i++)
    bntree->weight_idx[i] = (int *)malloc(bntree->num_states * sizeof(int));

  if (model_type == R0)
    init_params_R0 (bntree);
  else if (model_type == G0)
    init_params_G0 (bntree);

  else if (model_type == R1)
    init_params_R1 (bntree);
  else if (model_type == G1)
    init_params_G1 (bntree);

  else if (model_type == G2)
    init_params_G2 (bntree);
  else {
    fprintf (stderr, "Model type %s not yet supported\n", modelstr);
    exit (-1);
  }


  /* create inverse parameter map */
  bntree->inv_param_map = (CellListItem **)malloc(bntree->num_params * sizeof (CellListItem *));
  for (i=0; i<bntree->num_params; i++)
    bntree->inv_param_map[i] = NewCellList();
  
  for (i=0; i<bntree->num_states; i++) {
    for (j=0; j<bntree->num_states; j++) {
      param_num = bntree->param_map[i][j];
      cell_list_append(bntree->inv_param_map[param_num], i, j);
    }
  }


  /* allocate local CPDs and expectations */
  for (u=0; u<bntree->num_nodes; u++) {
    node = bntree->preorder[u];
    if (node->parent == NULL) continue;
    
    node->E = (double **)malloc(bntree->num_states * sizeof(double *));
    for (i=0; i<bntree->num_states; i++)
      node->E[i] = (double *)malloc(bntree->num_states * sizeof(double));

    node->P = (float **)malloc(bntree->num_states * sizeof(float *));
    for (i=0; i<bntree->num_states; i++)
      node->P[i] = (float *)malloc(bntree->num_states * sizeof(float));

    node->W = (double **)malloc(bntree->num_states * sizeof(double *));
    for (i=0; i<bntree->num_states; i++)
      node->W[i] = (double *)malloc(bntree->num_states * sizeof(double));
  }

  /* randomly initialize the parameters in a way that guarantees all 
     probabilities will be within 0 and 1, and all rows sum to 1 */
  /* also, self substitution probs will be at least 0.5 */
  for (u=0; u<bntree->num_nodes; u++) {
    node = bntree->preorder[u];
    if (node->parent == NULL) continue;

    srand(time(NULL));
    for (k=0; k<bntree->num_params - bntree->num_states; k++) {
      param_val = bounded_random (0.05, 0.50) / (bntree->num_states - 1);
      for (cli = bntree->inv_param_map[k]; cli != NULL; cli = cli->next) {
	i = cli->row;
	j = cli->col;
	node->P[i][j] = param_val;
      }
    }
    for (i=0; i<bntree->num_states; i++) {
      row_sum = 0;
      for (j=0; j<bntree->num_states; j++) {
	if (j != i) {
	  row_sum += node->P[i][j];
	}
      }
      node->P[i][i] = 1 - row_sum;
    }
  }

  /* run EM algorithm */
  /* We'll output the current tree after every iteration */
  iteration = 0;
  converged = 0;
  ll = calc_log_likelihood (bntree, tuple_arr, count_arr, num_tuples, 0);
  fprintf (stderr, "Initial log likelihood: %f\n", ll);
  while (! converged) {
    iteration++;

    /* improve likelihood */
    Estep(bntree, tuple_arr, count_arr, num_tuples);
    Mstep(bntree);

    /* output new model */
    out_fp = fopen (output_filename, "w");
    print_tree (out_fp, bntree);
    fclose (out_fp);

    /* check for convergence */
    last_ll = ll;
    ll = calc_log_likelihood (bntree, tuple_arr, count_arr, num_tuples, 0);
    fprintf (stderr, "EM iteration %d, log likelihood: %f\n", iteration, ll);
    if ( (last_ll - ll) / last_ll < EM_TOL || iteration >= MAX_EM_ITERS)
      converged = 1;
  }
  
  return bntree;
}

void Estep (BNTree *bntree, char **tuple_arr, double *count_arr, int num_tuples) {
  int u;
  int i, j;
  Node *node;

  /* first calculate the straight-up expectations */
  calc_expected_num_subst (bntree, tuple_arr, count_arr, num_tuples);

  /* add pseudocounts */
  for (u=0; u<bntree->num_nodes; u++) {
    node = bntree->preorder[u];
    if (node->parent == NULL) continue;
    for (i=0; i<bntree->num_states; i++) {
      for (j=0; j<bntree->num_states; j++) {
	node->E[i][j] += PSEUDOCOUNT;
      }
    }
  }
}

void Mstep (BNTree *bntree) {
  int i, j;
  int u;
  Node *node;
  double row_sum;
  double *row_sums;
  double count_sum;
  CellListItem *cli;
  int k;
  double weight_sum;
  double prob_sum;

  /* for general models, the M-step is trivial */
  /*
  if (model_type == G0 || model_type == G1 || model_type == G2) {
    for (u=0; u<bntree->num_nodes; u++) {
      node = bntree->preorder[u];
      if (node->parent == NULL) continue;

      for (i=0; i<bntree->num_states; i++) {
	row_sum = 0;
	for (j=0; j<bntree->num_states; j++)
	  row_sum += node->E[i][j];
	for (j=0; j<bntree->num_states; j++)
	  node->P[i][j] = node->E[i][j] / row_sum;
      }
    }
    return;
  }
  */

  /* analytic solution (approximation?) */
  for (u=0; u<bntree->num_nodes; u++) {
    node = bntree->preorder[u];
    if (node->parent == NULL) continue;
    
    /* calculate row sums */
    count_sum = 0;
    row_sums = (double *)malloc(bntree->num_states * sizeof(double));
    for (i=0; i<bntree->num_states; i++) {
      row_sums[i] = 0;
      for (j=0; j<bntree->num_states; j++) {
	row_sums[i] += node->E[i][j];
	count_sum += node->E[i][j];
      }
    }

    /* figure out parameter weights */
    for (i=0; i<bntree->num_states; i++) {
      for (j=0; j<bntree->num_states; j++) {
	if (bntree->weight_idx[i][j] == -1)
	  node->W[i][j] = 1;
	else
	  node->W[i][j] = row_sums[bntree->weight_idx[i][j]] / count_sum;
      }
    }

    for (i=0; i<bntree->num_states; i++) {
      for (j=0; j<bntree->num_states; j++) {
	/* calculate weight sum and prob sum for the 
	   parameter in this cell */
	k = bntree->param_map[i][j];
	weight_sum = 0;
	prob_sum = 0;
	
	for (cli = bntree->inv_param_map[k]; cli != NULL; cli = cli->next) {
	  weight_sum += node->W[cli->row][cli->col];
	  prob_sum += node->E[cli->row][cli->col] / 
	    row_sums[cli->row];
	}

	node->P[i][j] = node->W[i][j] / weight_sum * prob_sum;
      }
    }

    /* normalize CPD */
    for (i=0; i<bntree->num_states; i++) {
      row_sum = 0;
      for (j=0; j<bntree->num_states; j++)
	row_sum += node->P[i][j];
      /* fprintf (stderr, "Normalization coefficient %f\n", row_sum); */
      for (j=0; j<bntree->num_states; j++)
	node->P[i][j] /= row_sum;
    }

    /* check constraints */
    for (i=0; i<bntree->num_states; i++) {
      row_sum = 0;
      for (j=0; j<bntree->num_states; j++) {
	row_sum += node->P[i][j];
	if (node->P[i][j] < 0 || node->P[i][j] > 1 ||
	    fabs (node->P[i][j] - 0.0) < 1e-12)
	  fprintf (stderr, "Node %s P[%d][%d] = %f\n", node->name, i, j, node->P[i][j]);
      }
      if (fabs(1.0 - row_sum) > 1e-6)
	fprintf (stderr, "Sum of entries on row %d is %g\n", i, row_sum);
    }

  }
}

void init_params_R0 (BNTree *bntree) {
  int i, j;
  int param_idx;
 
  param_idx = 0;

  /* base to base substitutions */
  for (i=0; i<4; i++) {
    for (j=i+1; j<4; j++) {
      bntree->param_map[i][j] = param_idx;
      bntree->param_map[j][i] = param_idx;

      bntree->weight_idx[i][j] = j;
      bntree->weight_idx[j][i] = i;

      param_idx++;
    }
  }

  /* deletions */
  for (i=0; i<4; i++) {
    bntree->param_map[i][4] = param_idx;
    bntree->weight_idx[i][4] = -1;
  }
  param_idx++;

  /* insertions */
  for (j=0; j<4; j++) {
    bntree->param_map[4][j] = param_idx;
    bntree->weight_idx[4][j] = j;
  }
  param_idx++;

  /* base unaligns */
  for (i=0; i<4; i++) {
    bntree->param_map[i][5] = param_idx;
    bntree->weight_idx[i][5] = -1;
  }
  param_idx++;

  /* gap unaligns */
  bntree->param_map[4][5] = param_idx;
  bntree->weight_idx[4][5] = -1;
  param_idx++;

  /* base realigns */
  for (j=0; j<4; j++) {
    bntree->param_map[5][j] = param_idx;
    bntree->weight_idx[5][j] = j;
  }
  param_idx++;

  /* gap realigns */
  bntree->param_map[5][4] = param_idx;
  bntree->weight_idx[5][4] = -1;
  param_idx++;

  /* self substitutions */
  for (i=0; i<6; i++) {
    bntree->param_map[i][i] = param_idx;
    bntree->weight_idx[i][i] = -1;
    param_idx++;
  }

  assert (param_idx == bntree->num_params);
}

void init_params_G0 (BNTree *bntree) {
  int i, j;
  int param_idx;

  param_idx = 0;

  /* general model, one param per cell */
  for (i=0; i<6; i++) {
    for (j=0; j<6; j++) {
      bntree->param_map[i][j] = param_idx;
      bntree->weight_idx[i][j] = 1;
      param_idx++;
    }
  }

  assert (param_idx == bntree->num_params);
}

void init_params_R1 (BNTree *bntree) {
  int i, j;
  int dimer_pair_num = 0;
  int param_idx = 0;
  int types[] = {0,   1,  2,  3, 16, 17,
		 4,   5,  6,  7, 16, 17,
		 8,   9, 10, 11, 16, 17,
		 12, 13, 14, 15, 16, 17,
		 18, 18, 18, 18, 19, 20,
		 21, 21, 21, 21, 22, 23};

  for (i=0; i<36; i++) {
    for (j=0; j<36; j++) {
      int type1 = types[i];
      int type2 = types[j];
      int char1 = i / 6;
      int char2 = i % 6;
      int char3 = j / 6;
      int char4 = j % 6;

      if (i == j) continue;

      if (type1 < 16 && type2 < 16) { /* base dimer to base dimer substitution */
	if (i > j) continue;  /* already took care of these */

	param_idx = dimer_pair_num++;
	bntree->param_map[i][j] = param_idx;
	bntree->param_map[j][i] = param_idx;
	
	bntree->weight_idx[i][j] = j;
	bntree->weight_idx[j][i] = i;
      }

      else if (type1 == type2 && (char1 != char3 || char2 != char4)) {
	/* mutation in a gap pattern base with no change in pattern */
	if (type1 == 16) {
	  if (is_transition (char1, char3))
	    param_idx = 192;
	  else
	    param_idx = 193;
	}
	else if (type1 == 17) {
	  if (is_transition (char1, char3))
	    param_idx = 194;
	  else
	    param_idx = 195;

	}
	else if (type1 == 18) {
	  if (is_transition (char2, char4))
	    param_idx = 196;
	  else
	    param_idx = 197;
	}
	else if (type1 == 21) {
	  if (is_transition (char2, char4))
	    param_idx = 198;
	  else
	    param_idx = 199;
	}
	bntree->param_map[i][j] = param_idx;	
	bntree->weight_idx[i][j] = -1;
      }

      else { /* gap pattern to different gap pattern */
	param_idx = 120 + 8 * int_max(type1 - 15, 0) + int_max(type2 - 15, 0); 
	if (type1 < type2)
	  param_idx--;

	bntree->param_map[i][j] = param_idx;	

	if (type2 < 16)
	  bntree->weight_idx[i][j] = j;
	else
	  bntree->weight_idx[i][j] = -1;
      }
    }
  }

  /* self substitutions */
  for (i=0; i<36; i++) {
    bntree->param_map[i][i] = 200 + i;
    bntree->weight_idx[i][i] = -1;
  }

  assert (dimer_pair_num == 120);
}


void init_params_G1 (BNTree *bntree) {
  int i, j;
  int param_idx;

  param_idx = 0;

  /* general model, one param per cell */
  for (i=0; i<36; i++) {
    for (j=0; j<36; j++) {
      if (i != j) {
	bntree->param_map[i][j] = param_idx;
	bntree->weight_idx[i][j] = -1;
	param_idx++;
      }
    }
  }

  for (i=0; i<36; i++) {
    bntree->param_map[i][i] = param_idx;
    bntree->weight_idx[i][i] = -1;
    param_idx++;
  }

  assert (param_idx == bntree->num_params);
}

void init_params_G2 (BNTree *bntree) {
  int i, j;
  int param_idx;

  param_idx = 0;

  /* general model, one param per cell */
  for (i=0; i<216; i++) {
    for (j=0; j<216; j++) {
      if (i != j) {
	bntree->param_map[i][j] = param_idx;
	bntree->weight_idx[i][j] = -1;
	param_idx++;
      }
    }
  }

  for (i=0; i<216; i++) {
    bntree->param_map[i][i] = param_idx;
    bntree->weight_idx[i][i] = -1;
    param_idx++;
  }

  assert (param_idx == bntree->num_params);
}


void parse_ss (BNTree *bntree, char *ss_filename, char **tuple_arr, 
	       double *count_arr, double *background_probs) {
  char nseqs_str[STRLEN], 
    length_str[STRLEN], 
    tuple_size_str[STRLEN], 
    ntuples_str[STRLEN],
    names_str[STRLEN], 
    alphabet_str[STRLEN], 
    ncats_str[STRLEN];
  int i, j, k;
  int num_tuples;
  int tuple_size;
  int num_seqs;
  double background_total;
  FILE *fp;

  fp = fopen(ss_filename, "r");
  fscanf (fp, "NSEQS = %s\n", nseqs_str);
  fscanf (fp, "LENGTH = %s\n", length_str);
  fscanf (fp, "TUPLE_SIZE = %s\n", tuple_size_str);
  fscanf (fp, "NTUPLES = %s\n", ntuples_str);
  fscanf (fp, "NAMES = %s\n", names_str);
  fscanf (fp, "ALPHABET = %s\n", alphabet_str);
  fscanf (fp, "NCATS = %s\n", ncats_str);

  num_tuples = atoi (ntuples_str);
  tuple_size = atoi (tuple_size_str);
  num_seqs = atoi (nseqs_str);

  if (tuple_size != bntree->order + 1) {
    printf ("Error, model order is %d but tuple size is %d\n", bntree->order, tuple_size);
    exit (-1);
  }
 
  for (i=0; i<bntree->num_states; i++)
    background_probs[i] = 0;

  for (i=0; i<num_tuples; i++) {
    if (fscanf (fp, "%*d\t%[ACGT_. ]\t%lf\n", tuple_arr[i], &(count_arr[i])) != 2) {
      fprintf (stderr, "Could not read tuple and count for tuple number %d\n", i);
      exit (-1);
    }

    for (j=0; j<num_seqs; j++) {
      char state[STRLEN];
      int state_num;

      for (k=0; k<tuple_size; k++) {
	state[k] = tuple_arr[i][k * (num_seqs + 1) + j];
      }
      state[tuple_size] = '\0';
      state_num = 0;
      for (k=0; k<tuple_size; k++)
	state_num += intpow (ALPHABET_SIZE, (tuple_size - k - 1)) * bntree->symbol_map[(int)state[k]];
      background_probs[state_num] += count_arr[i];
    }
  }

  fclose (fp);

  /* normalize the background prob distribution */
  background_total = 0;
  for (i=0; i<bntree->num_states; i++)
    background_total += background_probs[i];

  for (i=0; i<bntree->num_states; i++)
    background_probs[i] /= background_total;
}

double calc_log_likelihood (BNTree *bntree, char **tuple_arr, 
			    double *count_arr, int num_tuples,
			    int conditional) {
  int i;
  double ll = 0;

  for (i=0; i<num_tuples; i++) {
    ll += count_arr[i] * log (calc_tuple_str_prob (bntree, tuple_arr[i], conditional));
  }

  return ll;
}

int is_transition (int base1, int base2) {
  return ( (base1 == 0 && base2 == 2) ||
	   (base1 == 2 && base2 == 0) ||
	   (base1 == 1 && base2 == 3) ||
	   (base1 == 3 && base2 == 1) );
}

int int_max (int x, int y) {
  if (x > y) return x;
  else return y;
}

/* returns a random number between x and y */
double bounded_random (double x, double y) {
  double val;
  
  if (x >= y) {
    fprintf (stderr, "bounded_random called with lower bound greater than or equal to upper bound\n");
    exit(-1);
  }

  val = x + ( (y - x) * rand() / RAND_MAX );

  return val;
}
