#ifndef TRAINING_H
#define TRAINING_H

#include "bntree.h"

typedef enum {
  TT0,
  TT1,
  TT2,
  R0,
  R1,
  R2,
  G0,
  G1,
  G2
} probability_model_type;

BNTree* train_bntree (char *treestr, char *modelstr, char *ss_filename,
		      char *model_name, char *output_filename);

/* The init_params_R? functions had the following signature, but the arg background_probs
   wasn't used. I removed them in R0 and R1. If it is necessary, here is the originial
   signature.
		void init_params_R0 (BNTree *bntree, double *background_probs);
*/

void init_params_R0 (BNTree *bntree);
void init_params_G0 (BNTree *bntree);

void init_params_R1 (BNTree *bntree);
void init_params_G1 (BNTree *bntree);

void init_params_G2 (BNTree *bntree);

void parse_ss (BNTree *bntree, char *ss_filename, char **tuple_arr, 
	       double *count_arr, double *background_probs);
double calc_log_likelihood (BNTree *bntree, char **tuple_arr, 
			    double *count_arr, int num_tuples,
			    int conditional);
void Estep (BNTree *bntree, char **tuple_arr, double *count_arr, int num_tuples);
/* The function Mstep had the following signature, but the arg model_type
   wasn't used. I removed it. If it is necessary, here is the originial
   signature.
	void Mstep (BNTree *bntree, probability_model_type model_type);
*/
void Mstep (BNTree *bntree);
int is_transition (int base1, int base2);
int int_max (int x, int y);
double bounded_random (double x, double y);

#endif
