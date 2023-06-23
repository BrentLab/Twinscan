#include <stdio.h>
#include <assert.h>
#include "bntree.h"
#include "training.h"

int main (int argc, char** argv) {
  char *treestr;
  char *modelstr;
  char *ss_filename;
  char *model_name;
  char *output_filename;

  if (argc != 6) {
    fprintf (stderr, "Usage:  %s <tree> <model type> <SS file> <model name> <output file>\n", argv[0]);
    exit(-1);
  }

  treestr = argv[1];
  modelstr = argv[2];
  ss_filename = argv[3];
  model_name = argv[4];
  output_filename = argv[5];
  
  train_bntree (treestr, modelstr, ss_filename, model_name, output_filename);

  return 0;
}
