#ifndef NODE_C
#define NODE_C

#include "node.h"

Node* NewNode(int id) {
  Node* node = (Node *)malloc(sizeof(Node));
  node->id = id;
  node->children = (Node **)malloc(2 * sizeof(Node *));
  node->num_children = 0;
  strcpy (node->name, "");
  strcpy (node->seq, "");

  return node;
}

#endif
