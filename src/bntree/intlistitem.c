#ifndef INTLISTITEM_C
#define INTLISTITEM_C

#include <stdlib.h>
#include "intlistitem.h"

IntListItem* NewIntList() {
  IntListItem *new_list = (IntListItem *)malloc(sizeof(IntListItem));
  
  new_list->next = NULL;
  new_list->tail = NULL;
  
  return new_list;
}

void int_list_append (IntListItem *head, int val) {

  if (head->tail == NULL) { /* empty list */
    head->val = val;
    head->tail = head;
  }

  else {
    IntListItem *new = (IntListItem *)malloc(sizeof(IntListItem));
    new->val = val;
    new->next = NULL;
    new->tail = NULL;
    
    head->tail->next = new;
    head->tail = new;
  }
}

#endif
