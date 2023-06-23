#ifndef CELLLISTITEM_C
#define CELLLISTITEM_C

#include <stdlib.h>
#include "celllistitem.h"

CellListItem* NewCellList() {
  CellListItem *new_list = (CellListItem *)malloc(sizeof(CellListItem));
  
  new_list->next = NULL;
  new_list->tail = NULL;
  
  return new_list;
}

void cell_list_append (CellListItem *head, int row, int col) {

  if (head->tail == NULL) { /* empty list */
    head->row = row;
    head->col = col;
    head->tail = head;
  }

  else {
    CellListItem *new = (CellListItem *)malloc(sizeof(CellListItem));
    new->row = row;
    new->col = col;
    new->next = NULL;
    new->tail = NULL;
    
    head->tail->next = new;
    head->tail = new;
  }
}

#endif
