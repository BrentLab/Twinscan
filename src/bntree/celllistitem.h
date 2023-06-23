#ifndef CELLLISTITEM_H
#define CELLLISTITEM_H

typedef struct cell_list_item CellListItem;

struct cell_list_item {
  int row;
  int col;
  
  CellListItem *next;
  CellListItem *tail;
};

CellListItem* NewCellList ();
void cell_list_append (CellListItem *head, int row, int col);

#endif
