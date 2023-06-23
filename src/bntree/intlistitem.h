#ifndef INTLISTITEM_H
#define INTLISTITEM_H

typedef struct int_list_item IntListItem;

struct int_list_item {
  int val;
  
  IntListItem *next;
  IntListItem *tail;
};

IntListItem* NewIntList ();
void int_list_append (IntListItem *head, int val);

#endif
