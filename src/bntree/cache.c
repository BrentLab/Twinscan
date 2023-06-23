#ifndef CACHE_C
#define CACHE_C

#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "cache.h"

/* taken from Ian Korf's zTools.c in the Zoe package */
static double HASH_MULTIPLIER[7] = {
	3.1415926536, /* PI */
	2.7182818285, /* e */
	1.6180339887, /* golden mean */
	1.7320508076, /* square root of 3 */
	2.2360679775, /* square root of 5 */
	2.6457513111, /* square root of 7 */
	3.3166247904, /* square root of 11 */
};

/* taken from Ian Korf's zTools.c in the Zoe package */
int hash_function (Cache *cache, char *key) {
  int i;
  double sum;

  sum = 0;
  for (i=0; i<cache->key_length; i++) {
    sum += key[i] * HASH_MULTIPLIER[i % 7];
  }
  
  return (int) (cache->num_slots * (sum - floor(sum)));
}

int key_match (Cache *cache, char *key1, char *key2) {
  if (memcmp (key1, key2, cache->key_length))
    return 0;
  else
    return 1;
}

Cache* NewCache (int num_slots, int max_size, int key_length, int array_length) {
  Cache *cache = (Cache *) malloc (sizeof (Cache));
  cache->num_slots = num_slots;
  cache->max_size = max_size;
  cache->key_length = key_length;
  cache->array_length = array_length;
  cache->current_size = 0;
  cache->slots = (CacheEntry **) calloc (num_slots, sizeof (CacheEntry *));
  cache->list_head = NULL;
  cache->list_tail = NULL;

  return cache;
}

CacheEntry* NewCacheEntry(Cache *cache, char *key, float *vals) {
  CacheEntry *entry = (CacheEntry *) malloc (sizeof (CacheEntry));
  entry->key = (char *) malloc (cache->key_length * sizeof (char));
  memcpy (entry->key, key, cache->key_length);
  entry->vals = vals;
  entry->next_in_slot = NULL;
  entry->prev_in_slot = NULL;
  entry->next_in_list = NULL;
  entry->prev_in_list = NULL;

  return entry;
}

void cache_insert (Cache *cache, char *key, float *vals) {
  CacheEntry *entry = NewCacheEntry (cache, key, vals);
  int slot = hash_function (cache, key);

  if (cache->slots[slot] == NULL) {
    /* this is the first entry in the slot */
    cache->slots[slot] = entry;
  }
  else {
    /* find the last entry in the given slot */
    CacheEntry *last_in_slot = cache->slots[slot];
    while (last_in_slot->next_in_slot != NULL) {
      last_in_slot = last_in_slot->next_in_slot;
    }
    /* append a new cache entry to it */
    last_in_slot->next_in_slot = entry;
    entry->prev_in_slot = last_in_slot;
  }
 
  if (cache->list_head == NULL) {
    /* this is the first entry in the cache */
    cache->list_head = entry;
    cache->list_tail = entry;
  }
  else {
    /* add the new entry to the back of the list */
    cache->list_tail->next_in_list = entry;
    entry->prev_in_list = cache->list_tail;
    cache->list_tail = entry;
  }

  /* if we're out of space, delete the first entry in the list */
  if (cache->current_size >= cache->max_size)
    cache_remove (cache, cache->list_head);
  else
    cache->current_size += sizeof (CacheEntry) + 
      cache->key_length * sizeof (char) +
      cache->array_length * sizeof (float);
}

float* cache_get (Cache *cache, char *key) {
  CacheEntry *entry;
  int slot = hash_function (cache, key);
  
  if (cache->slots[slot] == NULL)
    return NULL; /* not found */

  entry = cache->slots[slot];
  while (entry->next_in_slot != NULL) {
    if (key_match (cache, key, entry->key)) {
      move_to_back_of_list (cache, entry);
      return entry->vals;
    }
    entry = entry->next_in_slot;
  }

  /* we're on the last entry in the slot */
  if (key_match (cache, key, entry->key)) {
    move_to_back_of_list (cache, entry);
    return entry->vals;
  }
  else
    return NULL; /* not found */
}

void cache_remove (Cache *cache, CacheEntry *entry) {
  
  /* remove from slot */
  int slot = hash_function (cache, entry->key);
  if (cache->slots[slot] == entry) {
    cache->slots[slot] = entry->next_in_slot;
  }
  if (entry->prev_in_slot != NULL) {
    entry->prev_in_slot->next_in_slot = entry->next_in_slot;
  }
  if (entry->next_in_slot != NULL) {
    entry->next_in_slot->prev_in_slot = entry->prev_in_slot;
  }

  /* remove from list */
  if (cache->list_head == entry) {
    cache->list_head = entry->next_in_list;
  }
  if (cache->list_tail == entry) {
    cache->list_tail = entry->prev_in_list;
  }
  if (entry->prev_in_list != NULL) {
    entry->prev_in_list->next_in_list = entry->next_in_list;
  }
  if (entry->next_in_list != NULL) {
    entry->next_in_list->prev_in_list = entry->prev_in_list;
  }
  
  /* free memory */
  cache_entry_free (entry);
}

void move_to_back_of_list (Cache *cache, CacheEntry *entry) {
  if (entry->next_in_list != NULL)
    entry->next_in_list->prev_in_list = entry->prev_in_list;
  else
    return;  /* this entry is already at the back of the list */
  
  if (entry->prev_in_list != NULL)
    entry->prev_in_list->next_in_list = entry->next_in_list;

  if (entry == cache->list_head)
    cache->list_head = entry->next_in_list;

  cache->list_tail->next_in_list = entry;
  entry->prev_in_list = cache->list_tail;
  entry->next_in_list = NULL;
  cache->list_tail = entry;
}

void cache_free (Cache *cache) {
  if (cache->list_head != NULL) {
    /* free all entries */
    CacheEntry *entry = cache->list_head;
    
    while (entry->next_in_list != NULL) {
      entry = entry->next_in_list;
      cache_entry_free (entry->prev_in_list);
    }
    cache_entry_free (entry);
  }    

  /* free pointer array */
  free (cache->slots);

  /* free the cache structure */
  free (cache);
}

void cache_entry_free (CacheEntry *entry) {
  free (entry->key);
  free (entry->vals);
  free (entry);
}

#endif
