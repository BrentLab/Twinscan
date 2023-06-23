#ifndef CACHE_H
#define CACHE_H

/* Implements a cache that uses byte arrays as keys
   and float arrays as their associated values.  Initialized
   with a maximum size in bytes.  Uses LRU algorithm
   to try to minimize cache misses given the limited
   amount of space. */

struct cache_entry {
  char *key;
  float *vals;
  struct cache_entry *next_in_slot;
  struct cache_entry *prev_in_slot;
  struct cache_entry *next_in_list;
  struct cache_entry *prev_in_list;
};
typedef struct cache_entry CacheEntry;

struct cache {
  int num_slots;
  int max_size;
  int current_size; /* in bytes */
  int key_length;
  int array_length;
  CacheEntry **slots;
  CacheEntry *list_head;
  CacheEntry *list_tail;
};
typedef struct cache Cache;

int hash_function (Cache *cache, char *key);
int key_match (Cache *cache, char *key1, char *key2);
CacheEntry* NewCacheEntry(Cache *cache, char *key, float *vals);
Cache* NewCache(int num_slots, int max_size, int key_length, int array_length);
void cache_insert (Cache *cache, char *key, float *vals);
float* cache_get (Cache *cache, char *key);
void cache_remove (Cache *cache, CacheEntry *entry);
void move_to_back_of_list (Cache *cache, CacheEntry *entry);
void cache_free(Cache *cache);
void cache_entry_free(CacheEntry *entry);

#endif
