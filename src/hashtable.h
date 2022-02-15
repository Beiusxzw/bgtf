#ifndef BGTF_HASHTABLE_H
#define BGTF_HASHTABLE_H

#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include "utils.h"
#include "list.h"




/* These structs should not be accessed directly from user code.
 * All access should be via the public functions declared below. */

typedef struct RedBlackTreeNode RedBlackTreeNode_t;
struct RedBlackTreeNode
{
    uint8_t red;
    RedBlackTreeNode_t *left;
    RedBlackTreeNode_t *right;
    RedBlackTreeNode_t *parent;
    void *key;
    void *val;
};

typedef struct RedBlackTree {
    RedBlackTreeNode_t *root;
    RedBlackTreeNode_t *nil;
    int8_t (*keycmp)(const void *key1, const void *key2);
    int8_t (*valuecmp)(const void *value1, const void *value2);
    void (*keyDeallocator)(void *key);
    void (*valueDeallocator)(void *value);
} RedBlackTree;

RedBlackTree *RedBlackTreeInit();
uint8_t RedBlackTreeEmpty(RedBlackTree *rbtree);
void RedBlackTreeInsert(RedBlackTree *rbtree, void *key, void *val);
void RedBlackTreeErase(RedBlackTree *rbtree, RedBlackTreeNode_t *node);
RedBlackTreeNode_t *RedBlackTreeFirst(RedBlackTree *rbtree);
RedBlackTreeNode_t *RedBlackTreeNext(RedBlackTree *rbtree, RedBlackTreeNode_t *it);
void RedBlackTreeSetKeyCompare(RedBlackTree *rbtree, int8_t (*keycmp)(const void *, const void *));
void RedBlackTreeSetValueCompare(RedBlackTree *rbtree, int8_t (*valuecmp)(const void *, const void *));
void RedBlackTreeDestroy(RedBlackTree *rbtree);
void RedBlackTreeTraversal(RedBlackTree *rbtree, void (*mapfunc)(void *key, void *value));

/* These structs should not be accessed directly from user code.
 * All access should be via the public functions declared below. */

typedef struct KeyValuePair
{
    const void *key;
    void *value;
    struct KeyValuePair *next;
} KeyValuePair;

typedef struct
{
    size_t numOfBuckets;
    size_t numOfElements;
    KeyValuePair **bucketArray;
    float idealRatio, lowerRehashThreshold, upperRehashThreshold;
    int (*keycmp)(const void *key1, const void *key2);
    int (*valuecmp)(const void *value1, const void *value2);
    uint32_t (*hashFunction)(const void *key);
    void (*keyDeallocator)(void *key);
    void (*valueDeallocator)(void *value);

    void *appendix1;
    void *appendix2;
    void *appendix3;
    uint64_t counter1;
    uint64_t counter2;
    uint64_t counter3;
} HashTable;

void HashTableIteration(HashTable *tab, void process_item(void *key, void *hashed_obj, HashTable *tab));
HashTable *HashTableCreate(size_t numOfBuckets);
void HashTableDestroy(HashTable *hashTable);
int HashTableContainsKey(const HashTable *hashTable, const void *key);
int HashTableContainsValue(const HashTable *hashTable, const void *value);
int HashTablePut(HashTable *hashTable, const void *key, void *value);
int HashTablePutReplace(HashTable *hashTable, const void *key, void *value, int replace_key);
void *HashTableGet(const HashTable *hashTable, const void *key);
void HashTableRemove(HashTable *hashTable, const void *key);
void HashTableRemoveAll(HashTable *hashTable);
int HashTableIsEmpty(const HashTable *hashTable);
size_t HashTableSize(const HashTable *hashTable);
size_t HashTableGetNumBuckets(const HashTable *hashTable);
void HashTableSetKeyComparisonFunction(HashTable *hashTable,
                                       int (*keycmp)(const void *key1, const void *key2));
void HashTableSetValueComparisonFunction(HashTable *hashTable,
                                         int (*valuecmp)(const void *value1, const void *value2));
void HashTableSetHashFunction(HashTable *hashTable,
                              uint32_t (*hashFunction)(const void *key));
void HashTableRehash(HashTable *hashTable, size_t numOfBuckets);
void HashTableSetIdealRatio(HashTable *hashTable, float idealRatio,
                            float lowerRehashThreshold,
                            float upperRehashThreshold);
void HashTableSetDeallocationFunctions(HashTable *hashTable,
                                       void (*keyDeallocator)(void *key),
                                       void (*valueDeallocator)(void *value));
uint32_t HashTableStringHashFunction(const void *key);

void free_values_destroy(HashTable *tab);
int HashTablePutReplaceEx(HashTable *hashTable, const void *key, void *value, int replace_key, int dealloc_key, int dealloc_value);

/**
 * @brief Macro for HashTable key-value iteration
 * 
 */
#define HashTableForEach(hashtable, kelem, velem, code) \
{\
    KeyValuePair *pair; \
    void *key, *value; \
    for (size_t i = 0; i < hashtable->numOfBuckets; ++i) { \
        KeyValuePair *pair = hashtable->bucketArray[i]; \
        while (pair != NULL) { \
            kelem = pair->key; \
            velem = pair->value; \
            code; \
            pair = pair->next; \
        } \
    } \
} while (0);\

#endif