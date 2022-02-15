#include "hashtable.h"

static int pointercmp(const void *pointer1, const void *pointer2) {
	return (pointer1 != pointer2);
}

static uint32_t pointerHashFunction(const void *pointer) {
	return ((uint32_t) pointer) ;
}

static int isProbablePrime(size_t oddNumber) {
	size_t i;

	for (i=3; i<51; i+=2)
		if (oddNumber == i)
			return 1;
		else if (oddNumber%i == 0)
			return 0;

	return 1; /* maybe */
}

RedBlackTree *RedBlackTreeInit()
{
	RedBlackTree *rbtree = (RedBlackTree *)malloc(sizeof(RedBlackTree));
	rbtree->root = NULL;
	rbtree->keycmp = NULL;
	rbtree->valuecmp = NULL;
	RedBlackTreeNode_t *tmp = (RedBlackTreeNode_t *)malloc(sizeof(RedBlackTreeNode_t));
	tmp->parent = tmp->left = tmp->right = tmp;
	tmp->val = NULL;
	tmp->key = NULL;
	tmp->red = 0;
	rbtree->nil = tmp;
	tmp = (RedBlackTreeNode_t *)malloc(sizeof(RedBlackTreeNode_t));
	tmp->val = NULL;
	tmp->key = NULL;
	tmp->red = 0;
	tmp->parent = tmp->left = tmp->right = rbtree->nil;
	rbtree->root = tmp;
	return rbtree;
}

static RedBlackTreeNode_t *RedBlackTreeNodeInit(RedBlackTree *rbtree)
{
	RedBlackTreeNode_t *node = (RedBlackTreeNode_t *)malloc(sizeof(RedBlackTreeNode_t));
	node->left = node->right = node->parent = rbtree->nil;
	return node;
}

uint8_t RedBlackTreeEmpty(RedBlackTree *rbtree)
{
	return rbtree->root == NULL;
}

static void RedBlackTreeLrotate(RedBlackTree *rbtree, RedBlackTreeNode_t *x)
{
	RedBlackTreeNode_t *nil = rbtree->nil;
	RedBlackTreeNode_t *y;
    y = x->right;
    x->right = y->left;
    if(y->left != nil) y->left->parent=x;
    y->parent = x->parent;   
    if(x == x->parent->left) {
        x->parent->left = y;
    } else {
        x->parent->right = y;
    }
    y->left = x;
    x->parent = y;
}

static void RedBlackTreeRrotate(RedBlackTree *rbtree, RedBlackTreeNode_t *y)
{	
	RedBlackTreeNode_t *nil = rbtree->nil;
	RedBlackTreeNode_t *x;
    x = y->left;
    y->left = x->right;
    if(nil != x->right) x->right->parent = y;
    x->parent = y->parent;
    if(y == y->parent->left) {
        y->parent->left = x;
    } else {
        y->parent->right = x;
    }
    x->right = y;
    y->parent = x;
}

static void __RedBlackTreeInsert(RedBlackTree *rbtree, RedBlackTreeNode_t *node)
{
	RedBlackTreeNode_t *nil = rbtree->nil;
	if (!rbtree->keycmp) {
		fatal("KeyCompare function should be set for RedBlackTree.");
	}
	RedBlackTreeNode_t *x, *y;
	node->left = node->right = nil;
	y = rbtree->root;
	x = rbtree->root->left;
	while (x != nil) {
		y = x;
		if ( rbtree->keycmp(x->key, node->key) > 0) {
			x = x->left;
		} else {
			x = x->right;
		}
	}
	node->parent = y;
	if ((y == rbtree->root) || rbtree->keycmp(y->key, node->key) > 0) {
		y->left = node;
	} else {
		y->right = node;
	}
}

static void _RedBlackTreeInsert(RedBlackTree *rbtree, void *key, void *val, RedBlackTreeNode_t *node)
{
	RedBlackTreeNode_t *x, *y;
	x = node;
	x->key = key;
	x->val = val;
	__RedBlackTreeInsert(rbtree, x);
	x->red = 1;
	
	while (x->parent->red) {
		if (x->parent == x->parent->parent->left) {
			y = x->parent->parent->right;
			if (y->red) {
				x->parent->red = 0;
				y->red = 0;
				x->parent->parent->red = 1;
				x = x->parent->parent;
			} else {
				if ( x == x->parent->right) {
					x = x->parent;
					RedBlackTreeLrotate(rbtree, x);
				}
				x->parent->red = 0;
				x->parent->parent->red = 1;
				RedBlackTreeRrotate(rbtree, x->parent->parent);
			}
		} else {
			y = x->parent->parent->left; 
			if (y->red) {
				x->parent->red = 0;
				y->red = 0;
				x->parent->parent->red = 1;
				x = x->parent->parent;
			} else {
				if (x == x->parent->left) {
					x = x->parent;
					RedBlackTreeRrotate(rbtree, x);
				}
				x->parent->red = 0;
				x->parent->parent->red = 1;
				RedBlackTreeLrotate(rbtree, x->parent->parent);
			}
		}
	}
	rbtree->root->left->red = 0;
}

void RedBlackTreeInsert(RedBlackTree *rbtree, void *key, void *val)
{
	RedBlackTreeNode_t *node = RedBlackTreeNodeInit(rbtree);
	_RedBlackTreeInsert(rbtree, key, val, node);
}

RedBlackTreeNode_t *RedBlackTreeFirst(RedBlackTree *rbtree)
{
	RedBlackTreeNode_t *nil = rbtree->nil;
	RedBlackTreeNode_t *x;
	if (!(x = rbtree->root)) {
		return NULL;
	}
	while (x->left != nil) x = x->left;
	return x;
}

static RedBlackTreeNode_t *_RedBlackTreeNext(RedBlackTree *rbtree, RedBlackTreeNode_t *x)
{
	RedBlackTreeNode_t *nil = rbtree->nil;
	RedBlackTreeNode_t *y;
	RedBlackTreeNode_t *root = rbtree->root;
	if ((y = x->right ) != nil) {
		while (y->left != nil) {
			y = y->left;
		}
		return y;
	} else {
		y = x->parent;
		while (x == y->right) {
			x = y;
			y = y->parent;
		}
		if (y == root) return nil;
		return y;
	}
}

RedBlackTreeNode_t *RedBlackTreeNext(RedBlackTree *rbtree, RedBlackTreeNode_t *it)
{
	RedBlackTreeNode_t *x = _RedBlackTreeNext(rbtree, it);
	return x;
}

static void _RedBlackTreeTraversal(RedBlackTree *rbtree, RedBlackTreeNode_t *node, void (*mapfunc)(void *key, void *value))
{
	RedBlackTreeNode_t *nil = rbtree->nil;
	if (node->left && node->left != nil) _RedBlackTreeTraversal(rbtree, node->left, mapfunc);
	if ((node != rbtree->root) && (node != rbtree->nil)) mapfunc(node->key, node->val);
	if (node->right && node->right != nil) _RedBlackTreeTraversal(rbtree, node->right, mapfunc);
}

void RedBlackTreeTraversal(RedBlackTree *rbtree, void (*mapfunc)(void *key, void *value))
{
	_RedBlackTreeTraversal(rbtree, rbtree->root->left, mapfunc);
}

static void RedBlackTreeFixup(RedBlackTree *rbtree, RedBlackTreeNode_t *x)
{
	RedBlackTreeNode_t *root = rbtree->root->left;
	RedBlackTreeNode_t *w;
	while ((!x->red) && (root != x)) {
		if (x == x->parent->left) {
			w = x->parent->right;
			if (w->red) {
				w->red = 0;
				x->parent->red = 1;
				RedBlackTreeLrotate(rbtree, x->parent);
				w = x->parent->right;
			}
		
			if ((!w->right->red) && (!w->left->red)) {
				w->red = 1;
				x = x->parent;
			} else {
				if (!w->right->red) {
					w->left->red = 0;
					w->red = 1;
					RedBlackTreeRrotate(rbtree, w);
					w = x->parent->right;
				}
				w->red = x->parent->red;
				x->parent->red = 0;
				w->right->red = 0;
				RedBlackTreeLrotate(rbtree, x->parent);
				x = root;
			}
		} else {
			w = x->parent->left;
			if (w->red) {
				w->red = 0;
				x->parent->red = 1;
				RedBlackTreeRrotate(rbtree, x->parent);
				w = x->parent->left;
			}
			if ((!w->right->red) && (!w->left->red)) {
				w->red = 1;
				x = x->parent;
			} else {
				if (!w->left->red) {
					w->right->red = 0;
					w->red = 1;
					RedBlackTreeLrotate(rbtree, w);
					w = x->parent->left;
				}
				w->red = x->parent->red;
				x->parent->red = 0;
				w->left->red = 0;
				RedBlackTreeRrotate(rbtree, x->parent);
				x = root;
			}	
		}
	}
	x->red = 0;
}

void RedBlackTreeErase(RedBlackTree *rbtree, RedBlackTreeNode_t *node)
{
	RedBlackTreeNode_t *x, *y;
	RedBlackTreeNode_t *nil = rbtree->nil;
	RedBlackTreeNode_t *root = rbtree->root;
	y = ((node->left == nil) || (node->right == nil)) ? node : _RedBlackTreeNext(rbtree, node);
	x = y->left == nil ? y->right : y->left;
	if (root == (x->parent = y->parent)) {
		root->left = x;
	} else {
		if (y == y->parent->left) {
			y->parent->left = x;
		} else {
			y->parent->right = x;
		}
	}
	if (y != node) {
		if (!(y->red)) RedBlackTreeFixup(rbtree, x);
		y->left = node->left;
		y->right = node->right;
		y->parent = node->parent;
		y->red = node->red;
		node->left->parent = node->right->parent = y;
		if (node == node->parent->left) {
			node->parent->left = y;
		} else {
			node->parent->right = y;
		}
	} else {
		if (!(y->red)) RedBlackTreeFixup(rbtree, x);
	}
}

void RedBlackTreeSetKeyCompare(RedBlackTree *rbtree, int8_t (*keycmp)(const void *, const void *))
{
	rbtree->keycmp = keycmp;
}

void RedBlackTreeSetValueCompare(RedBlackTree *rbtree, int8_t (*valcmp)(const void *, const void *))
{
	rbtree->valuecmp = valcmp;
}

static void _RedBlackTreeDestroy(RedBlackTree *rbtree, RedBlackTreeNode_t *node)
{
	if (node == rbtree->nil) return;
	_RedBlackTreeDestroy(rbtree, node->left);
	_RedBlackTreeDestroy(rbtree, node->right);
	free(node);
}

void RedBlackTreeDestroy(RedBlackTree *rbtree)
{
	RedBlackTreeNode_t *node = rbtree->root;
	_RedBlackTreeDestroy(rbtree, node);
	free(rbtree->nil);
	free(rbtree);
}

static size_t calculateIdealNumOfBuckets(HashTable *hashTable) {
	size_t idealNumOfBuckets = hashTable->numOfElements / hashTable->idealRatio;
	if (idealNumOfBuckets < 5)
		idealNumOfBuckets = 5;
	else
		idealNumOfBuckets |= 0x01; /* make it an odd number */
	while (!isProbablePrime(idealNumOfBuckets))
		idealNumOfBuckets += 2;

	return idealNumOfBuckets;
}

HashTable *HashTableCreate(size_t numOfBuckets) {
	HashTable *hashTable;
	int i;
	assert(numOfBuckets > 0);
	hashTable = (HashTable *) malloc(sizeof(HashTable));
	if (hashTable == NULL)
		return NULL;

	hashTable->appendix1=NULL;
	hashTable->appendix2=NULL;
	hashTable->appendix3=NULL;

	hashTable->counter1=0;
	hashTable->counter2=0;
	hashTable->counter3=0;

	hashTable->bucketArray = (KeyValuePair **)
						malloc(numOfBuckets * sizeof(KeyValuePair *));
	if (hashTable->bucketArray == NULL) {
		free(hashTable);
		return NULL;
	}
	
	hashTable->numOfBuckets = numOfBuckets;
	hashTable->numOfElements = 0;

	for (i=0; i<numOfBuckets; i++)
		hashTable->bucketArray[i] = NULL;

	hashTable->idealRatio = 3.0;
	hashTable->lowerRehashThreshold = 0.0;
	hashTable->upperRehashThreshold = 15.0;

	hashTable->keycmp = pointercmp;
	hashTable->valuecmp = pointercmp;
	hashTable->hashFunction = pointerHashFunction;
	hashTable->keyDeallocator = NULL;
	hashTable->valueDeallocator = NULL;

	return hashTable;
}

void HashTableDestroy(HashTable *hashTable) {
	int i;

	for (i=0; i<hashTable->numOfBuckets; i++) {
		KeyValuePair *pair = hashTable->bucketArray[i];
		while (pair != NULL) {
			KeyValuePair *nextPair = pair->next;

			if (hashTable->keyDeallocator != NULL)
				hashTable->keyDeallocator((void *) pair->key);
			if (hashTable->valueDeallocator != NULL){
//			fprintf(stderr,"FREE %p\n", pair->value);
			hashTable->valueDeallocator(pair->value);
		}
			free(pair);
			pair = nextPair;
		}
	}

	free(hashTable->bucketArray);
	free(hashTable);
}


int HashTableContainsKey(const HashTable *hashTable, const void *key) {
	return (HashTableGet(hashTable, key) != NULL);
}


int HashTableContainsValue(const HashTable *hashTable, const void *value) {
	int i;

	for (i=0; i<hashTable->numOfBuckets; i++) {
		KeyValuePair *pair = hashTable->bucketArray[i];
		while (pair != NULL) {
			if (hashTable->valuecmp(value, pair->value) == 0)
				return 1;
			pair = pair->next;
		}
	}

	return 0;
}


int HashTablePut(HashTable *hashTable, const void *key, void *value) {
	return HashTablePutReplace(hashTable, key, value, 1);
}

int HashTablePutReplaceEx(HashTable *hashTable, const void *key, void *value, int replace_key, int dealloc_key, int dealloc_value) {
	size_t hashValue;
	KeyValuePair *pair;

	assert(key != NULL);
	assert(value != NULL);

	hashValue = hashTable->hashFunction(key) % hashTable->numOfBuckets;


	pair = hashTable->bucketArray[hashValue];

	while (pair != NULL && hashTable->keycmp(key, pair->key) != 0)
		pair = pair->next;

	if (pair) {
		if (pair->key != key) {
			if(replace_key) {
				if(hashTable->keyDeallocator && dealloc_key)
					hashTable->keyDeallocator((void *) pair->key);
				pair->key = key;
			}
		}
		if (pair->value != value) {
			if (hashTable->valueDeallocator != NULL && dealloc_value)
				hashTable->valueDeallocator(pair->value);
			pair->value = value;
		}
	}
	else {
		KeyValuePair *newPair = (KeyValuePair *) malloc(sizeof(KeyValuePair));
		if (newPair == NULL) {
		return -1;
		}
		else {
			newPair->key = key;
			newPair->value = value;
			newPair->next = hashTable->bucketArray[hashValue];
			hashTable->bucketArray[hashValue] = newPair;
			hashTable->numOfElements++;

			if (hashTable->upperRehashThreshold > hashTable->idealRatio) {
				float elementToBucketRatio = (float) hashTable->numOfElements /
											 (float) hashTable->numOfBuckets;
				if (elementToBucketRatio > hashTable->upperRehashThreshold)
					HashTableRehash(hashTable, 0);
			}
		}
	}

	return 0;
}
int HashTablePutReplace(HashTable *hashTable, const void *key, void *value, int replace_key) {
	return HashTablePutReplaceEx(hashTable, key, value,  replace_key, 1, 1);
}


void *HashTableGet(const HashTable *hashTable, const void *key) {
	size_t hashValue = hashTable->hashFunction(key) % hashTable->numOfBuckets;
	KeyValuePair *pair = hashTable->bucketArray[hashValue];
	while (pair != NULL && hashTable->keycmp(key, pair->key) != 0)
		pair = pair->next;
	return (pair == NULL)? NULL : pair->value;
}

void HashTableRemove(HashTable *hashTable, const void *key) {
	size_t hashValue = hashTable->hashFunction(key) % hashTable->numOfBuckets;


	KeyValuePair *pair = hashTable->bucketArray[hashValue];
	KeyValuePair *previousPair = NULL;

	while (pair != NULL && hashTable->keycmp(key, pair->key) != 0) {
		previousPair = pair;
		pair = pair->next;
	}

	if (pair != NULL) {
		if (hashTable->keyDeallocator != NULL)
			hashTable->keyDeallocator((void *) pair->key);
		if (hashTable->valueDeallocator != NULL)
			hashTable->valueDeallocator(pair->value);
		if (previousPair != NULL)
			previousPair->next = pair->next;
		else
			hashTable->bucketArray[hashValue] = pair->next;
		free(pair);
		hashTable->numOfElements--;

		if (hashTable->lowerRehashThreshold > 0.0) {
			float elementToBucketRatio = (float) hashTable->numOfElements /
										 (float) hashTable->numOfBuckets;
			if (elementToBucketRatio < hashTable->lowerRehashThreshold)
				HashTableRehash(hashTable, 0);
		}
	}


}

void HashTableRemoveAll(HashTable *hashTable) {
	int i;

	for (i=0; i<hashTable->numOfBuckets; i++) {
		KeyValuePair *pair = hashTable->bucketArray[i];
		while (pair != NULL) {
			KeyValuePair *nextPair = pair->next;
			if (hashTable->keyDeallocator != NULL)
				hashTable->keyDeallocator((void *) pair->key);
			if (hashTable->valueDeallocator != NULL)
				hashTable->valueDeallocator(pair->value);
			free(pair);
			pair = nextPair;
		}
		hashTable->bucketArray[i] = NULL;
	}

	hashTable->numOfElements = 0;
	HashTableRehash(hashTable, 5);
}


int HashTableIsEmpty(const HashTable *hashTable) {
	return (hashTable->numOfElements == 0);
}


size_t HashTableSize(const HashTable *hashTable) {
	return hashTable->numOfElements;
}


size_t HashTableGetNumBuckets(const HashTable *hashTable) {
	return hashTable->numOfBuckets;
}

void HashTableSetKeyComparisonFunction(HashTable *hashTable,
		int (*keycmp)(const void *key1, const void *key2)) {
	assert(keycmp != NULL);
	hashTable->keycmp = keycmp;
}


void HashTableSetValueComparisonFunction(HashTable *hashTable,
		int (*valuecmp)(const void *value1, const void *value2)) {
	assert(valuecmp != NULL);
	hashTable->valuecmp = valuecmp;
}


void HashTableSetHashFunction(HashTable *hashTable,
		uint32_t (*hashFunction)(const void *key))
{
	assert(hashFunction != NULL);
	hashTable->hashFunction = hashFunction;
}


void HashTableRehash(HashTable *hashTable, size_t numOfBuckets) {
	KeyValuePair **newBucketArray;
	int i;

	assert(numOfBuckets >= 0);
	if (numOfBuckets == 0)
		numOfBuckets = calculateIdealNumOfBuckets(hashTable);

	if (numOfBuckets == hashTable->numOfBuckets)
		return; /* already the right size! */

	newBucketArray = (KeyValuePair **)
								malloc(numOfBuckets * sizeof(KeyValuePair *));
	if (newBucketArray == NULL) {
		/* Couldn't allocate memory for the new array.  This isn't a fatal
		 * error; we just can't perform the rehash. */
		return;
	}

	for (i=0; i<numOfBuckets; i++)
		newBucketArray[i] = NULL;

	for (i=0; i<hashTable->numOfBuckets; i++) {
		KeyValuePair *pair = hashTable->bucketArray[i];
		while (pair != NULL) {
			KeyValuePair *nextPair = pair->next;
			size_t hashValue = hashTable->hashFunction(pair->key) % numOfBuckets;
			pair->next = newBucketArray[hashValue];
			newBucketArray[hashValue] = pair;
			pair = nextPair;
		}
	}

	free(hashTable->bucketArray);
	hashTable->bucketArray = newBucketArray;
	hashTable->numOfBuckets = numOfBuckets;
}


void HashTableSetIdealRatio(HashTable *hashTable, float idealRatio,
		float lowerRehashThreshold, float upperRehashThreshold) {
	assert(idealRatio > 0.0);
	assert(lowerRehashThreshold < idealRatio);
	assert(upperRehashThreshold == 0.0 || upperRehashThreshold > idealRatio);

	hashTable->idealRatio = idealRatio;
	hashTable->lowerRehashThreshold = lowerRehashThreshold;
	hashTable->upperRehashThreshold = upperRehashThreshold;
}

void HashTableSetDeallocationFunctions(HashTable *hashTable,
		void (*keyDeallocator)(void *key),
		void (*valueDeallocator)(void *value)) {
	hashTable->keyDeallocator = keyDeallocator;
	hashTable->valueDeallocator = valueDeallocator;
}

/**
 * @brief __ac_X31_hash_string
 * 
 * @param key 
 * @return uint32_t 
 */
uint32_t HashTableStringHashFunction(const void *key) {
    char *s = (char *)key;
	uint32_t h = (uint32_t)*s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + (uint32_t)*s;
	return h;

}

void free_values_destroy(HashTable * tab)
{
	KeyValuePair * cursor;
	int bucket;
	
	for(bucket=0; bucket< tab -> numOfBuckets; bucket++)
	{
		cursor = tab -> bucketArray[bucket];
		while (1)
		{
			if(!cursor) break;
			char * read_txt = (char *) cursor ->value;
			free(read_txt);
			cursor = cursor->next;
		}
	}
	HashTableDestroy(tab);
}

void HashTableIteration(HashTable * tab, void process_item(void * key, void * hashed_obj, HashTable * tab) )
{
	int i;
	for (i=0; i< tab ->numOfBuckets; i++) {
		KeyValuePair *pair = tab ->bucketArray[i];
		while (pair != NULL) {
			process_item(( void * )pair -> key, pair -> value, tab);
			KeyValuePair *nextPair = pair->next;
			pair = nextPair;
		}
	}
}