// Copyright 2020 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

#ifndef RTREE_H
#define RTREE_H


#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "utils.h"


#define ALLOW_REINSERTS
#define MAXITEMS 32      // 32 max items per node
#define MINFILL  20      // 20% min fill


#define rtmalloc (_malloc?_malloc:malloc)
#define rtfree (_free?_free:free)



typedef struct rtree RTree;

uint8_t RTreeInsert(struct rtree *rtree, uint64_t *rect, void *item);
struct rtree *RTreeInit(size_t elsize, int dims);
void RTreeDestroy(struct rtree *rtree);
size_t RTreeCount(struct rtree *rtree);
uint8_t RTreeDelete(struct rtree *rtree, uint64_t *rect, void *item);
uint8_t RTreeSearch(struct rtree *rtree, uint64_t *rect, 
                  uint8_t (*iter)(const uint64_t *rect, const void *item, 
                               void *udata), 
                  void *udata);

void RTreeSetAllocator(void *(malloc)(size_t), void (*free)(void*));

#endif

