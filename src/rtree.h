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



typedef struct RTree RTree_t;

uint8_t RTreeInsert(RTree_t *rtree, uint64_t *rect, void *item);
RTree_t *RTreeInit(size_t elsize, int dims);
void RTreeDestroy(RTree_t *rtree);
size_t RTreeCount(RTree_t *rtree);
uint8_t RTreeDelete(RTree_t *rtree, uint64_t *rect, void *item);
uint8_t RTreeSearch(
    RTree_t *rtree, 
    uint64_t *rect, 
    uint8_t (*iter)(
        const uint64_t *rect, 
        const void *item, 
        void *udata
    ), 
    void *udata);

void RTreeSetAllocator(void *(malloc)(size_t), void (*free)(void*));

#endif

