/*

  Copyright (c) 2016 Martin Sustrik

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"),
  to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom
  the Software is furnished to do so, subject to the following conditions:
  The above copyright notice and this permission notice shall be included
  in all copies or substantial portions of the Software.
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
  IN THE SOFTWARE.

*/

#ifndef RBTREE_INCLUDED
#define RBTREE_INCLUDED

#include <stddef.h>
#include <stdint.h>

struct rbtree_node {
    int red;
    struct rbtree_node *left;
    struct rbtree_node *right;
    struct rbtree_node *up;
    int64_t val;
};

struct rbtree {
    struct rbtree_node root;
    struct rbtree_node nil;
};

/* Initialize the tree. */
void rbtree_init(struct rbtree *self);

/* Returns 1 if there are no items in the tree. 0 otherwise. */
int rbtree_empty(struct rbtree *self);

/* Insert an item into the tree & set its value to 'val' */
void rbtree_insert(struct rbtree *self, int64_t val,
    struct rbtree_node *item);

/* Remove an item from a tree. */
void rbtree_erase(struct rbtree *self, struct rbtree_node *item);

/* Return an item with the lowest value. If there are no items in the tree, NULL
   is returned. */
struct rbtree_node *rbtree_first(struct rbtree *self);

/* Iterate through the tree. Items are returned starting with those with
   the lowest values and ending with those with the highest values. Items with
   equal values are returned in no particular order. If 'it' points to the
   last item, NULL is returned. */
struct rbtree_node *rbtree_next(struct rbtree *self,
    struct rbtree_node *it);

#endif
