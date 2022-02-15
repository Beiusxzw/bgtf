#ifndef BGTF_LIST_H
#define BGTF_LIST_H

#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include "utils.h"


/* These structs should not be accessed directly from user code.
 * All access should be via the public functions declared below. */


/**
 * @brief ArrayList Data structure
 * 
 */
typedef struct
{
    void **elementList;
    size_t numOfElements;
    size_t capacityOfElements;
    void (*elemDeallocator)(void *elem);
} ArrayList;

/**
 * @brief Initialization of an array list
 * 
 * @param init_capacity the initial capacity
 * @return ArrayList* 
 */
ArrayList *ArrayListCreate(int init_capacity);

/**
 * @brief Destory an array list
 * 
 * @param list 
 */
void ArrayListDestroy(ArrayList *list);

/**
 * @brief Get item from an array list by index
 * 
 * @param list 
 * @param n 
 * @return void* 
 */
void *ArrayListGet(ArrayList *list, size_t n);

/**
 * @brief Add item to an array list
 * 
 * @param list 
 * @param new_elem 
 * @return int 
 */
int ArrayListPush(ArrayList *list, void *new_elem);
int ArrayListPush_NoRepeatedPtr(ArrayList *list, void *new_elem);

/**
 * @brief Deallocator function for the array list 
 * 
 * @param list 
 * @param elem_deallocator 
 */
void ArrayListSetDeallocationFunction(ArrayList *list, void (*elem_deallocator)(void *elem));
void ArrayListSort(ArrayList *list, int compare_L_minus_R(void *L_elem, void *R_elem));

#define ArrayListForEach(list, elem, code) \
{ \
	for (size_t i = 0; i < list->numOfElements; ++i) { \
		elem = (list->elementList[i]); \
		code; \
	} \
} while (0); \

#define ArrayListPartition(list, n_partition, elem_ptr, partition_size, code) \
{ \
    size_t psz = list->numOfElements / n_partition; \
   	for (size_t i = 0; i < list->numOfElements; i += psz) { \
        if (i + psz <= list->numOfElements) { \
            partition_size =  psz; \
        } else { \
            partition_size = list->numOfElements - i; \
        } \
		elem_ptr = (list->elementList + i); \
		code; \
	} \ 
} while (0); \

/**
 * @brief General LinkedList data structure
 *        Elements in the linked-list should be defined like
 *        struct Element {
 *            typename data;
 *            struct Element *next; // the `next` member can be defined 
 *                                  // anywhere in the structure
 *        }
 *        
 */
typedef struct {
    void *elem_head;
    void *elem_tail;
    size_t sz;
    void (*dealloc)(void *elem);
} LinkedList;

LinkedList *LinkedListInit();

/**
 * @brief Implementation of appending an item to the linked-list
 *        This function is wrapped by a macro function.
 * 
 * @param list The linked-list
 * @param elem element to append
 * @param offset offset to the `next` member in the struct
 * @return uint8_t 0 on success and 1 on failure
 */
uint8_t __LinkedListAppend(LinkedList *list, void *elem, uint64_t offset);

uint8_t __LinkedListRemove(LinkedList *list, void *elem, uint64_t offset);
void __LinkedListDestroy(LinkedList *list, uint64_t offset, uint8_t free_elem);
void *__LinkedListGet(LinkedList *list, size_t idx, uint64_t offset);
void LinkedListSetDealloc(LinkedList *list, void (*dealloc)(void *elem));
void __LinkedListMergeSort(LinkedList *list, int (*compare)(void *, void *), uint64_t offset);
#define LinkedListAppend(typename, list, elem) { \
    uint64_t offset = offsetof(typename, next); \
    if ((__LinkedListAppend(list, elem, offset))) { \
        fatalf("LinkedList at <0x%x> append failed", list); \
    }; \
} while (0); \

#define LinkedListRemove(typename, list, elem) ({ \
    void *ret; \
    do { \
        uint64_t offset = offsetof(typename, next); \
        ret = __LinkedListRemove(list, elem, offset); \
    } while (0); \
    ret; \
}) \

#define LinkedListDestroy(typename, list, free_elem) { \
    uint64_t offset = offsetof(typename, next); \
    __LinkedListDestroy(list, offset, free_elem); \
} while (0); \

#define LinkedListGet(typename, list, idx) ({ \
    void *ret; \
    do { \
    uint64_t offset = offsetof(typename, next); \
    ret = (void *)__LinkedListGet(list, idx, offset); \
    } while (0); \
    ret; \
}) \

#define LinkedListMergeSort(typename, list, compare)  { \
    uint64_t offset = offsetof(typename, next); \
    __LinkedListMergeSort(list, compare, offset); \
} while (0); \

#define LinkedListForEach(typename, list, vvar, code) { \
    uint64_t offset = offsetof(typename, next); \
    vvar = (typename *) list->elem_head; \
    while (vvar) { \
        code; \
        vvar = (typename *)(*(void **)((void *)vvar + offset)); \
    } \
} \

#define LinkedListForEachInRange(typename, list, vvar, from, to, code) { \
    uint64_t offset = offsetof(typename, next); \
    vvar = (typename *) list->elem_head; \
    while (vvar) { \
        code; \
        vvar = (typename *)(*(void **)((void *)vvar + offset)); \
    } \
} \

#define LinkedListNext(typename, vvar) ((typename *)(*(void **)((void *)vvar + offsetof(typename, vvar))))

#endif