#include "list.h"

void quickSortRun(void *arr, int spot_low, int spot_high, int compare(void *arr, int l, int r), void exchange(void *arr, int l, int r));

void quickSort(void *arr, int arr_size, int compare(void *arr, int l, int r), void exchange(void *arr, int l, int r))
{
    quickSortRun(arr, 0, arr_size - 1, compare, exchange);
}

void quickSortRun(void *arr, int spot_low, int spot_high, int compare(void *arr, int l, int r), void exchange(void *arr, int l, int r))
{
    int pivot, j, i;

    if (spot_high <= spot_low)
        return;
    pivot = spot_high;
    i = spot_low;

    for (j = spot_low; j < spot_high; j++)
        if (compare(arr, j, pivot) <= 0)
        {
            exchange(arr, i, j);
            i++;
        }

    exchange(arr, i, spot_high);

    quickSortRun(arr, spot_low, i - 1, compare, exchange);
    quickSortRun(arr, i + 1, spot_high, compare, exchange);
}
void basicSort_run(void *arr, int start, int items, int compare(void *arr, int l, int r), void exchange(void *arr, int l, int r))
{
    int i, j;
    for (i = start; i < start + items - 1; i++)
    {
        int min_j = i;
        for (j = i + 1; j < start + items; j++)
        {
            if (compare(arr, min_j, j) > 0)
                min_j = j;
        }
        if (i != min_j)
            exchange(arr, i, min_j);
    }
}

void mergeSortRun(void *arr, int start, int items, int compare(void *arr, int l, int r), void exchange(void *arr, int l, int r), void merge(void *arr, int start, int items, int items2))
{
    if (items > 11)
    {
        int half_point = items / 2;
        mergeSortRun(arr, start, half_point, compare, exchange, merge);
        mergeSortRun(arr, start + half_point, items - half_point, compare, exchange, merge);
        merge(arr, start, half_point, items - half_point);
    }
    else
    {
        basicSort_run(arr, start, items, compare, exchange);
    }
}
void mergeSort(void *arr, int arr_size, int compare(void *arr, int l, int r), void exchange(void *arr, int l, int r), void merge(void *arr, int start, int items, int items2))
{
    mergeSortRun(arr, 0, arr_size, compare, exchange, merge);
}

ArrayList *ArrayListCreate(int init_capacity)
{
    ArrayList *ret = malloc(sizeof(ArrayList));
    memset(ret, 0, sizeof(ArrayList));
    ret->capacityOfElements = init_capacity;
    ret->elementList = malloc(sizeof(void *) * init_capacity);
    return ret;
}

void ArrayListDestroy(ArrayList *list)
{
    size_t x1;
    if (list->elemDeallocator)
        for (x1 = 0; x1 < list->numOfElements; x1++)
            list->elemDeallocator(list->elementList[x1]);

    free(list->elementList);
    free(list);
}

void *ArrayListGet(ArrayList *list, size_t n)
{
    if (n < 0 || n >= list->numOfElements)
        return NULL;
    return list->elementList[n];
}

int ArrayListPush_NoRepeatedPtr(ArrayList *list, void *new_elem)
{
    size_t ii;
    for (ii = 0; ii < list->numOfElements; ii++)
    {
        if (list->elementList[ii] == new_elem)
            return -1;
    }
    return ArrayListPush(list, new_elem);
}

int ArrayListPush(ArrayList *list, void *new_elem)
{
    if (list->capacityOfElements <= list->numOfElements)
    {
        if (list->capacityOfElements * 1.3 > list->capacityOfElements + 10)
            list->capacityOfElements = list->capacityOfElements * 1.3;
        else
            list->capacityOfElements = list->capacityOfElements + 10;

        list->elementList = realloc(list->elementList, sizeof(void *) * list->capacityOfElements);
        assert(list->elementList);
    }
    list->elementList[list->numOfElements++] = new_elem;
    return list->numOfElements;
}
void ArrayListSetDeallocationFunction(ArrayList *list, void (*elem_deallocator)(void *elem))
{
    list->elemDeallocator = elem_deallocator;
}

int ArrayListSort_compare(void *sortdata0, int L, int R)
{
    void **sortdata = sortdata0;
    ArrayList *list = sortdata[0];
    int (*comp_elems)(void *L_elem, void *R_elem) = sortdata[1];
	void *L_elem = ArrayListGet(list, L);
    void *R_elem = ArrayListGet(list, R);
    return comp_elems(L_elem, R_elem);
}

void ArrayListSort_exchange(void *sortdata0, int L, int R)
{
    void **sortdata = sortdata0;
    ArrayList *list = sortdata[0];

    void *tmpp = list->elementList[L];
    list->elementList[L] = list->elementList[R];
    list->elementList[R] = tmpp;
}

void ArrayListSort_merge(void *sortdata0, int start, int items, int items2)
{
    void **sortdata = sortdata0;
    ArrayList *list = sortdata[0];
    int (*comp_elems)(void *L_elem, void *R_elem) = sortdata[1];

    void **merged = malloc(sizeof(void *) * (items + items2));
    int write_cursor, read1 = start, read2 = start + items;

    for (write_cursor = 0; write_cursor < items + items2; write_cursor++)
    {
        void *Elm1 = list->elementList[read1];
        void *Elm2 = list->elementList[read2];

        int select_1 = (read1 == start + items) ? 0 : (read2 == start + items + items2 || comp_elems(Elm1, Elm2) < 0);
        if (select_1)
            merged[write_cursor] = list->elementList[read1++];
        else
            merged[write_cursor] = list->elementList[read2++];
    }

    memcpy(list->elementList + start, merged, sizeof(void *) * (items + items2));
    free(merged);
}

void ArrayListSort(ArrayList *list, int compare_L_minus_R(void *L_elem, void *R_elem))
{
    void *sortdata[2];
    sortdata[0] = list;
    sortdata[1] = compare_L_minus_R;

    mergeSort(sortdata, list->numOfElements, ArrayListSort_compare, ArrayListSort_exchange, ArrayListSort_merge);
}

LinkedList *LinkedListInit()
{
    LinkedList *list = malloc(sizeof(LinkedList));
    memset(list, 0, sizeof(LinkedList));
    return list;
}

uint8_t __LinkedListAppend(LinkedList *list, void *elem, uint64_t offset)
{
    void *curr, **next;
    if (list->elem_head == NULL) {
        list->elem_head = list->elem_tail = elem;
        list->sz = 1;
        return 0;
    } else {
        next = (list->elem_tail + offset);
        *next = elem;
        list->elem_tail = elem;
        list->sz++;
        return 0;
    }
    return 1;
}

uint8_t __LinkedListRemove(LinkedList *list, void *elem, uint64_t offset)
{
    void **curr = &list->elem_head;
    while ((*curr) != elem) {
        curr = (void **)(*curr + offset);
    }
    *curr = (*(void **)(elem + offset));
}


void __LinkedListDestroy(LinkedList *list, uint64_t offset, uint8_t free_elem)
{
    if (free_elem && (list->sz > 0)) {
        if (!list->dealloc) {
            warnf("No deallocation function found. Skip")
        } else {
            void *curr = list->elem_head, *next;
            while (curr != list->elem_tail) {
                next = (*(void **)(curr + offset));
                list->dealloc(curr);
                    
                curr = next;
            }
            list->dealloc(curr);
        }
    }
    free(list);
}

void *__LinkedListGet(LinkedList *list, size_t idx, uint64_t offset) 
{
    void *curr = list->elem_head;
    size_t i = 0;
    while (i < idx) {
        curr = (*(void **)(curr + offset));
        i++;
    }
    return curr;
}

void LinkedListSetDealloc(LinkedList *list, void (*dealloc)(void *elem))
{
    list->dealloc = dealloc;
}

static void *LinkedListMerge(void *a, void *b, int (*compare)(void *, void *), uint64_t offset)
{
    void *result = NULL;
    void **next;
    if (a == NULL) return b;
    else if (b == NULL) return a;
    if (compare(a, b) <= 0) {
        result = a;
        next = (void **)(result + offset);
        *next = LinkedListMerge((*(void **)(a + offset)), b, compare, offset);
    } else {
        result = b;
        next = (void **)(result + offset);
        *next = LinkedListMerge(a, ((*(void **)(b + offset))), compare, offset);
    }
    return result;
}

static void LinkedListSplit(void *source, void **front_ref, void **back_ref, uint64_t offset)
{
    void *a, *b;
    void **next;
    a = source;
    b = (*(void **)(a + offset));
    while (b != NULL) {
        b = (*(void **)(b + offset));
        if (b != NULL) {
            a = (*(void **)(a + offset));
            b = (*(void **)(b + offset));
        }
    }
    *front_ref = source;
    *back_ref = (*(void **)(a + offset));
    next = (void **)(a + offset);
    *next = NULL;
}

static void __LinkedListMergeSort_helper(void **head_ref, int (*compare)(void *, void *), uint64_t offset)
{
    void *head = *head_ref;
    void *a, *b;
    if ((head == NULL) || ((*(void **)(head + offset)) == NULL)) {
        return;
    }
    LinkedListSplit(head, &a, &b, offset);
    __LinkedListMergeSort_helper(&a, compare, offset);
    __LinkedListMergeSort_helper(&b, compare, offset);
    *head_ref = LinkedListMerge(a, b, compare, offset);
}

void __LinkedListMergeSort(LinkedList *list, int (*compare)(void *, void *), uint64_t offset)
{
    __LinkedListMergeSort_helper(&list->elem_head, compare, offset);
    void *curr = list->elem_head;
    for (size_t i = 0; i < list->sz - 1; ++i) {
        curr = (*(void **)(curr + offset));
    }
    list->elem_tail = curr;
}

/*
typedef struct IntElem {
    int data;
    struct IntElem *next;
} IntElem_t;

uint8_t IntElemCmp(IntElem_t *a, IntElem_t *b)
{
    if (a->data == b->data) return 0;
    else return a->data > b->data ? 1 : 0;
}

int main(int argc, char const *argv[])
{

    LinkedList *list = LinkedListInit();
    IntElem_t *iter, *rm;
    for (int i = 10; i > 0; --i) {
        IntElem_t *item = malloc(sizeof(IntElem_t));
        item->data = i;
        if (i == 5) rm = item;
        LinkedListAppend(struct IntElem, list, item);
    }
    for (int i = 20; i > 10; --i) {
        IntElem_t *item = malloc(sizeof(IntElem_t));
        item->data = i;
        if (i == 5) rm = item;
        LinkedListAppend(struct IntElem, list, item);
    }
    LinkedListForEach(IntElem_t, list, iter, {
        printf("%d\n", iter->data);
    })
    LinkedListRemove(IntElem_t, list, rm);
    LinkedListMergeSort(IntElem_t, list, IntElemCmp);
    IntElem_t *item =  LinkedListGet(IntElem_t, list, 3);
    printf("array[3]: %d\n", item->data);
    LinkedListForEach(IntElem_t, list, iter, {
        printf("%d\n", iter->data);
    })
    return 0;
}
*/