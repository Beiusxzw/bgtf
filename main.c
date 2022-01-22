#include "utils.h"
#include "bgtf.h"
#include "bgtfIO.h"
#include "rtree.h"

static void map(void *key, void *value)
{
    logf("%s-%s", key, value);
}

static void BGTFprint(const uint64_t *rect, void *item, void *udata)
{
    BGTFRecord_t *record = (BGTFRecord_t *)item;
    BGTFfprint(stdout, record, "gene_name");
}


int main(int argc, char *argv[])
{
#ifdef SIGSEGV
    signal(SIGSEGV, crash_handler);
#endif /* SIGSEGV */
#ifdef SIGFPE
    signal(SIGFPE, crash_handler);
#endif /* SIGFPE */
#ifdef SIGILL
    signal(SIGILL, crash_handler);
#endif /* SIGILL */
#ifdef SIGPIPE
    signal(SIGPIPE, crash_handler);
#endif /* SIGPIPE */
#ifdef SIGBUS
    signal(SIGBUS, crash_handler);
#endif /* SIGBUS */

    /*
    BGTF *fp = GTFRead(argv[1]);
    BGTFSave(fp, "/Users/snow/Downloads/test.bgtf", "wb+");
    char src[] = "8:47839148-47839246";
    struct GenomicPosition *loc = parseGenomicLocation(src);
    BGTFGetRecord(fp, loc, BGTFprint);
    free(loc);
    BGTFClose(fp);
    */


    BGTF *rfp = BGTFInit();
    BGTFLoad(rfp, "/Users/snow/Downloads/test.bgtf");
    BGTFflushRecord(stdout, rfp, "gene_id");


    BGTFClose(rfp);
    // Testing redblack tree;
    char *key, *val;
    RedBlackTree *rbtree =  RedBlackTreeInit();
    RedBlackTreeSetKeyCompare(rbtree, strcmp);
    RedBlackTreeTraversal(rbtree, map);
    RedBlackTreeDestroy(rbtree);

    return 0;
}
