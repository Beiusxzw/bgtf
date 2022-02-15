#include "utils.h"
#include "bgtf.h"
#include "bgtfIO.h"
#include "rtree.h"
#include "molcount.h"

static void usage(FILE *fp)
{
    fprintf(fp,    
"Usage:   bgtf <command> [options]\n"
"\n"
"Commands:\n"
"     convert          create a BGTF File from GTF file\n"
"     count            count \n"
"For detailed usage of subcommand, use bgtf <command> -h\n"
"For detailed version of subcommand, use bgtf <command> -v\n"
"\n");
}

int convert_main(int argc, char *argv[]);
int count_main(int argc, char *argv[]);

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
    char src[] = "1:3002016-4102125";
    struct GenomicPosition *loc = parseGenomicLocation(src);
    BGTF *fp = GTFRead("/Users/snowxue/Documents/refData/Mus_musculus.GRCm38.97.gtf");
    BGTFSave(fp, "/Users/Snowxue/Downloads/mm10.bgtf", "wb+");
    BGTFGetRecord(fp, loc, BGTFprint, NULL);
    BGTFClose(fp);
*/
    if (argc < 2) { usage(stderr); return 1; }
    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { usage(stdout); return 0; }
    }
    int ret = 0;
    if (strcmp(argv[1], "convert") == 0) ret = convert_main(argc-1, argv+1);
    if (strcmp(argv[1], "count") == 0) ret = count_main(argc-1, argv+1);
    return ret;
}
