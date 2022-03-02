#include <getopt.h>

#include "bgtf.h"
#include "bgtfIO.h"

uint8_t setRecordAttrKV_helper(RecordAttrKV_t *attr, uint32_t magic, char *v)
{
    switch (magic)
    {
    case GENE_ID:
        attr->gene_id = malloc(strlen(v)+1);
        memset(attr->gene_id, 0, strlen(v)+1);
        memcpy(attr->gene_id, v, strlen(v));
        break;
    case GENE_VERSION:
        attr->gene_version = atol(v);
        break;
    case GENE_NAME:
        attr->gene_name = malloc(strlen(v)+1);
        memset(attr->gene_name, 0, strlen(v)+1);
        memcpy(attr->gene_name, v, strlen(v));
        break;
    case TRANSCRIPT_ID:
        attr->transcript_id = malloc(strlen(v)+1);
        memset(attr->transcript_id, 0, strlen(v)+1);
        memcpy(attr->transcript_id, v, strlen(v));
        break;
    case TRANSCRIPT_VERSION:
        attr->transcript_version = atol(v);
        break;
    case TRANSCRIPT_NAME:
        attr->transcript_name = malloc(strlen(v)+1);
        memset(attr->transcript_name, 0, strlen(v)+1);
        memcpy(attr->transcript_name, v, strlen(v));
        break;
    case EXON_NUMBER:
        attr->exon_number = atol(v);
        break;
    case EXON_ID:
        attr->exon_id = malloc(strlen(v)+1);
        memset(attr->exon_id, 0, strlen(v)+1);
        memcpy(attr->exon_id, v, strlen(v));
        break;
    case FAMILY_ID:
        attr->family_id = malloc(strlen(v)+1);
        memset(attr->family_id, 0, strlen(v)+1);
        memcpy(attr->family_id, v, strlen(v));
        break;
    case CLASS_ID:
        attr->class_id = malloc(strlen(v)+1);
        memset(attr->class_id, 0, strlen(v)+1);
        memcpy(attr->class_id, v, strlen(v));
        break;
    default:
        break;
    }
    return 0;
}


uint8_t setRecordAttrKV(RecordAttrKV_t *attr, char *k, char *v)
{
    enum BGTF_ATTR_MAGIC magic = HashTableStringHashFunction(k);
    return setRecordAttrKV_helper(attr, magic, v);
}

char *getRecordAttrKV_helper(RecordAttrKV_t *attr, uint32_t magic)
{
    switch (magic)
    {
    case GENE_ID:
        return attr->gene_id;
    case GENE_NAME:
        return attr->gene_name;
    case TRANSCRIPT_ID:
        return attr->transcript_id;
    case TRANSCRIPT_NAME:
        return attr->transcript_name;
    case EXON_ID:
        return attr->exon_id;
    case FAMILY_ID:
        return attr->family_id;
    case CLASS_ID:
        return attr->class_id;
    default:
        break;
    }
    return NULL;
}

char *getRecordAttrKV(RecordAttrKV_t *attr, char *k)
{
    enum BGTF_ATTR_MAGIC magic = HashTableStringHashFunction(k);
    return getRecordAttrKV_helper(attr, magic);
}


void RecordAttrKVDestory(RecordAttrKV_t *attr)
{
    if (attr->gene_id) free(attr->gene_id);
    if (attr->gene_name) free(attr->gene_name);
    if (attr->transcript_id) free(attr->transcript_id);
    if (attr->transcript_name) free(attr->transcript_name);
    if (attr->exon_id) free(attr->exon_id);
    if (attr->family_id) free(attr->family_id);
    if (attr->class_id) free(attr->class_id);
    free(attr);
}

void BGTFRecordDestroy(BGTFRecord_t *record)
{
    RecordAttrKVDestory(record->attrs);
    free(record->chrom);
    free(record->source);
    free(record->type);
    free(record);
}

BGTFHdr_t *BGTFHdrInit()
{
    BGTFHdr_t *hdr = (BGTFHdr_t *)malloc(sizeof(BGTFHdr_t));
    hdr->version = LIBBGTF_VERSION;
    return hdr;
}

void BGTFHdrDestroy(BGTFHdr_t *hdr)
{
    free(hdr);
}

BGTFIdx *BGTFIdxInit()
{
    BGTFIdx *index = (BGTFIdx *)malloc(sizeof(BGTFIdx));
    index->chrom_idx = HashTableCreate(32);
    HashTableSetHashFunction(index->chrom_idx, HashTableStringHashFunction);
    HashTableSetKeyComparisonFunction(index->chrom_idx, strcmp);
    return index;
}

static BGTFRecordIdx *BGTFRecordIdxInit()
{
    BGTFRecordIdx *idx = (BGTFRecordIdx *)malloc(sizeof(BGTFRecordIdx));
    memset(idx, 0, sizeof(BGTFRecordIdx));
    idx->rt = RTreeInit(sizeof(BGTFRecord_t), 1);
    return idx;
}

static void BGTFRecordIdxDestroy(BGTFRecordIdx *idx)
{
    RTreeDestroy(idx->rt);
    free(idx);
}

uint8_t BGTFBuildIndex(BGTF *file, char *attr_name)
{
    uint32_t attr_idx = 1; char *k;
    if (attr_name) {
        file->attr_table = HashTableCreate(0x1000);
        HashTableSetKeyComparisonFunction(file->attr_table, strcmp);
        HashTableSetHashFunction(file->attr_table, HashTableStringHashFunction);
    }
    BGTFRecord_t *curr = NULL, *head = NULL;
    char *last = NULL; size_t c = 0;
    BGTFRecordIdx *idx;
    log("Generating in-memory chromosome index");
    ArrayListForEach(file->records, curr, {
        if (HashTableGet(file->index->chrom_idx, curr->chrom) == NULL) {
            if (c > 0) {
                idx = HashTableGet(file->index->chrom_idx, last);
                idx->head = head;
                idx->tail = curr;
            }
            idx = BGTFRecordIdxInit();
            HashTablePut(file->index->chrom_idx, curr->chrom, idx);
            head = curr;
            last = curr->chrom;
        }
        if (attr_name) {
            k = getRecordAttrKV(curr->attrs, attr_name);
            if (!(HashTableContainsKey(file->attr_table, k))) {
                HashTablePut(file->attr_table, k, attr_idx);
                attr_idx++;
            }
        }
        idx = HashTableGet(file->index->chrom_idx, last);
        RTreeInsert(idx->rt, (uint64_t[]){ curr->start, curr->end}, curr);
        c++;
    })
    if (attr_name) {
        ArrayList *sort_arr = ArrayListCreate(HashTableSize(file->attr_table));
        HashTableForEach(file->attr_table, k, attr_idx, {
            ArrayListPush(sort_arr, k);
        })
        ArrayListSort(sort_arr, strcmp);
        attr_idx = 1;
        ArrayListForEach(sort_arr, k, {
            HashTablePut(file->attr_table, k, attr_idx);
            attr_idx++;
        })
        ArrayListDestroy(sort_arr);
    }
    // finish up
    idx = HashTableGet(file->index->chrom_idx, last);
    idx->head = head;
    idx->tail = curr;
    file->n_records = file->records->numOfElements;
    file->n_chroms = HashTableSize(file->index->chrom_idx);
    if (attr_name) {
        logf("Total number of %s: %d", attr_name, HashTableSize(file->attr_table));
    }
}

void BGTFIdxDestroy(BGTFIdx *idx)
{
    char *chrom;
    BGTFRecordIdx *value;
    HashTableForEach(idx->chrom_idx, chrom, value, {
        BGTFRecordIdxDestroy(value);
    })
    HashTableDestroy(idx->chrom_idx);
}

BGTF *BGTFInit()
{
    BGTF *file = (BGTF *)malloc(sizeof(BGTF));
    memset(file, 0, sizeof(BGTF));
    file->index = BGTFIdxInit();
    file->hdr = BGTFHdrInit();
    file->write = NULL;
    file->is_write = 0;
    return file;
}

void BGTFClose(BGTF *fp)
{
    if (fp->attr_table) HashTableDestroy(fp->attr_table);
    ArrayListDestroy(fp->records);
    BGTFIdxDestroy(fp->index);
    BGTFHdrDestroy(fp->hdr);
    free(fp);
}

int8_t RecordCompare(BGTFRecord_t *a, BGTFRecord_t *b)
{
    int8_t r;
    r = strcmp(a->chrom, b->chrom);
    if (r == 0) goto compare_pos;
    return r;
compare_pos:
    if ((a->start == b->start) && (a->end == b->end)) return 0;
    if (a->start == b->start) {
        return a->end > b->end ? 1 : -1;
    } else {
        return a->start > b->start ? 1 : -1;
    }
}

void BGTFSortRecord(BGTF *fp)
{
    log("Sorting the GTF file in memory");
    ArrayList *records = fp->records;
    ArrayListSort(records, RecordCompare);
}


uint8_t BGTFGetRecord(BGTF *fp, struct GenomicLocation *loc, void (*func)(BGTFRecord_t *record), void *ustruct)
{
    char *chrom = loc->chrom;
    uint64_t start = loc->start;
    uint64_t end = loc->end;
    if (chrom == NULL) return 1;
    BGTFRecordIdx *idx = HashTableGet(fp->index->chrom_idx, chrom);
    if (idx == NULL) return 1;
    uint64_t interval[] = {start, end};
    RTreeSearch(idx->rt, interval, func, ustruct);
    return 0;
}

void destroyGenomicLocation(struct GenomicLocation *loc, void(*dealloc_data)(char *data))
{
    if (dealloc_data && loc->data) {
        dealloc_data(loc->data);

    } else {
        free(loc->data);
    }
    free(loc);
}


/* Command-line functions */
void convert_usage()
{
    fprintf(stdout,    
"Usage:   bgtf count [options]\n"
"\n"
"Arguments:\n"
"   Required arguments:\n"
"     -g          [PATH]   Input GTF file \n"
"     -o          [PATH]   Output BGTF file  \n"
"\n");
    return;
}
void convert_version()
{
    fprintf(stdout,    
"bgtf convert version %s\n", BGTF_CONVERT_VERSION
    );
    return;
}

int convert_main(int argc, char *argv[])
{
    int c;
    char *file_path = NULL;
    char *output_path = NULL;
    while (1)
    {
        static struct option long_options[] =
            {
                {"gtf", required_argument, 0, 'g'},
                {"output", required_argument, 0, 'o'},
                {"help", no_argument, NULL, 'h'},
                {"version", no_argument, NULL, 'v'},
                {0, 0, 0, 0}};
        /* getopt_long stores the option index here. */
        int option_index = 0;
        c = getopt_long(argc, argv, "g:ho:v", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf("option %s", long_options[option_index].name);
            if (optarg)
                printf(" with arg %s", optarg);
            printf("\n");
            break;

        case 'g':
            file_path = optarg;
            break;  
        
        case 'o':
            output_path = optarg;
            break;

        case 'h':
            convert_usage();
            exit(0);

        case 'v':
            convert_version();
            exit(0);
        

        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            fatal("Argument capture failed\n");
        }
    }
    BGTF *fp = GTFRead(file_path, NULL);
    BGTFSave(fp, output_path, "wb+");
    BGTFClose(fp);
}