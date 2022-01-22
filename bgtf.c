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
    hdr->version = LIBBIGWIG_VERSION;
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

uint8_t BGTFBuildIndex(BGTF *file)
{
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
            // logf("%s", curr->chrom);
            idx = BGTFRecordIdxInit();
            HashTablePut(file->index->chrom_idx, curr->chrom, idx);
            head = curr;
            last = curr->chrom;
        }
        idx = HashTableGet(file->index->chrom_idx, last);
        RTreeInsert(idx->rt, (uint64_t[]){ curr->start, curr->end}, curr);
        c++;
    })
    // finish up
    idx = HashTableGet(file->index->chrom_idx, last);
    idx->head = head;
    idx->tail = curr;
    file->n_records = HashTableSize(file->index->chrom_idx);
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
    file->index = BGTFIdxInit();
    file->hdr = BGTFHdrInit();
    file->write = NULL;
    file->is_write = 0;
    return file;
}

void BGTFClose(BGTF *fp)
{
    ArrayListDestroy(fp->records);
    BGTFIdxDestroy(fp->index);
    BGTFHdrDestroy(fp->hdr);
    free(fp);
}

int8_t RecordCompare(BGTFRecord_t *record1, BGTFRecord_t *record2)
{
    int8_t r;
    r = strcmp(record1->chrom, record2->chrom);
    if (r == 0) goto compare_pos;
    return r;
compare_pos:
    if (record1->start == record2->end) 
        return 0;
    else return record1->start > record1->end ? 1 : -1;
}

void BGTFSortRecord(BGTF *fp)
{
    log("Sorting the GTF file in memory");
    ArrayList *records = fp->records;
    ArrayListSort(records, RecordCompare);
}


void BGTFGetRecord(BGTF *fp, struct GenomicLocation *loc, void (*func)(BGTFRecord_t *record))
{
    char *chrom = loc->chrom;
    uint64_t start = loc->start;
    uint64_t end = loc->end;
    BGTFRecordIdx *idx = HashTableGet(fp->index->chrom_idx, chrom);
    uint64_t interval[] = {start, end};
    RTreeSearch(idx->rt, interval, func, NULL);
}