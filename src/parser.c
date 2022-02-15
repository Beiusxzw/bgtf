#include "bgtf.h"

struct GenomicLocation *parseGenomicLocation(char *src)
{
    struct GenomicLocation *loc = malloc(sizeof(struct GenomicLocation));
    char *chrom = NULL;
    uint64_t start,end;
    char *endtok;
    char *tok = strtok_r(src, ":", &endtok);
    while (tok != NULL) {
        chrom = tok;
        tok = strtok_r(NULL, "-", &endtok);
        start = atoll(tok);
        end = atoll(endtok);
        break;
    }
    if (!chrom) return NULL;
    loc->chrom = chrom; loc->start = start; loc->end = end;
    return loc;
}

static uint8_t _parseRecordAttrKV(RecordAttrKV_t *attr, char *src)
{
    char *endtok;
    char *tok = strtok_r(src, " ", &endtok);
    char *k,*v;
    uint8_t c = 0;
    while (tok != NULL) {
        if (c > 0) {
            return 1;
        }
        k = tok;
        v = tok = strtok_r(NULL, " ", &endtok);
        int vlen = strlen(v);
        v[vlen-1] = '\0';
        v = v + 1;
        tok = strtok_r(NULL, " ", &endtok);
        c++;
    }
    setRecordAttrKV(attr, k, v);
}

RecordAttrKV_t *parseRecordAttrKV(char *src)
{
    char *endtok;
    RecordAttrKV_t *attr = (RecordAttrKV_t *)malloc(sizeof(RecordAttrKV_t));
    memset(attr, 0, sizeof(RecordAttrKV_t));
    char *tok = strtok_r(src, ";", &endtok);
    while (tok != NULL) {
        if (!_parseRecordAttrKV(attr, tok)) {
            // warn("Record is corrupted");
            // This is because some attributes from GTF 
            // contains multiple spaces.
        };
        tok = strtok_r(NULL, ";", &endtok);
    }
    return attr;
}

BGTFRecord_t *parseOneRecord(char *buf, char *line) {
    if (*line == '#') return NULL;
    BGTFRecord_t *record = (BGTFRecord_t *)malloc(sizeof(BGTFRecord_t));
    int len = strlen(line);
    char *c = line;
    char *ptr = &buf[0];
    for (size_t i = 0; i < len; i++) {
        if (*c == '\t' || *c == '\n') {
            *ptr = '\0';
            ptr ++;
            c++;
        } else {
            *ptr = *c;
            ptr++;
            c++;
        }
    }
    ptr = &buf[0];
    record->chrom = malloc(strlen(ptr)+1);
    memset(record->chrom, 0, strlen(ptr)+1);
    memcpy(record->chrom, ptr, strlen(ptr));
    ptr += strlen(ptr) + 1;
    record->source = malloc(strlen(ptr)+1);
    memset(record->source, 0, strlen(ptr)+1);
    memcpy(record->source, ptr, strlen(ptr));
    ptr += strlen(ptr) + 1;
    record->type = malloc(strlen(ptr)+1);
    memset(record->type, 0, strlen(ptr)+1);
    memcpy(record->type, ptr, strlen(ptr));
    ptr += strlen(ptr) + 1;
    record->start = atoll(ptr);
    ptr += strlen(ptr) + 1;
    record->end = atoll(ptr);
    ptr += strlen(ptr) + 1;
    record->score = 0; // TODO: fill this field
    ptr += strlen(ptr) + 1;
    record->strand = 1 ? (*ptr == '+') : 0;
    ptr += strlen(ptr) + 1;
    record->phase = 0; // TODO: fill this field
    ptr += strlen(ptr) + 1;
    record->attrs = parseRecordAttrKV(ptr);
    return record;
}
