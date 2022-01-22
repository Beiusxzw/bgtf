#ifndef LIBBGTF_H
#define LIBBGTF_H

#include "hashtable.h"
#include "utils.h"
#include "rtree.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <zlib.h>

#ifdef __cplusplus
extern "C" {
#endif


/*!
 * The library version number.
 */
#define LIBBIGWIG_VERSION 0x000001;
#define BGTF_MAGIC 0x8518B4

enum RTreeMode {
    RTREE_DISK = 0x1,
    RTREE_MEM  = 0x2
};

/**
 * @brief The header section of a binary GTF file.
 *
 */
typedef struct {
    uint64_t version;         /* The version of the file */
    uint16_t n_levels;        /* number of zoom levels */
    uint64_t cl_offset;       /* offset to the on-disk chromosome list */
    uint64_t record_offset;   /* offset to the on-disk GTF record */
    uint64_t index_offset;    /* offset to the on-disk data index */
    uint16_t record_count;    /* total number of records */
    uint64_t summary_offset;  /* on disk summary offset */
    uint32_t buf_size;        /* The compression buffer size (if the data is compressed) */
} BGTFHdr_t;

BGTFHdr_t *BGTFHdrInit();
void BGTFHdrDestroy(BGTFHdr_t *hdr);


/**
 * @brief The magic number is derived from `__ac_X31_hash_string` in khash.
 *        The same hash function is also used in HashMap storing attributes 
 *        from the GTF file
 * 
 */
enum BGTF_ATTR_MAGIC
{
    GENE_ID = 0xFB38C285,
    GENE_VERSION = 0xEEC739AE,
    GENE_NAME = 0x10147D75,
    TRANSCRIPT_ID = 0xF196D644,
    TRANSCRIPT_VERSION = 0xF52AA54F,
    TRANSCRIPT_NAME = 0xE73C9D74,
    EXON_NUMBER = 0x90DDCE56,
    EXON_ID = 0xB1EA8EA8,
    FAMILY_ID = 0xAEF4ABD6,
    CLASS_ID = 0x2945B542,
};

/**
 * @brief Additional information as key-value pairs for GTF filesca
 */
typedef struct RecordAttrKV {
    char     *gene_id;              /* can represent either an gene id (e.g. ensembl) or TE name */
    uint16_t  gene_version;         /* gene version */
    char     *gene_name;            /* name of the gene */
    char     *transcript_id;        /* can represent either an transcript id (e.g. ensembl) or TE duplicate */
    uint16_t  transcript_version;   /* transcript version */
    char     *transcript_name;      /* name of the transcript */
    uint16_t  exon_number;          /* name of the exon */
    char     *exon_id;              /* exon id (e.g. ensembl) */
    char     *family_id;            /* TE record from repeatmasker */
    char     *class_id;             /* TE record from repeatmasker */
} RecordAttrKV_t; 

RecordAttrKV_t *parseRecordAttrKV(char *src);
// TODO: Should move to internal API
uint8_t setRecordAttrKV_helper(RecordAttrKV_t *attr, uint32_t magic, char *v);
uint8_t setRecordAttrKV(RecordAttrKV_t *attr, char *k, char *v);

/**
 * @brief record for gene annotation
 * 
 */
typedef struct BGTFRecord BGTFRecord_t;
struct BGTFRecord {
    char *chrom;           /* chromosome name */
    char *source;          /* annotation source */
    char *type;            /* feature type */
    uint64_t start;        /* start of the genome coordinate */
    uint64_t end;          /* end of the genome coordinate */
    uint8_t  score;        /* this field is not used in GTF, however we keep it */
    uint8_t  strand;       /* strand of the record. 1 for `+` and 0 for `-` */
    uint32_t phase;        /* genomic phase (for CDS features) */
    RecordAttrKV_t *attrs; /* extra GTF attributes */
};

BGTFRecord_t *parseOneRecord(char *buf, char *line);
void BGTFRecordDestroy(BGTFRecord_t *record);
int8_t RecordCompare(BGTFRecord_t *record1, BGTFRecord_t *record2);

typedef struct BGTFRecordIdx
{
    BGTFRecord_t *head;
    BGTFRecord_t *tail;
    RTree *rt;
} BGTFRecordIdx;

typedef struct BGTFIdx
{
    HashTable *chrom_idx;
} BGTFIdx;

BGTFIdx *BGTFIdxInit();
void BGTFIdxDestroy(BGTFIdx *idx);



typedef struct {
	int errcode:16, is_write:2, is_be:2, compress_level:9, is_compressed:2, is_gzip:1;
	int cache_size;
    int block_length, block_offset;
    int64_t block_address, uncompressed_address;
    void *uncompressed_block, *compressed_block;
    z_stream *gz_stream; /* for gzip-compressed files */
} BGTFRWBuffer_t;

typedef struct BGTFRW
{
    FILE *fp;
    BGTFRWBuffer_t *buffer;
} BGTFRW;

typedef struct BGTF {
    ArrayList *records; // array list containing all records
    BGTFIdx   *index;   // the in-memory index for the BGTF file
    BGTFHdr_t *hdr;     // header of the BGTF file
    BGTFRW *write;   // Write handler
    uint8_t is_write;   // 
    size_t n_records;
} BGTF;

struct GenomicLocation
{
    char *chrom;
    uint64_t start;
    uint64_t end;
};

struct GenomicLocation *parseGenomicLocation(char *src);
void BGTFGetRecord(BGTF *fp, struct GenomicLocation *loc, void (*func)(BGTFRecord_t *record));
BGTF *BGTFInit();
void BGTFClose(BGTF *fp);
void BGTFSortRecord(BGTF *fp);
void BGTFflushRecord(FILE *__stream__, BGTF *fp, char *k);
uint8_t BGTFBuildIndex(BGTF *file);


#endif 