#ifndef BGTF_MCOUNT
#define BGTF_MCOUNT

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include "bgtf.h"
#include "list.h"
#include "bgtfIO.h"

typedef struct BamFile {
    samFile *fp;
    bam_hdr_t *bam_hdr;
    bam1_t *aln;
} BamFile_t;

typedef struct BamAlignment {
    GenomicLocation *loc;
    uint8_t strand;
    uint8_t spliced_alignment;
    uint8_t spliced_annotation;
    LinkedList *spliced_loc;
    char *annotation;
    char *barcode;
} BamAlignment;


static inline BamFile_t *BamOpen(char *path)
{
    BamFile_t *bamfile = malloc(sizeof(BamFile_t));
    bamfile->fp = hts_open(path, "r");    // open bam file
    bamfile->bam_hdr = sam_hdr_read(bamfile->fp); // read header
    bamfile->aln = bam_init1();               // initialize an alignment
    return bamfile;
}


static inline void BamClose(BamFile_t *bamfile)
{
	bam_destroy1(bamfile->aln);
	sam_close(bamfile->fp);
    free(bamfile);
}

BamAlignment *BamAlignmentInit();
void BamAlignmentDestroy(BamAlignment *baln);

BamAlignment *parseOneAlignment(BamFile_t *bamfile);
uint8_t assignFeature(const uint64_t *rect, void *item, void *udata);
uint8_t annotateBam(BGTF *fp, BamAlignment *baln, char *attr_name);


typedef struct {
    char *feature_name;
    char *feature_type;
} AssignFeature_t;


typedef struct RecordLocs {
    LinkedList *locs;
    uint8_t strand;
    char *attr_name;
} RecordLocs;

typedef struct ExpressionMatrix {
    uint8_t n_threads;
    pthread_mutex_t mu;
    uint8_t is_sparse;
    size_t n_col;
    size_t n_row;
    HashTable *col_index;
    HashTable *row_index;
    void *spliced;
    void *unspliced;
    size_t n_non_zero;
} ExpressionMatrix_t;

ExpressionMatrix_t *ExpressionMatrixInit(size_t n_col, size_t n_row, uint8_t sparse);
void ExpressionMatrixInitSetRowIndex(ExpressionMatrix_t *matrix, HashTable *row_index);
void ExpressionMatrixInitSetColIndex(ExpressionMatrix_t *matrix, HashTable *col_index);
void ExpressionMatrixSet(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced, uint8_t trylock, uint16_t value);
void ExpressionMatrixInitInc(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced, uint8_t trylock);
void ExpressionMatrixDestroy(ExpressionMatrix_t *matrix);
void ExpressionMatrixSetMt(ExpressionMatrix_t *matrix);
void ExpressionMatrixUnsetMt(ExpressionMatrix_t *matrix);

uint8_t molecularCount(char *bgtf_path, 
                    char *attr_name, 
                    char *bam_path,
                    char *bam_tag,
                    char *barcode_path,
                    char *output_path,
                    char *bam_output_path,
                    int n_threads);

#endif