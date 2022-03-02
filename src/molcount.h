#ifndef BGTF_MCOUNT
#define BGTF_MCOUNT

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include "bgtf.h"
#include "list.h"
#include "bgtfIO.h"


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
    uint8_t splice_only;
    size_t n_col;
    size_t n_row;
    HashTable *col_index;
    HashTable *row_index;
    void *spliced;
    void *unspliced;
    size_t n_non_zero;
    uint8_t _non_zero_cache;
} ExpressionMatrix_t;

ExpressionMatrix_t *ExpressionMatrixInit(size_t n_col, size_t n_row, uint8_t sparse, uint8_t splice_only);
void ExpressionMatrixInitSetRowIndex(ExpressionMatrix_t *matrix, HashTable *row_index);
void ExpressionMatrixInitSetColIndex(ExpressionMatrix_t *matrix, HashTable *col_index);
void ExpressionMatrixSet(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced, uint8_t lock, uint16_t value);
void ExpressionMatrixInitInc(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced, uint8_t lock);
void ExpressionMatrixDestroy(ExpressionMatrix_t *matrix);
void ExpressionMatrixSetMt(ExpressionMatrix_t *matrix);
void ExpressionMatrixUnsetMt(ExpressionMatrix_t *matrix);

uint8_t molecularCount(char *bgtf_path, 
                    char *rmsk_path, 
                    char *attr_name,
                    char *rmsk_attr_name, 
                    char *bam_path,
                    char *bam_dir,
                    char *bam_tag,
                    uint8_t bam_paired,
                    char *barcode_path,
                    char *output_path,
                    char *bam_output_path,
                    int n_threads,
                    uint32_t n_alignment_per_thread);

#endif