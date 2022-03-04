#include <getopt.h>
#include <htslib/bgzf.h>

#include "tpool.h"
#include "molcount.h"
#include "sparse.h"
#include "bgtf-bam.h"

uint8_t assignFeature(const uint64_t *rect, void *item, void *udata)
{
    BGTFRecord_t *record = (BGTFRecord_t *)item;
   
        GenomicLocation *loc = GenomicLocationNew();
        loc->chrom = record->chrom;
        loc->start = record->start;
        loc->end = record->end;
        loc->data = malloc(sizeof(AssignFeature_t));
        ((AssignFeature_t *)loc->data)->feature_name = getRecordAttrKV(record->attrs, ((RecordLocs *)udata)->attr_name);
        ((AssignFeature_t *)loc->data)->feature_type = record->type;
        LinkedListAppend(GenomicLocation, ((RecordLocs *)udata)->locs, loc);
    
    return 1;
}


void AnnotateLocDestroyHelper(void *e)
{
    destroyGenomicLocation(e, free);
}


RecordLocs *RecordLocsInit(uint8_t strand, char *attr_name)
{
    RecordLocs *record_locs = malloc(sizeof(RecordLocs));
    memset(record_locs, 0, sizeof(RecordLocs));
    record_locs->locs = LinkedListInit();;
    record_locs->strand = strand;
    record_locs->attr_name = attr_name;
    LinkedListSetDealloc(record_locs->locs, AnnotateLocDestroyHelper);
    return record_locs;
}


void RecordLocsDestroy(RecordLocs *record_locs)
{
    LinkedListDestroy(GenomicLocation, record_locs->locs, 1);
    free(record_locs);
}

uint8_t annotateBam(BGTF *fp, BamAlignment *baln, char *attr_name)
{

    RecordLocs *record_locs = RecordLocsInit(baln->strand, attr_name);
    GenomicLocation *record_loc, *loc;
    if (BGTFGetRecord(fp, baln->loc, assignFeature, record_locs) > 0 || (record_locs->locs->sz == 0))
    {
        RecordLocsDestroy(record_locs);
        return 1;
    };

    HashTable *counter = HashTableCreate(0x4);
    HashTableSetKeyComparisonFunction(counter, strncmp);
    HashTableSetHashFunction(counter, HashTableStringHashFunction);
    
    uint8_t success_flag = 0;
    uint8_t spliced_flag = 0;
    char *feature_name;
    if (baln->spliced_alignment > 0) {
        /* The CIGAR suggests that this read is spliced. We need to tell 
           Whether this is an perfectly spliced read or unspliced read */
        baln->spliced_annotation = 1;
        GenomicLocation *first = baln->spliced_loc->elem_head;
        GenomicLocation *last = baln->spliced_loc->elem_tail;
        LinkedListMergeSort(GenomicLocation, record_locs->locs, GenomicLocationCompare);

        // The record locations are already sorted
        loc = NULL;
        LinkedListForEach(GenomicLocation, record_locs->locs, loc, {
            if ((loc->start <= first->start) && ((loc->end == first->end) || (loc->end > first->end))) {
                feature_name = ((AssignFeature_t *)loc->data)->feature_name;
                counterAdd(counter, feature_name);
                spliced_flag = 1;
            }
            if ((loc->start >= first->start) && (!spliced_flag)) {
                // Early stopping
                baln->spliced_annotation = 0;
                break;
            }
        });
        // We cannot find a perfect or approximate spliced junction of this alignment
        if (!spliced_flag) baln->spliced_annotation = 0;
        if (baln->spliced_alignment > 1) {
            // More than one splice junctions
            spliced_flag = 0;
            GenomicLocation *curr = first->next;
            while (curr != last) {
                LinkedListForEach(GenomicLocation, record_locs->locs, loc, {
                    if (((loc->start == curr->start) && (loc->end >= curr->end)) || 
                        ((loc->end == curr->end) && (loc->start <= curr->start))
                    ) {
                        // We find one or more accurate splice junctions.
                        spliced_flag = 1;
                        feature_name = ((AssignFeature_t *)loc->data)->feature_name;
                    } 
                    if ((loc->start >= curr->start) && (!spliced_flag)) {
                        // We cannot find a perfect spliced junction of this alignment
                        baln->spliced_annotation = 0;
                    }
                })                
                curr = curr->next;
            }
        }
        spliced_flag = 0;
        LinkedListForEach(GenomicLocation, record_locs->locs, loc, {
            if ((loc->start == last->start) && (loc->end >= last->end)) {
                feature_name = ((AssignFeature_t *)loc->data)->feature_name;
                counterAdd(counter, feature_name);
                spliced_flag = 1;
            }
            if ((loc->start > last->start) && (!spliced_flag)) {
                // early stopping for unspliced
                baln->spliced_annotation = 0;
                break;
            }
        });
        // We cannot find a perfect spliced junction of this alignment
        if (!spliced_flag) baln->spliced_annotation = 0; 
       
    } else {
        /* The read is not aligned to be spliced. We need to tell whether a record can
           tell that this read is perfectly spliced. Otherwise it will be annotated as 
           as an unspliced read. */
        uint8_t spliced = 0;
        LinkedListForEach(GenomicLocation, record_locs->locs, loc, {
            if ((strcmp( (((AssignFeature_t *)loc->data)->feature_type), "gene") != 0) &&
                 (strcmp( (((AssignFeature_t *)loc->data)->feature_type), "transcript") != 0)
            ) {
                if (loc->start <= baln->loc->start && loc->end >= baln->loc->end) {
                    spliced = 1;
                }
            }
            feature_name = ((AssignFeature_t *)loc->data)->feature_name;
            counterAdd(counter, feature_name);
        });
        if (spliced) baln->spliced_annotation = 1;
    }

    char *key; int value; char *max_key; int max_value = 0;
    if (HashTableIsEmpty(counter)) {
        success_flag = 1;
        goto cleanup;
    }
    HashTableForEach(counter, key, value, {
        if (value > max_value) {
            max_key = key; 
            max_value = value; // max record supporting the assignment
        }
    })
    baln->annotation = max_key;
    goto cleanup;
cleanup:
    RecordLocsDestroy(record_locs);
    HashTableDestroy(counter);
    return success_flag;
}

static uint8_t annotateBamRmsk(BGTF *gfp, BGTF* tfp, BamAlignment *baln, char *gene_attr_name, char *rmsk_attr_name)
{
    uint8_t ret = annotateBam(gfp, baln, gene_attr_name);
    if ((ret == 0) && (baln->annotation)) return ret;
    if ((!baln->spliced_alignment) && (!baln->spliced_annotation)) {
        // We cannot find an Gene record to support this read
        // Also, this read is not aligned to be spliced
        RecordLocs *record_locs = RecordLocsInit(baln->strand, rmsk_attr_name);
        GenomicLocation *record_loc, *loc;
        if (BGTFGetRecord(tfp, baln->loc, assignFeature, record_locs) > 0 || (record_locs->locs->sz == 0))
        {
            RecordLocsDestroy(record_locs);
            return 1;
        };
        uint8_t success_flag = 0;
        HashTable *counter = HashTableCreate(0x4);
        HashTableSetKeyComparisonFunction(counter, strncmp);
        HashTableSetHashFunction(counter, HashTableStringHashFunction);
        char *feature_name;
        LinkedListForEach(GenomicLocation, record_locs->locs, loc, {
            if (loc->start <= baln->loc->start && loc->end >= baln->loc->end) {
                feature_name = ((AssignFeature_t *)loc->data)->feature_name;
                counterAdd(counter, feature_name);  
            }

        });
        char *key; int value; char *max_key; int max_value = 0;
        if (HashTableIsEmpty(counter)) {
            success_flag = 1;
            RecordLocsDestroy(record_locs);
            HashTableDestroy(counter);
            return success_flag & ret;
        }
        HashTableForEach(counter, key, value, {
            if (value > max_value) {
                max_key = key; 
                max_value = value; // max record supporting the assignment
            }
        })
        baln->annotation = max_key;
        baln->is_rmsk = 1;
        RecordLocsDestroy(record_locs);
        HashTableDestroy(counter);
        return success_flag;
    } else {
        return 1;
    }
}

uint8_t annotateBamPaired(BGTF *fp, PairedAlignment *paln, char *attr_name)
{
    // TODO: Rewrite this function
    if (!paln->paired) {
        // The read is not properly paired
        return annotateBam(fp, paln->mate1, attr_name);
    } else {
        uint8_t success1 = 0, success2 = 0;
        success1 = annotateBam(fp, paln->mate1, attr_name);
        success2 = annotateBam(fp, paln->mate2, attr_name);

        if (success1 == 0) {
            paln->annotation = paln->mate1->annotation;
        } else if (success2 == 0) {
            paln->annotation = paln->mate2->annotation;
        }
        paln->spliced_annotation = paln->mate1->spliced_annotation & paln->mate2->spliced_annotation;
        return !((!success1) | (!success2));
    }
}

uint8_t annotateBamPairedRmsk(BGTF *gfp, BGTF* tfp, PairedAlignment *paln, char *gene_attr_name, char *rmsk_attr_name)
{
    // TODO: Rewrite this function
    if (!paln->paired) {
        // The read is not properly paired
        return annotateBamRmsk(gfp, tfp, paln->mate1, gene_attr_name, rmsk_attr_name);
    } else {
        uint8_t ret = annotateBamPaired(gfp, paln, gene_attr_name);
        if ((ret == 0) && (paln->annotation)) return ret;
        if ((!paln->spliced_annotation) && (!paln->mate1->spliced_alignment) && (!paln->mate2->spliced_alignment))
        {
            annotateBamRmsk(gfp, tfp, paln->mate1, gene_attr_name, rmsk_attr_name);
            annotateBamRmsk(gfp, tfp, paln->mate2, gene_attr_name, rmsk_attr_name);
            if (strcmp(paln->mate1->annotation, paln->mate2->annotation) == 0) {
                paln->annotation = paln->mate1->annotation;
                paln->is_rmsk = 1;
                return 0;
            }
        } else {
            return 1;
        }
    }
}

static ExpressionMatrix_t *ExpressionMatrixDenseInit(size_t n_row, size_t n_col, uint8_t splice_only)
{
    logf("Initialize dense matrix %d × %d", n_row, n_col);
    ExpressionMatrix_t *matrix = malloc(sizeof(ExpressionMatrix_t));
    memset(matrix, 0, sizeof(ExpressionMatrix_t));
    matrix->spliced = calloc(n_col * n_row, sizeof(uint16_t));
    memset(matrix->spliced, 0, n_col * n_row * sizeof(uint16_t));
    if (!splice_only) {
        matrix->unspliced = calloc(n_col * n_row, sizeof(uint16_t));
        memset(matrix->unspliced, 0, n_col * n_row * sizeof(uint16_t));
    }
    return matrix;
}

static ExpressionMatrix_t *ExpressionMatrixSparseInit(size_t n_row, size_t n_col, uint8_t splice_only)
{
    logf("Initialize sparse matrix %d × %d", n_row, n_col);
    ExpressionMatrix_t *matrix = malloc(sizeof(ExpressionMatrix_t));
    memset(matrix, 0, sizeof(ExpressionMatrix_t));
    matrix->spliced = newMatrix(n_row, n_col);
    if (!splice_only) matrix->unspliced = newMatrix(n_row, n_col);
    return matrix;
}

ExpressionMatrix_t *ExpressionMatrixInit(size_t n_row, size_t n_col, uint8_t sparse, uint8_t splice_only)
{
    ExpressionMatrix_t *matrix;
    if (sparse) {
        matrix = ExpressionMatrixSparseInit(n_row, n_col, splice_only);
    } else {
        matrix = ExpressionMatrixDenseInit(n_row, n_col, splice_only);
    }
    matrix->is_sparse = sparse;
    matrix->n_row = n_row;
    matrix->n_col = n_col;
    matrix->splice_only = splice_only;    
    matrix->n_non_zero = 0;
    return matrix;
}

void ExpressionMatrixInitSetRowIndex(ExpressionMatrix_t *matrix, HashTable *row_index)
{
    matrix->row_index = row_index;
}

void ExpressionMatrixInitSetColIndex(ExpressionMatrix_t *matrix, HashTable *col_index)
{
    matrix->col_index = col_index;
}

uint16_t ExpressionMatrixInitGet(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced) {
    int64_t col = HashTableGet(matrix->col_index, col_name) - 1;
    int64_t row = HashTableGet(matrix->row_index, row_name) - 1; 
    if (col < 0 || row < 0) return 0; 
    if (matrix->is_sparse) {
        Matrix exp;
        if (spliced) {
            exp = (Matrix) matrix->spliced;
        } else {
            exp = (Matrix) matrix->unspliced;
        }
        return getEntryMatrix(exp, row, col);
    } else {
        uint16_t *exp;
        if (spliced) {
            exp = (uint16_t *) matrix->spliced;
        } else {
            exp = (uint16_t *) matrix->unspliced;
        }
        return exp[row * matrix->n_col + col];
    }
}

void ExpressionMatrixSet(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced, uint8_t lock, uint16_t value)
{
    matrix->_non_zero_cache = 0;
    int64_t col = HashTableGet(matrix->col_index, col_name) - 1;
    int64_t row = HashTableGet(matrix->row_index, row_name) - 1; 
    if (col < 0 || row < 0) return;
    if (matrix->is_sparse) {
        // TODO: Sparse implementation
        Matrix exp;
        if (spliced) {
            exp = (Matrix) matrix->spliced;
        } else {
            exp = (Matrix) matrix->unspliced;
        }
        if (matrix->n_threads && lock) {
            while (pthread_mutex_lock(&matrix->mu) != 0);
        }
        addEntryToIMatrix(exp, row, col, value);
        if (matrix->n_threads && lock) {
            pthread_mutex_unlock(&matrix->mu);
        }
    } else {
        uint16_t *spliced_exp, *unspliced_exp;;
        spliced_exp = (uint16_t *) matrix->spliced;
        if (matrix->n_threads && lock) {
            while (pthread_mutex_lock(&matrix->mu) != 0);
        }
        if (!matrix->splice_only)
        {
            unspliced_exp = (uint16_t *) matrix->unspliced;
            if (spliced) spliced_exp[row * matrix->n_col + col] = value;
            else unspliced_exp[row * matrix->n_col + col] = value;
        } else {
            spliced_exp[row * matrix->n_col + col] = value;
        }
        if (matrix->n_threads) {
            pthread_mutex_unlock(&matrix->mu);
        }
    }
}



void ExpressionMatrixInc(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced, uint8_t lock)
{
    matrix->_non_zero_cache = 0;
    if (row_name == NULL) {
        warnf("The row_name is %s for this read", row_name);
        return;
    }
    if (col_name == NULL) {
        warnf("The col_name is %s for this read", col_name);
        return;
    }
    int64_t col = HashTableGet(matrix->col_index, col_name) - 1;
    int64_t row = HashTableGet(matrix->row_index, row_name) - 1;
    if (col < 0 || row < 0) return;
    if (matrix->is_sparse) {
        // TODO: Sparse implementation
        Matrix exp;
        if (spliced) {
            exp = (Matrix) matrix->spliced;
        } else {
            exp = (Matrix) matrix->unspliced;
        }
        if (matrix->n_threads && lock) {
            while (pthread_mutex_lock(&matrix->mu) != 0);
        }
        accumulateEntryInIMatrix(exp, row, col, 1);
        if (matrix->n_threads && lock) {
            pthread_mutex_unlock(&matrix->mu);
        }
    } else {
        uint16_t *spliced_exp, *unspliced_exp;;
        spliced_exp = (uint16_t *) matrix->spliced;
        
        if (matrix->n_threads && lock) {
            while (pthread_mutex_lock(&matrix->mu) != 0);
        }
        
        if (!matrix->splice_only) {
            unspliced_exp = (uint16_t *) matrix->unspliced;
            if (spliced) {
                spliced_exp[row * matrix->n_col + col]++;
            } else {
                unspliced_exp[row * matrix->n_col + col]++;
            }
            
        } else {
            spliced_exp[row * matrix->n_col + col]++;
        }
        if (matrix->n_threads && lock) {
            pthread_mutex_unlock(&matrix->mu);
        }
    }
}

void ExpressionMatrixDestroy(ExpressionMatrix_t *matrix)
{   
    if (matrix->is_sparse) {

    } else {
        if (matrix->spliced) free(matrix->spliced);
        if (matrix->unspliced) free(matrix->unspliced);
        free(matrix);
    }
}
 
static size_t ExpressionMatrixnNonZero(ExpressionMatrix_t *matrix)
{
    if (matrix->_non_zero_cache) return matrix->n_non_zero;
    size_t nonzeros = 0;
    size_t i,j;
    uint16_t v,u=0;
    if (matrix->is_sparse) {
    } else {
        uint16_t *spliced_exp, *unspliced_exp;
        spliced_exp = (uint16_t *) matrix->spliced;
        if (!matrix->splice_only) unspliced_exp = (uint16_t *) matrix->unspliced;
        for (i = 0; i < matrix->n_row; ++i) {
            for (j = 0; j < matrix->n_col; ++j) {
                v = spliced_exp[i * matrix->n_col + j];
                if (!matrix->splice_only) u = unspliced_exp[i * matrix->n_col + j];
                if (v > 0 || u > 0) {
                    nonzeros++;
                }
            }
        }
    }
    matrix->n_non_zero = nonzeros;
    matrix->_non_zero_cache = 1;
    return nonzeros;
}

uint8_t ExpressionMatrixToFile(ExpressionMatrix_t *matrix, FILE *spliced_stream, FILE *unspliced_stream)
{
    ExpressionMatrixnNonZero(matrix);
    fprintf(spliced_stream, "%s", "\%\%MatrixMarket matrix coordinate integer general\n");
    fprintf(spliced_stream, "%s", "\%\n");
    fprintf(spliced_stream, "%llu %llu %llu\n", matrix->n_col, matrix->n_row, matrix->n_non_zero);
    fprintf(unspliced_stream, "%s", "\%\%MatrixMarket matrix coordinate integer general\n");
    fprintf(unspliced_stream, "%s", "\%\n");
    fprintf(unspliced_stream, "%llu %llu %llu\n", matrix->n_col, matrix->n_row, matrix->n_non_zero);
    size_t i,j;
    uint16_t v,u;
    if (matrix->is_sparse) {
    } else {
        uint16_t *spliced_exp, *unspliced_exp;
        spliced_exp = (uint16_t *) matrix->spliced;
        unspliced_exp = (uint16_t *) matrix->unspliced;
        for (i = 0; i < matrix->n_row; ++i) {
            for (j = 0; j < matrix->n_col; ++j) {
                v = spliced_exp[i * matrix->n_col + j];
                u = unspliced_exp[i * matrix->n_col + j];
                if (v > 0 || u > 0) {
                    fprintf(spliced_stream, "%llu %llu %d\n", j+1, i+1, v);
                    fprintf(unspliced_stream, "%llu %llu %d\n", j+1, i+1, u);
                }
            }
        }
    }
    return 0;
}

uint8_t ExpressionMatrixToFile_(ExpressionMatrix_t *matrix, FILE *spliced_stream)
{
    ExpressionMatrixnNonZero(matrix);
    fprintf(spliced_stream, "%s", "\%\%MatrixMarket matrix coordinate integer general\n");
    fprintf(spliced_stream, "%s", "\%\n");
    fprintf(spliced_stream, "%llu %llu %llu\n", matrix->n_col, matrix->n_row, matrix->n_non_zero);
    uint16_t *spliced_exp = (uint16_t *)matrix->spliced;
    size_t i,j;
    uint16_t v;
    for (i = 0; i < matrix->n_row; ++i)
    {
        for (j = 0; j < matrix->n_col; ++j)
        {
            v = spliced_exp[i * matrix->n_col + j];
            if (v > 0)
            {
                fprintf(spliced_stream, "%llu %llu %d\n", j + 1, i + 1, v);
            }
        }
    }
}

void ExpressionMatrixSetMt(ExpressionMatrix_t *matrix)
{
    pthread_mutexattr_t attr;
    pthread_mutexattr_init(&attr);
    pthread_mutex_init(&matrix->mu, &attr);
    matrix->n_threads = 1;
}

void ExpressionMatrixUnsetMt(ExpressionMatrix_t *matrix)
{
    pthread_mutex_destroy(&matrix->mu);
    matrix->n_threads = 0;
}

typedef struct {
    void **aln;                  /* Either BamAlignment or PairedAlignment */
    size_t n_bam_alignments;      /* Number of BamAlignments or PairedAlignments */
    BGTF *fp;                     /* GTF file */
    BGTF *rmsk;                   /* repeat masker GTF file */
    ExpressionMatrix_t *mtx;      /* Output matrix */
    ExpressionMatrix_t *rmsk_mtx; /* Repeat masker output matrix */
    char *attr_name;              /* Attribute name to be count */
    char *rmsk_attr_name;         /* Repeat masker attribute name to be count */
    int64_t rid;                  /* ArgID. Not useful now  */
    uint8_t bam_paired;           /* indicate whether member baln is BamAlignments or PairedAlignments */
    void *free_ptr;               /* Pointer to destroy the ArrayList wrapper */
} ParallelizerArg;

typedef struct {
    char *barcode;
    char *annotation;
    uint8_t spliced_annotation;
    uint8_t is_rmsk;
} ParallelizerResult;

void ParallelizerArgClean(void *arg)
{
    ParallelizerArg *pa = (ParallelizerArg *)arg;
    if (pa->bam_paired) {
        for (size_t i = 0; i < pa->n_bam_alignments; ++i) {
            PairedAlignmentDestroy(pa->aln[i]);
        }
    } else {
        for (size_t i = 0; i < pa->n_bam_alignments; ++i) {
            BamAlignmentDestroy(pa->aln[i]);
        }
    }
    ArrayListDestroy(pa->free_ptr);
    free(pa);
}

void molecularCountParallelizer(void *arg) {
    ParallelizerArg *pa = (ParallelizerArg *)arg;
    ParallelizerResult *result;
    ArrayList *par_result = ArrayListCreate(pa->n_bam_alignments);
    ArrayListSetDeallocationFunction(par_result, free);
    uint8_t ret;
    if (pa->bam_paired) {
        for (size_t i = 0; i < pa->n_bam_alignments; ++i) {
            if (pa->rmsk_mtx) {
                ret = annotateBamPairedRmsk(pa->fp,
                                    pa->rmsk, 
                                    pa->aln[i],
                                    pa->attr_name,
                                    pa->rmsk_attr_name);
            } else {
                ret = annotateBamPaired(pa->fp, 
                                pa->aln[i], 
                                pa->attr_name);
            }
            if (ret == 0) {
                char *annotation = ((PairedAlignment *) pa->aln[i])->annotation;
                result = malloc(sizeof(ParallelizerResult));
                result->barcode = ((PairedAlignment *) pa->aln[i])->barcode;
                result->annotation = annotation;
                result->spliced_annotation = ((PairedAlignment *) pa->aln[i])->spliced_annotation;
                result->is_rmsk = ((PairedAlignment *) pa->aln[i])->is_rmsk;
                ArrayListPush(par_result, result);
            }
        }
    } else {
        for (size_t i = 0; i < pa->n_bam_alignments; ++i) {
            
            if (pa->rmsk != NULL) {
                ret = annotateBamRmsk(pa->fp, 
                                    pa->rmsk, 
                                    pa->aln[i],
                                    pa->attr_name,
                                    pa->rmsk_attr_name);
            } else {
                ret = annotateBam(pa->fp, 
                                pa->aln[i], 
                                pa->attr_name);
            }
            if (ret == 0) {
                char *annotation = ((BamAlignment *) pa->aln[i])->annotation;
                result = malloc(sizeof(ParallelizerResult));
                result->barcode = ((BamAlignment *) pa->aln[i])->barcode;
                result->annotation = annotation;
                result->spliced_annotation = ((BamAlignment *) pa->aln[i])->spliced_annotation;
                result->is_rmsk = ((BamAlignment *) pa->aln[i])->is_rmsk;
                ArrayListPush(par_result, result);
            }
        }
    }
    pthread_mutex_lock(&pa->mtx->mu);
    if (pa->rmsk_mtx) pthread_mutex_lock(&pa->rmsk_mtx->mu);

    ArrayListForEach(par_result, result, {
        if (result->is_rmsk) {
            ExpressionMatrixInc(pa->rmsk_mtx, result->barcode, result->annotation, 1, 0);
        } else {
            ExpressionMatrixInc(pa->mtx, result->barcode, result->annotation,   result->spliced_annotation, 0);
        }
    })

    if (pa->rmsk_mtx) pthread_mutex_unlock(&pa->rmsk_mtx->mu);
    pthread_mutex_unlock(&pa->mtx->mu);
    
    ArrayListDestroy(par_result);
    ParallelizerArgClean(arg);
}

uint8_t molecularCountInternal(
                    BGTF *rfp,
                    BGTF *tfp,
                    BGZF *ofp,
                    char *attr_name,
                    char *rmsk_attr_name, 
                    ExpressionMatrix_t *matrix,
                    ExpressionMatrix_t *rmsk_matrix,
                    char *bam_tag,
                    uint8_t bam_paired,
                    char *bam_barcode,
                    int n_threads, 
                    tpool_t *p,
                    tpool_process_t *q,
                    uint32_t n_alignment_per_thread,
                    BamFile_t *bamfile
                );

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
                    uint32_t n_alignment_per_thread)
{
    log("Using arguments:")

    logf("BAM CellBarcode tag: %s", bam_tag);
    logf("Using (%d/%d) threads", n_threads, MAX_THREADS);

    BGTF *rfp = NULL;
    BGTF *tfp = NULL;
    BGZF *ofp = NULL;
    
    if (bgtf_str_endswith(bgtf_path, "bgtf")) {
        rfp = BGTFInit();
        BGTFLoad(rfp, bgtf_path, attr_name);
    } else if (bgtf_str_endswith(bgtf_path, "gtf")) {
        rfp = GTFRead(bgtf_path, attr_name);
    } else {
        fatalf("The input reference file %s has wrong extension. \n"
               "Please check whether the file ends with .gtf or .bgtf", bgtf_path);
    }

    
    if (rmsk_path) {
        if (bgtf_str_endswith(rmsk_path, "bgtf")) {
            tfp = BGTFInit();
            BGTFLoad(tfp, rmsk_path, rmsk_attr_name);
        } else if (bgtf_str_endswith(rmsk_path, "gtf")) {
            tfp = GTFRead(rmsk_path, rmsk_attr_name);
        } else {
            fatalf("The input reference file %s has wrong extension. \n "
            "Please check whether the file ends with .gtf or .bgtf", rmsk_path);
        }
    }

    BamFile_t *bamfile;

    if (bam_path) {
        bamfile = BamOpen(bam_path);
        if (bam_output_path) {
            ofp = bgzf_open(bam_output_path, "wb+");
            bam_hdr_write(ofp, bamfile->bam_hdr);
        }
    }

    HashTable *barcodeIndex;

    uint8_t free_bam_dir = 0;
    if (bam_dir && (bam_dir[strlen(bam_dir)-1]) != '/') {
        char *bam_dir_ = bam_dir;
        bam_dir =  malloc(strlen(bam_dir_) + 2);
        memset(bam_dir, 0, strlen(bam_dir_) + 2);
        memcpy(bam_dir, bam_dir_, strlen(bam_dir_));
        bam_dir[strlen(bam_dir_)] = '/';
        free_bam_dir = 1;
    }

    if (barcode_path) {
        barcodeIndex = readBarcodeFile(barcode_path);
    } else {
        barcodeIndex = bgtf_listdir(bam_dir, NULL, "bam");
    }

    ExpressionMatrix_t *matrix = ExpressionMatrixInit(HashTableSize(barcodeIndex), HashTableSize(rfp->attr_table), 0, 0);
    ExpressionMatrix_t *rmsk_matrix = NULL;
    ExpressionMatrixInitSetColIndex(matrix, rfp->attr_table);
    ExpressionMatrixInitSetRowIndex(matrix, barcodeIndex);

    if (rmsk_path) {
        rmsk_matrix = ExpressionMatrixInit(HashTableSize(barcodeIndex), HashTableSize(tfp->attr_table), 0, 1);
        ExpressionMatrixInitSetColIndex(rmsk_matrix, tfp->attr_table);
        ExpressionMatrixInitSetRowIndex(rmsk_matrix, barcodeIndex);
        if (n_threads > 1) {
            ExpressionMatrixSetMt(rmsk_matrix);
        }
    }


    
    char *k;
    size_t attr_idx;
    ArrayList *sort_arr;
    sort_arr = ArrayListCreate(HashTableSize(rfp->attr_table));
    HashTableForEach(rfp->attr_table, k, attr_idx, {
        ArrayListPush(sort_arr, k);
    })
    
    uint8_t free_output_path = 0;
    if (output_path && (output_path[strlen(output_path)-1]) != '/') {
        char *output_path_ = output_path;
        output_path =  malloc(strlen(output_path_) + 2);
        memset(output_path, 0, strlen(output_path_) + 2);
        memcpy(output_path, output_path_, strlen(output_path_));
        output_path[strlen(output_path_)] = '/';
        free_output_path = 1;
    }

    char *gene_name_output_path = bgtf_str_concat(output_path, "features.tsv");
    FILE *gene_name_out;
    if (!(gene_name_out = fopen(gene_name_output_path, "w+"))) {
        fatalf("failed to open %s", gene_name_output_path);
    }
    logf("Output features name to %s", gene_name_output_path);
    ArrayListSort(sort_arr, strcmp);
    ArrayListForEach(sort_arr, k, {
        fprintf(gene_name_out, "%s\n", k);
    })
    ArrayListDestroy(sort_arr);
    fclose(gene_name_out);
    free(gene_name_output_path);
    
    if (rmsk_path) {
        char *rmsk_name_output_path;
        rmsk_name_output_path = bgtf_str_concat(output_path, "rmsk_features.tsv");        
        logf("Output features name to %s", rmsk_name_output_path);
        FILE *rmsk_name_out;
        if (!(rmsk_name_out = fopen(rmsk_name_output_path, "w+"))) {
            fatalf("failed to open %s", rmsk_name_output_path);
        }
        sort_arr = ArrayListCreate(HashTableSize(tfp->attr_table));
        HashTableForEach(tfp->attr_table, k, attr_idx, {
            ArrayListPush(sort_arr, k);
        })
        ArrayListSort(sort_arr, strcmp);
        ArrayListForEach(sort_arr, k, {
            fprintf(rmsk_name_out, "%s\n", k);
        })
        ArrayListDestroy(sort_arr);
        fclose(rmsk_name_out);
        free(rmsk_name_output_path);
    }
    sort_arr = ArrayListCreate(HashTableSize(barcodeIndex));
    HashTableForEach(barcodeIndex, k, attr_idx, {
        ArrayListPush(sort_arr, k);
    })
    char *barcode_output_path = bgtf_str_concat(output_path, "barcodes.tsv");
    FILE *barcode_out;
    if (!(barcode_out = fopen(barcode_output_path, "w+"))) {
        fatalf("Failed to open %s", barcode_output_path);
    }

    logf("Output barcodes to %s", barcode_output_path);
    ArrayListSort(sort_arr, strcmp);
    ArrayListForEach(sort_arr, k, {
        fprintf(barcode_out, "%s\n", k);
    })
    fclose(barcode_out);
    free(barcode_output_path);
    ArrayListDestroy(sort_arr);

    tpool_t *p = NULL;
    tpool_process_t *q = NULL;
    if (n_threads > 1) {
        pthread_setconcurrency(2);
        p = tpool_init(n_threads);
        q = tpool_process_init(p, 16, true);
        ExpressionMatrixSetMt(matrix);
    }

    log("Start counting molecules");
    if (bam_path) {
        molecularCountInternal(
            rfp, 
            tfp, 
            ofp,
            attr_name,
            rmsk_attr_name,
            matrix,
            rmsk_matrix,
            bam_tag,
            bam_paired,
            NULL,
            n_threads, p, q,
            n_alignment_per_thread,
            bamfile);
    } else {
        size_t idx;
        HashTableForEach(barcodeIndex, bam_path, idx, {
            logf("Processing %s", bam_path);
            bamfile = BamOpen(bam_path);
                molecularCountInternal(
                    rfp, 
                    tfp, 
                    ofp,
                    attr_name,
                    rmsk_attr_name,
                    matrix,
                    rmsk_matrix,
                    bam_tag,
                    bam_paired,
                    bam_path, // cell barcode provided here
                    n_threads, p, q,
                    n_alignment_per_thread,
                    bamfile
                );
            BamClose(bamfile);
        });
    }

    if (n_threads > 1) {
        ExpressionMatrixUnsetMt(matrix);
        tpool_process_destroy(q);
        tpool_destroy(p);
    }

    log("Counting finished");
    char *mtx_output_path = calloc(strlen(output_path) + 12, 1);
    memcpy(mtx_output_path, output_path, strlen(output_path));
    memcpy(mtx_output_path + strlen(output_path), "spliced.mtx", 11);
    FILE *spliced;
    if (!(spliced = fopen(mtx_output_path, "w+"))) {
        fatalf("Failed to open %s", mtx_output_path);
    };
    free(mtx_output_path);
    mtx_output_path = calloc(strlen(output_path) + 14, 1);
    memcpy(mtx_output_path, output_path, strlen(output_path));
    memcpy(mtx_output_path + strlen(output_path), "unspliced.mtx", 13);
    FILE *unspliced;
    if (!(unspliced = fopen(mtx_output_path, "w+"))) {
        fatalf("Failed to open %s", mtx_output_path);
    };
    free(mtx_output_path);
    mtx_output_path = calloc(strlen(output_path) + 9, 1);
    memcpy(mtx_output_path, output_path, strlen(output_path));
    memcpy(mtx_output_path + strlen(output_path), "rmsk.mtx", 8);
    FILE *rmsk;
    if (rmsk_path) {
        if (!(rmsk = fopen(mtx_output_path, "w+"))) {
            fatalf("Failed to open %s", mtx_output_path);
        };
    }
    free(mtx_output_path);
    if (free_output_path) free(output_path);
    if (free_bam_dir) free(bam_dir);
    log("Writing count matrix");
    ExpressionMatrixToFile(matrix, spliced, unspliced);
    if (rmsk_path) {
        ExpressionMatrixToFile_(rmsk_matrix, rmsk);
    }
    fclose(spliced);
    fclose(unspliced);
    if (rmsk_path) fclose(rmsk);
    ExpressionMatrixDestroy(matrix);
    if (rmsk_path) {
        ExpressionMatrixDestroy(rmsk_matrix);
        if (n_threads > 1) ExpressionMatrixUnsetMt(rmsk_matrix);
    }
    destroyBarcodeTable(barcodeIndex);
    if (!bam_dir) BamClose(bamfile);
    if (bam_output_path) {
        bgzf_close(ofp);
    }
    return 0;
}


/* Command-line function */
void count_usage()
{
    fprintf(stdout,    
"Usage:   bgtf count [options]\n"
"\n"
"Arguments:\n"
"   Required arguments:\n"
"     -c          [PATH]   Cell barcode file \n"
"     -i          [STRING] BAM tag of barcode. Default is CB \n"
"     -g          [STRING] feature name to be counted \n"
"     -r          [PATH]   Gene reference file with extension .gtf or .bgtf \n"
"   Optional arguments:\n"
"     -b          [PATH]   Input BAM file \n"
"     -B          [PATH]   output gene annotation to a BAM file if set \n"
"     -d          [PATH]   Input BAM directory \n"
"     -m          [PATH]   RMSK GTF file. Count TE if set.  with extension .gtf or .bgtf \n"
"     -t          [INT]    n parallel thread if setted \n"
"     -T          [INT]    n alignment to process in one thread. Default is 1000000 \n"
"\n");
    return;
}
void count_version()
{
    fprintf(stdout,    
"bgtf count version %s\n", BGTF_COUNT_VERSION
    );
    return;
}

int count_main(int argc, char *argv[])
{
    int c;
    char *file_path = NULL;
    char *rmsk_path = NULL;
    char *output_path = NULL;
    char *bam_path = NULL;
    char *bam_dir = NULL;
    char *bam_output_path = NULL;
    char *bam_tag = NULL;
    char *barcode_path = NULL;
    char *attr;
    int ret;
    uint8_t paired = 0;
    int n_threads = 1;
    uint32_t n_alignment_per_thread = 1000000;
    while (1)
    {
        static struct option long_options[] =
            {   
                {"bamOutput", optional_argument, 0, 'B'},
                {"bamFile", optional_argument, 0, 'b'},
                {"cellBarcode", optional_argument, 0, 'c'},
                {"bamDir", optional_argument, 0, 'd'},
                {"attrName", required_argument, 0, 'g'},
                {"bamTag", required_argument, 0, 'i'},
                {"rmsk", optional_argument, 0, 'm'},
                {"gtf", required_argument, 0, 'r'},
                {"output", required_argument, 0, 'o'},
                {"nthreads", optional_argument, 0, 't'},
                {"nalignment", optional_argument, 0, 'T'},
                {"paired", no_argument, 0, 'p'},
                {"help", no_argument, NULL, 'h'},
                {"version", no_argument, NULL, 'v'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        c = getopt_long(argc, 
            argv, 
            "B:b:c:d:g:hi:m:o:pt:T:r:v", 
            long_options, 
            &option_index
        );

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
        case 'B':
            bam_output_path = optarg;
            break;
        case 'b':
            bam_path = optarg;
            break;  
        case 'c':
            barcode_path = optarg;
            break;
        case 'd':
            bam_dir = optarg;
            break;  
        case 'g':
            attr = optarg;
            break; 
        case 'i':
            bam_tag = optarg;
            break;
        case 'm':
            rmsk_path = optarg;
            break;
        case 'r':
            file_path = optarg;
            break; 
        case 'o':
            output_path = optarg;
            break;
        case 't':
            n_threads = MIN(strtol(optarg, NULL, 10), MAX_THREADS);
            break;
        case 'T':
            n_alignment_per_thread = MAX(strtol(optarg, NULL, 10), 1000000);
            break;
        case 'p':
            paired = 1;
            break;
        case 'h':
            count_usage();
            exit(0);
        case 'v':
            count_version();
            exit(0);
        }
    }
    
    if ((bam_output_path != NULL) && (n_threads > 1)) {
        warnf("Outputing BAM with multi-thread option is not yet supported");
        bam_output_path = NULL;
    }

    if ((bam_dir != NULL) && (bam_output_path != NULL)) {
        warnf("Outputing BAM with BAM file directory is not yet supported");
        bam_output_path = NULL;
    }

    if (!bam_path && !bam_dir) {
        fatalf("Either bam_path or bam_dir should be specified.")
    }

    if (bam_path && bam_dir) {
        warn("both bam_path and bam_dir specified. Only considering bam_path");
        bam_dir = NULL;
    }

    ret = molecularCount(
        file_path, 
        rmsk_path, 
        attr, 
        "gene_id",  // this argument is hard-coded for now
        bam_path, 
        bam_dir,
        bam_tag,
        paired,
        barcode_path, 
        output_path, 
        bam_output_path, 
        n_threads,
        n_alignment_per_thread
    );
    return ret;
}

uint8_t molecularCountInternal(
    BGTF *rfp,
    BGTF *tfp,
    BGZF *ofp,
    char *attr_name,
    char *rmsk_attr_name,
    ExpressionMatrix_t *matrix,
    ExpressionMatrix_t *rmsk_matrix,
    char *bam_tag,
    uint8_t bam_paired,
    char *bam_barcode,
    int n_threads,
    tpool_t *p,
    tpool_process_t *q,
    uint32_t n_alignment_per_thread,
    BamFile_t *bamfile)
{
    /* Begin: These declarations are used for multi-thread argument */
    size_t n_par_data = 0;
    int blk;
    ParallelizerResult *result = NULL;
    ParallelizerArg *arg = NULL;
    void **aln_ptr = NULL;
    ArrayList *par_data = NULL;
    /* End: These declarations are used for multi-thread argument */

    BamAlignment *baln = NULL;
    PairedAlignment *paln = NULL;
    char *barcode, *annotation;
    size_t n_alignment = 0, n_annotation = 0;
    if (n_threads == 1)
    {
        if (!bam_paired)
        {
            while (baln = parseOneAlignment(bamfile))
            {
                uint8_t ret;
                if (tfp)
                {
                    ret = annotateBamRmsk(rfp, tfp, baln, attr_name, rmsk_attr_name);
                }
                else
                {
                    ret = annotateBam(rfp, baln, attr_name);
                }
                if (ret == 0)
                {
                    // Successful annotation
                    n_annotation++;

                    if (bam_barcode)
                    {
                        barcode = bam_barcode;
                    }
                    else
                    {
                        barcode = bam_aux_get(bamfile->aln, bam_tag) + 1;
                    }

                    annotation = baln->annotation;
                    if (baln->is_rmsk)
                    {
                        ExpressionMatrixInc(rmsk_matrix, barcode, annotation, 1, 0);
                    }
                    else
                    {
                        ExpressionMatrixInc(matrix, barcode, annotation, baln->spliced_annotation, 0);
                    }
                    if (ofp)
                    {
                        bam_aux_update_str(bamfile->aln, "CT", strlen(annotation), annotation);
                        if (baln->spliced_annotation)
                        {
                            bam_aux_update_str(bamfile->aln, "PT", 8, "spliced");
                        }
                        else
                        {
                            bam_aux_update_str(bamfile->aln, "PT", 10, "unspliced");
                        }
                        bam_write1(ofp, bamfile->aln);
                    }
                }
                else
                {
                    if (ofp)
                    {
                        bam_aux_update_str(bamfile->aln, "PT", 10, "unassigned");
                        bam_write1(ofp, bamfile->aln);
                    }
                };
                n_alignment++;
                if (n_alignment % 1000000 == 0)
                {
                    logf("Count %llu reads", n_alignment);
                }
                BamAlignmentDestroy(baln);
            }
        }
        else
        {
            PairedAlignmentIterator_t *it = PairedAlignmentIteratorInit();

            while (paln = parsePairedAlignment(it, bamfile, bam_tag, bam_barcode))
            {
                uint8_t ret;
                if (rmsk_matrix)
                {
                    ret = annotateBamPairedRmsk(rfp, tfp, paln, attr_name, rmsk_attr_name);
                }
                else
                {
                    ret = annotateBamPaired(rfp, paln, attr_name);
                }
                if (ret == 0)
                {
                    // Successful annotation
                    n_annotation++;
                    barcode = paln->barcode;
                    annotation = paln->annotation;

                    if (paln->is_rmsk)
                    {
                        ExpressionMatrixInc(rmsk_matrix, barcode, annotation, 1, 0);
                    }
                    else
                    {
                        ExpressionMatrixInc(matrix, barcode, annotation, paln->spliced_annotation, 0);
                    }
                }
                n_alignment++;
                if (n_alignment % 1000000 == 0)
                {
                    logf("Count %llu reads", n_alignment);
                }
                HashTableRemove(it->read_table, paln->qname);
                PairedAlignmentDestroy(paln);
            }
            PairedAlignmentIteratorDestroy(it);
        }
    }
    else
    {

        par_data = ArrayListCreate(1000000);
        if (bam_paired)
        {
            PairedAlignmentIterator_t *it = PairedAlignmentIteratorInit();
            while (paln = parsePairedAlignment(it, bamfile, bam_tag, bam_barcode))
            {
                barcode = paln->barcode;
                ArrayListPush(par_data, paln);
                n_par_data++;
                n_alignment++;
                if (n_par_data % n_alignment_per_thread == 0)
                {
                    logf("Count %llu reads", n_alignment);
                    // Begin flushing
                    aln_ptr = par_data->elementList;
                    arg = malloc(sizeof(ParallelizerArg));
                    memset(arg, 0, sizeof(ParallelizerArg));
                    arg->attr_name = attr_name;
                    arg->rmsk_attr_name = rmsk_attr_name;
                    arg->fp = rfp;
                    arg->rmsk = tfp;
                    arg->mtx = matrix;
                    arg->rmsk_mtx = rmsk_matrix;
                    arg->n_bam_alignments = n_alignment_per_thread;
                    arg->aln = aln_ptr;
                    arg->rid = n_alignment;
                    arg->free_ptr = par_data;
                    arg->bam_paired = 1;
                    blk = tpool_dispatch(p, q, molecularCountParallelizer, (void *)arg, NULL, NULL, true);
                    par_data = ArrayListCreate(n_alignment_per_thread);
                    n_par_data = 0;
                }
            }
            logf("Count %llu reads", n_alignment);
            // Begin flushing
            aln_ptr = par_data->elementList;
            arg = malloc(sizeof(ParallelizerArg));
            memset(arg, 0, sizeof(ParallelizerArg));
            arg->attr_name = attr_name;
            arg->rmsk_attr_name = rmsk_attr_name;
            arg->fp = rfp;
            arg->rmsk = tfp;
            arg->mtx = matrix;
            arg->rmsk_mtx = rmsk_matrix;
            arg->n_bam_alignments = n_par_data;
            arg->aln = aln_ptr;
            arg->free_ptr = par_data;
            arg->rid = n_alignment;
            arg->bam_paired = 1;
            blk = tpool_dispatch(p, q, molecularCountParallelizer, (void *)arg, NULL, NULL, true);
            tpool_process_flush(q);
            PairedAlignmentIteratorDestroy(it);
        }
        else
        {
            while (baln = parseOneAlignment(bamfile))
            {
                if (bam_barcode)
                {
                    barcode = bam_barcode;
                }
                else
                {
                    barcode = bam_aux_get(bamfile->aln, bam_tag) + 1;
                }
                baln->barcode = malloc(strlen(barcode) + 1);
                memset(baln->barcode, 0, strlen(barcode) + 1);
                memcpy(baln->barcode, barcode, strlen(barcode));
                ArrayListPush(par_data, baln);
                n_par_data++;
                n_alignment++;

                if (n_par_data % n_alignment_per_thread == 0)
                {

                    logf("Count %llu reads", n_alignment);
                    // Begin flushing
                    aln_ptr = par_data->elementList;

                    arg = malloc(sizeof(ParallelizerArg));
                    memset(arg, 0, sizeof(ParallelizerArg));
                    arg->attr_name = attr_name;
                    arg->rmsk_attr_name = rmsk_attr_name;
                    arg->fp = rfp;
                    arg->rmsk = tfp;
                    arg->mtx = matrix;
                    arg->rmsk_mtx = rmsk_matrix;
                    arg->n_bam_alignments = n_alignment_per_thread;
                    arg->aln = aln_ptr;
                    arg->rid = n_alignment;
                    arg->bam_paired = 0;
                    arg->free_ptr = par_data;
                    /* Block if there are maximum jobs ! */
                    blk = tpool_dispatch(p, q, molecularCountParallelizer, (void *)arg, NULL, NULL, true);
                    par_data = ArrayListCreate(n_alignment_per_thread);
                    n_par_data = 0;
                    // End flushing
                }
            }
            logf("Count %llu reads", n_alignment);
            // Begin flushing
            aln_ptr = par_data->elementList;
            arg = malloc(sizeof(ParallelizerArg));
            memset(arg, 0, sizeof(ParallelizerArg));
            arg->attr_name = attr_name;
            arg->rmsk_attr_name = rmsk_attr_name;
            arg->fp = rfp;
            arg->rmsk = tfp;
            arg->mtx = matrix;
            arg->rmsk_mtx = rmsk_matrix;
            arg->n_bam_alignments = n_par_data;
            arg->aln = aln_ptr;
            arg->rid = n_alignment;
            arg->free_ptr = par_data;
            arg->bam_paired = 0;
            blk = tpool_dispatch(p, q, molecularCountParallelizer, (void *)arg, NULL, NULL, true);
            tpool_process_flush(q);
        }
    }
    return 0;
}