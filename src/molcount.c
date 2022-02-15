#include <getopt.h>
#include <htslib/bgzf.h>
#include "tpool.h"
#include "molcount.h"
#include "sparse.h"

static void bamToLoc( bam_hdr_t *bamHdr, bam1_t *aln, GenomicLocation *loc)
{
    loc->start = aln->core.pos + 1; // left most position of alignment in zero based coordianate (+1)
	loc->chrom = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
    loc->end = loc->start + bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln)); //length of the read
}

void BamAlignmentDestroyHelper(void *e)
{
    destroyGenomicLocation(e, 0);
}

BamAlignment *BamAlignmentInit()
{
    BamAlignment *baln = malloc(sizeof(BamAlignment));
    memset(baln, 0, sizeof(BamAlignment));
    baln->loc = malloc(sizeof(GenomicLocation));
    memset(baln->loc, 0, sizeof(GenomicLocation));
    baln->spliced_loc = LinkedListInit();
    LinkedListSetDealloc(baln->spliced_loc, BamAlignmentDestroyHelper);
    return baln;
}

void BamAlignmentDestroy(BamAlignment *baln)
{
    LinkedListDestroy(GenomicLocation, baln->spliced_loc, 1);
    if (baln->barcode) free(baln->barcode);
    free(baln->loc);
    free(baln);
}

BamAlignment *parseOneAlignment(BamFile_t *bamfile)
{
    if (!(sam_read1(bamfile->fp, bamfile->bam_hdr, bamfile->aln) > 0)) {
        return NULL;
    }
    BamAlignment *baln = BamAlignmentInit();
    baln->spliced_alignment = 0;
    baln->strand = bamfile->aln->core.flag & BAM_FREVERSE ? 0 : 1;
    bamToLoc(bamfile->bam_hdr, bamfile->aln, baln->loc);
    uint32_t *cigar = bam_get_cigar(bamfile->aln);
    uint32_t n_cigar = bamfile->aln->core.n_cigar;
    uint64_t start = bamfile->aln->core.pos + 1;
    uint64_t end = start;
    for (uint32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int ol = bam_cigar_oplen(cigar[i]);
        switch (op)
        {
        case BAM_CREF_SKIP:
            /* the bam alignment is likely to be a splice read. */
            baln->spliced_alignment++;
            GenomicLocation *cur = GenomicLocationNew();
            cur->chrom = bamfile->bam_hdr->target_name[bamfile->aln->core.tid];
            cur->start = start;
            cur->end = end - 1;
            LinkedListAppend(GenomicLocation, baln->spliced_loc, cur);
            start = end + ol;
            end = start;
            break;
        default:
            end = end + ol;
            break;
        }
    }
    if (baln->spliced_alignment > 0) {
        GenomicLocation *cur = GenomicLocationNew();
        cur->chrom = bamfile->bam_hdr->target_name[bamfile->aln->core.tid];
        cur->start = start;
        cur->end = end;
        LinkedListAppend(GenomicLocation, baln->spliced_loc, cur)
    }
    return baln;
}

// TODO: change "gene_name" to user arguments and strandess. 
// This requires passing argument to this callback function
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

static void counterAdd(HashTable *counter, char *key)
{
    if (!(HashTableContainsKey(counter, key)))
    {
        HashTablePut(counter, key, 1);
    }
    int c = HashTableGet(counter, key);
    HashTablePut(counter, key, c + 1);
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
        
        /*
        printf("HEAD ");printGenomicLocation(first); printf("\n");
        printf("TAIL ");printGenomicLocation(last); printf("\n");
        LinkedListForEach(GenomicLocation, record_locs->locs, loc, {
                printGenomicLocation(loc);printf("\n");
        });
        */

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

static ExpressionMatrix_t *ExpressionMatrixDenseInit(size_t n_row, size_t n_col)
{
    logf("Initialize dense matrix %d × %d", n_row, n_col);
    ExpressionMatrix_t *matrix = malloc(sizeof(ExpressionMatrix_t));
    memset(matrix, 0, sizeof(ExpressionMatrix_t));
    matrix->spliced = calloc(n_col * n_row, sizeof(uint16_t));
    matrix->unspliced = calloc(n_col * n_row, sizeof(uint16_t));
    return matrix;
}

static ExpressionMatrix_t *ExpressionMatrixSparseInit(size_t n_row, size_t n_col)
{
    logf("Initialize sparse matrix %d × %d", n_row, n_col);
    ExpressionMatrix_t *matrix = malloc(sizeof(ExpressionMatrix_t));
    memset(matrix, 0, sizeof(ExpressionMatrix_t));
    matrix->spliced = newMatrix(n_row, n_col);
    matrix->unspliced = newMatrix(n_row, n_col);
    return matrix;
}

ExpressionMatrix_t *ExpressionMatrixInit(size_t n_row, size_t n_col, uint8_t sparse)
{
    ExpressionMatrix_t *matrix;
    if (sparse) {
        matrix = ExpressionMatrixSparseInit(n_row, n_col);
    } else {
        matrix = ExpressionMatrixDenseInit(n_row, n_col);
    }
    matrix->is_sparse = sparse;
    matrix->n_row = n_row;
    matrix->n_col = n_col;
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

void ExpressionMatrixSet(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced, uint8_t trylock, uint16_t value)
{
    int64_t col = HashTableGet(matrix->col_index, col_name) - 1;
    int64_t row = HashTableGet(matrix->row_index, row_name) - 1; 
    if (col < 0 || row < 0) return;
    if (matrix->is_sparse) {
        Matrix exp;
        if (spliced) {
            exp = (Matrix) matrix->spliced;
        } else {
            exp = (Matrix) matrix->unspliced;
        }
        if (matrix->n_threads && trylock) {
            while (pthread_mutex_trylock(&matrix->mu) != 0);
        }
        addEntryToIMatrix(exp, row, col, value);
        if (matrix->n_threads && trylock) {
            pthread_mutex_unlock(&matrix->mu);
        }
    } else {
        uint16_t *spliced_exp, *unspliced_exp;;
        spliced_exp = (uint16_t *) matrix->spliced;
        unspliced_exp = (uint16_t *) matrix->unspliced;
        if (matrix->n_threads && trylock) {
            while (pthread_mutex_trylock(&matrix->mu) != 0);
        }
        if ((spliced_exp[row * matrix->n_col + col] == 0) && 
            (unspliced_exp[row * matrix->n_col + col] == 0) && 
            (value != 0)) {
            matrix->n_non_zero++;
        }
        if ((spliced_exp[row * matrix->n_col + col] != 0) && 
            (unspliced_exp[row * matrix->n_col + col] == 0) && 
            (value != 0) && spliced) {
            matrix->n_non_zero--;
            }
        if ((spliced_exp[row * matrix->n_col + col] == 0) && 
            (unspliced_exp[row * matrix->n_col + col] != 0) && 
            (value != 0) && !spliced) {
            matrix->n_non_zero--;
            }

        if (spliced) spliced_exp[row * matrix->n_col + col] = value;
        else unspliced_exp[row * matrix->n_col + col] = value;
        if (matrix->n_threads) {
            pthread_mutex_unlock(&matrix->mu);
        }
    }
}



void ExpressionMatrixInc(ExpressionMatrix_t *matrix, char *row_name, char *col_name, uint8_t spliced, uint8_t trylock)
{
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
        Matrix exp;
        if (spliced) {
            exp = (Matrix) matrix->spliced;
        } else {
            exp = (Matrix) matrix->unspliced;
        }
        if (matrix->n_threads && trylock) {
            while (pthread_mutex_trylock(&matrix->mu) != 0);
        }
        accumulateEntryInIMatrix(exp, row, col, 1);
        if (matrix->n_threads && trylock) {
            pthread_mutex_unlock(&matrix->mu);
        }
    } else {
        uint16_t *spliced_exp, *unspliced_exp;;
        spliced_exp = (uint16_t *) matrix->spliced;
        unspliced_exp = (uint16_t *) matrix->unspliced;
        if (matrix->n_threads && trylock) {
            while (pthread_mutex_trylock(&matrix->mu) != 0);
        }
        if ((spliced_exp[row * matrix->n_col + col] == 0) && 
            (unspliced_exp[row * matrix->n_col + col] == 0)) {
            matrix->n_non_zero++;
        }
        if (spliced) {
            spliced_exp[row * matrix->n_col + col]++;
        } else {
            unspliced_exp[row * matrix->n_col + col]++;
        }
        if (matrix->n_threads && trylock) {
            pthread_mutex_unlock(&matrix->mu);
        }
    }
}

void ExpressionMatrixDestroy(ExpressionMatrix_t *matrix)
{   
    if (matrix->is_sparse) {

    } else {
        free(matrix->spliced);
        free(matrix->unspliced);
        free(matrix);
    }
}

uint8_t ExpressionMatrixToFile(ExpressionMatrix_t *matrix, FILE *spliced_stream, FILE *unspliced_stream)
{
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
    BamAlignment **baln;
    size_t n_bam_alignments;
    BGTF *fp;
    ExpressionMatrix_t *mtx;
    char *attr_name;
} ParallelizerArg;

typedef struct {
    char *barcode;
    char *annotation;
    uint8_t spliced_annotation;
} ParallelizerResult;

void ParallelizerArgClean(void *arg)
{
    ParallelizerArg *pa = (ParallelizerArg *)arg;
    for (size_t i = 0; i < pa->n_bam_alignments; ++i) {
        BamAlignmentDestroy(pa->baln[i]);
    }
    free(pa);
}

void molecularCountParallelizer(void *arg) {
    ParallelizerArg *pa = (ParallelizerArg *)arg;
    ParallelizerResult *result;
    ArrayList *par_result = ArrayListCreate(pa->n_bam_alignments);
    ArrayListSetDeallocationFunction(par_result, free);
    for (size_t i = 0; i < pa->n_bam_alignments; ++i) {
        if (annotateBam(pa->fp, pa->baln[i], pa->attr_name) == 0) {
            char *annotation = pa->baln[i]->annotation;
            result = malloc(sizeof(ParallelizerResult));
            result->barcode = pa->baln[i]->barcode;
            result->annotation = annotation;
            result->spliced_annotation = pa->baln[i]->spliced_annotation;
            ArrayListPush(par_result, result);
        }
    }
    while (pthread_mutex_trylock(&pa->mtx->mu) != 0);
    ArrayListForEach(par_result, result, {
        ExpressionMatrixInc(pa->mtx, result->barcode, result->annotation, result->spliced_annotation, 0);
    })
    pthread_mutex_unlock(&pa->mtx->mu);
    ArrayListDestroy(par_result);
    ParallelizerArgClean(arg);
}

uint8_t molecularCount(char *bgtf_path, 
                    char *attr_name, 
                    char *bam_path,
                    char *bam_tag,
                    char *barcode_path,
                    char *output_path,
                    char *bam_output_path,
                    int n_threads)
{
    /* Begin: These declarations are used for multi-thread argument */
    tpool_t *p;
    tpool_process_t *q;
    if (n_threads > 1) {
        pthread_setconcurrency(2);
        p = tpool_init(n_threads);
        q = tpool_process_init(p, 16, true);
    }
    size_t n_par_data = 0;
    int blk;
    ParallelizerResult *result;
    ParallelizerArg *arg;
    BamAlignment **baln_ptr;
    size_t partition_size;
    ArrayList *par_data;
    /* End: These declarations are used for multi-thread argument */

    BGTF *rfp = BGTFInit();
    BGTFLoad(rfp, bgtf_path, attr_name);
    BamFile_t *bamfile = BamOpen(bam_path);
    BamAlignment *baln;
    HashTable *barcodeIndex = readBarcodeFile(barcode_path);
    BGZF *ofp;
    if (bam_output_path) {
        ofp = bgzf_open(bam_output_path, "wb+");
        bam_hdr_write(ofp, bamfile->bam_hdr);
    }
    ExpressionMatrix_t *matrix = ExpressionMatrixInit(HashTableSize(barcodeIndex), HashTableSize(rfp->attr_table), 0);
    ExpressionMatrixInitSetColIndex(matrix, rfp->attr_table);
    ExpressionMatrixInitSetRowIndex(matrix, barcodeIndex);
    if (n_threads > 1) {
        ExpressionMatrixSetMt(matrix);
    }
    char *barcode, *annotation;
    char *k;
    size_t attr_idx;
    ArrayList *sort_arr;
    sort_arr = ArrayListCreate(HashTableSize(rfp->attr_table));
    HashTableForEach(rfp->attr_table, k, attr_idx, {
        ArrayListPush(sort_arr, k);
    })
    if ((output_path[strlen(output_path)-1]) != '/') {
        char *output_path_ = output_path;
        output_path =  malloc(strlen(output_path_) + 2);
        memset(output_path, 0, strlen(output_path_) + 2);
        memcpy(output_path, output_path_, strlen(output_path_));
        output_path[strlen(output_path_)] = '/';
    }
    char *gene_name_output_path = malloc(strlen(output_path) + 13);
    memset(gene_name_output_path, 0, strlen(output_path) + 13);
    memcpy(gene_name_output_path, output_path, strlen(output_path));
    memcpy(gene_name_output_path + strlen(output_path), "features.tsv", 12);
    FILE *gene_name_out;
    if (!(gene_name_out = fopen(gene_name_output_path, "w+"))) {
        fatalf("failed to open %s", gene_name_output_path);
    }
    
    logf("Output features names to %s", gene_name_output_path);
    ArrayListSort(sort_arr, strcmp);
    ArrayListForEach(sort_arr, k, {
        fprintf(gene_name_out, "%s\n", k);
    })
    ArrayListDestroy(sort_arr);
    fclose(gene_name_out);
    free(gene_name_output_path);

    sort_arr = ArrayListCreate(HashTableSize(rfp->attr_table));
    HashTableForEach(barcodeIndex, k, attr_idx, {
        ArrayListPush(sort_arr, k);
    })
    char *barcode_output_path = malloc(strlen(output_path) + 13);
    memset(barcode_output_path, 0, strlen(output_path) + 13);
    memcpy(barcode_output_path, output_path, strlen(output_path));
    memcpy(barcode_output_path + strlen(output_path), "barcodes.tsv", 12);
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
    log("Start counting molecules");
    size_t n_alignment = 0, n_annotation = 0;
    if (n_threads == 1) {
        while (baln = parseOneAlignment(bamfile)) {
            // printf("read name %s\n", bam_get_qname(bamfile->aln));
            if (annotateBam(rfp, baln, attr_name) == 0) {
                // Successful annotation
                n_annotation++;
                barcode = bam_aux_get(bamfile->aln, bam_tag);
                barcode = barcode + 1;
                annotation = baln->annotation;
                ExpressionMatrixInc(matrix, barcode, annotation, baln->spliced_annotation, 0);
                if (bam_output_path) {
                    bam_aux_update_str(bamfile->aln, "CT", strlen(annotation), annotation);
                    if (baln->spliced_annotation) {
                        bam_aux_update_str(bamfile->aln, "PT", 8, "spliced");
                    } else {
                        bam_aux_update_str(bamfile->aln, "PT", 10, "unspliced");
                    }
                    bam_write1(ofp, bamfile->aln);
                }
            } else {
                if (bam_output_path) {
                    bam_aux_update_str(bamfile->aln, "PT", 10, "unassigned");
                    bam_write1(ofp, bamfile->aln);
                }
            }; 
            n_alignment++;
            if (n_alignment % 1000000 == 0) {
                logf("Count %llu reads", n_alignment);
            }
            BamAlignmentDestroy(baln);
        }
    } else {
        par_data = ArrayListCreate(1000000);
        while (baln = parseOneAlignment(bamfile)) {
            barcode = bam_aux_get(bamfile->aln, bam_tag);
            barcode = barcode + 1;
            baln->barcode = malloc(strlen(barcode)+1);
            memset(baln->barcode, 0, strlen(barcode)+1);
            memcpy(baln->barcode, barcode, strlen(barcode));
            ArrayListPush(par_data, baln);
            n_par_data++;
            n_alignment++;
            
            if (n_par_data % (1000000 * n_threads) == 0) {
                
                logf("Count %llu reads", n_alignment);
                // Begin flushing
                ArrayListPartition(par_data, n_threads, baln_ptr, partition_size, {
                    arg = malloc(sizeof(ParallelizerArg));
                    arg->attr_name = attr_name;
                    arg->fp = rfp;
                    arg->mtx = matrix;
                    arg->n_bam_alignments = partition_size;
                    arg->baln = baln_ptr;
                    blk = tpool_dispatch(p, q, molecularCountParallelizer, (void *)arg, NULL, NULL, true);
                });
                tpool_process_flush(q);
                ArrayListDestroy(par_data);
                par_data = ArrayListCreate(1000000);
                n_par_data = 0;
                // End flushing
            }
        }
    }
    logf("Count %llu reads", n_alignment);
    // Begin flushing
    ArrayListPartition(par_data, n_threads, baln_ptr, partition_size, {
        logf("partition_size = %llu", partition_size);
        arg = malloc(sizeof(ParallelizerArg));
        arg->attr_name = attr_name;
        arg->fp = rfp;
        arg->mtx = matrix;
        arg->n_bam_alignments = partition_size;
        arg->baln = baln_ptr;
        blk = tpool_dispatch(p, q, molecularCountParallelizer, (void *)arg, NULL, NULL, true);
    });
    tpool_process_flush(q);
    ArrayListDestroy(par_data);
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
    log("Writing count matrix");
    ExpressionMatrixToFile(matrix, spliced, unspliced);
    free(mtx_output_path);
    fclose(spliced);
    fclose(unspliced);
    if (n_threads > 1) {
        ExpressionMatrixUnsetMt(matrix);
        tpool_process_destroy(q);
        tpool_destroy(p);
        pthread_exit(NULL);
    }
    ExpressionMatrixDestroy(matrix);
    destroyBarcodeTable(barcodeIndex);
    BamClose(bamfile);
    if (bam_output_path) {
        bgzf_close(ofp);
    }
    return 0;
}


/* Command-line functions */
void count_usage()
{
    
    return;
}
void count_version()
{
    return;
}

int count_main(int argc, char *argv[])
{
    int c;
    char *file_path = NULL;
    char *output_path = NULL;
    char *bam_path = NULL;
    char *bam_output_path = NULL;
    char *bam_tag = NULL;
    char *barcode_path = NULL;
    char *attr;
    int ret;
    int n_threads = 1;
    while (1)
    {
        static struct option long_options[] =
            {   
                {"bamOutput", required_argument, 0, 'B'},
                {"bamFile", required_argument, 0, 'b'},
                {"cellBarcode", required_argument, 0, 'c'},
                {"bamTag", required_argument, 0, 'i'},
                {"attrName", required_argument, 0, 'g'},
                {"gtf", required_argument, 0, 'r'},
                {"output", required_argument, 0, 'o'},
                {"nthreads", optional_argument, 0, 'p'},
                {"help", no_argument, NULL, 'h'},
                {"version", no_argument, NULL, 'v'},
                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;
        c = getopt_long(argc, argv, "B:b:c:g:hi:o:p:r:v", long_options, &option_index);

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
        case 'g':
            attr = optarg;
            break; 
        case 'i':
            bam_tag = optarg;
            break;
        case 'r':
            file_path = optarg;
            break; 
        case 'o':
            output_path = optarg;
            break;
        case 'p':
            n_threads = MIN(strtol(optarg, NULL, 10), MAX_THREADS);
            break;
        case 'h':
            count_usage();
            exit(0);
        case 'v':
            count_version();
            exit(0);
        }
    }
    logf("Using (%d/%d) threads", n_threads, MAX_THREADS);
    if ((bam_output_path != NULL) && (n_threads > 1)) {
        warnf("Outputing BAM with multi-thread option is not yet supported");
        bam_output_path = NULL;
    }
    ret = molecularCount(file_path, attr, bam_path, bam_tag, barcode_path, output_path, bam_output_path, n_threads);
    return ret;
}
