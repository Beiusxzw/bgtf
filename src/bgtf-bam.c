#include "bgtf-bam.h"

PairedAlignmentIterator_t *PairedAlignmentIteratorInit()
{
    PairedAlignmentIterator_t *peIter = (PairedAlignmentIterator_t *)malloc(sizeof(PairedAlignmentIterator_t));
    memset(peIter, 0, sizeof(PairedAlignmentIterator_t));
    peIter->read_table = HashTableCreate(0x10000);
    HashTableSetKeyComparisonFunction(peIter->read_table, strcmp);
    HashTableSetHashFunction(peIter->read_table, HashTableStringHashFunction);
    return peIter;
}

void PairedAlignmentIteratorDestroy(PairedAlignmentIterator_t *it)
{
    HashTableDestroy(it->read_table);
    free(it);
}

static PairedAlignment *PairedAlignmentInit(uint8_t paired)
{
   PairedAlignment *paired_aln = (PairedAlignment *)malloc(sizeof(PairedAlignment));
   memset(paired_aln, 0, sizeof(PairedAlignment));
   paired_aln->paired = paired;
   paired_aln->_free_barcode = 1;
   return paired_aln;
}

void PairedAlignmentDestroy(PairedAlignment *paired_aln)
{
    if (paired_aln->mate1) BamAlignmentDestroy(paired_aln->mate1);
    if (paired_aln->mate2) BamAlignmentDestroy(paired_aln->mate2);
    free(paired_aln->qname);
    if (paired_aln->_free_barcode) free(paired_aln->barcode);
    free(paired_aln);
}

PairedAlignment *parsePairedAlignment(
    PairedAlignmentIterator_t *peIter, 
    BamFile_t *bamfile, 
    char *bam_tag,
    char *bam_barcode)
{
    PairedAlignment *paired_aln;
    uint8_t z = 1;
    BamAlignment *aln;
    while (z) {
        aln = parseOneAlignment(bamfile);
        if (!aln) {
            return NULL;
        }
        if ((bamfile->aln->core.flag & BAM_FPROPER_PAIR) == 0) {
            BamAlignmentDestroy(aln);
        } else {
            z = 0;
        }
    }

    char *_read_name = bam_get_qname(bamfile->aln);
    char *read_name = malloc(strlen(_read_name)+1);
    memset(read_name, 0, strlen(_read_name)+1);
    memcpy(read_name, _read_name, strlen(_read_name));
    char *barcode;

    if (bam_barcode) {
        /* if the bam barcode is provided, we will not extract
        barcode from bam*/
        barcode = bam_barcode;
    } else {
        char *_barcode = bam_aux_get(bamfile->aln, bam_tag) + 1;
        barcode = malloc(strlen(_barcode)+1);
        memset(barcode, 0, strlen(_barcode)+1);
        memcpy(barcode, _barcode, strlen(_barcode));
    }


    if (!(HashTableContainsKey(peIter->read_table, read_name))) {
        
        paired_aln = PairedAlignmentInit(1);
        if ((bamfile->aln->core.flag & BAM_FREAD1) != 0) {
            paired_aln->mate1 = aln;
        } else {
            paired_aln->mate2 = aln;
        }
        paired_aln->qname = read_name;
        paired_aln->barcode = barcode;
        
        if (bam_barcode) paired_aln->_free_barcode = 0;
        HashTablePut(peIter->read_table, read_name, paired_aln);

        return parsePairedAlignment(peIter, bamfile, bam_tag, bam_barcode);
    } else {
       
        paired_aln = (PairedAlignment *)HashTableGet(peIter->read_table, read_name);

        if ((bamfile->aln->core.flag & BAM_FREAD1) != 0) {
            assert(!paired_aln->mate1);
            paired_aln->mate1 = aln;
        } else {
            assert(!paired_aln->mate2);
            paired_aln->mate2 = aln;
        }
        free(read_name);
        if (!bam_barcode) free(barcode);
        return paired_aln;
    }
}

BamAlignment *parseOneAlignment(BamFile_t *bamfile)
{
    uint8_t z = 1;
    while (z) {
        if (!(sam_read1(bamfile->fp, bamfile->bam_hdr, bamfile->aln) > 0)) {
            return NULL;
        }
        if (bamfile->aln->core.tid < 0) {
            /* Check whether the BAM contig is valid. This should not happen.*/
            // warnf("The contig id from this BamRecord is invalid. Please check your BAM file!");
        } else {
            z = 0;
        }
    }
    BamAlignment *baln = BamAlignmentInit();
    baln->flag = bamfile->aln->core.flag;
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

void bamToLoc( bam_hdr_t *bamHdr, bam1_t *aln, GenomicLocation *loc)
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
