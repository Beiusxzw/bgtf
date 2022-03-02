#ifndef BGTF_BAM
#define BGTF_BAM 

#include <htslib/sam.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bgtf.h"
#include "list.h"

typedef struct BamFile {
    samFile *fp;
    bam_hdr_t *bam_hdr;
    bam1_t *aln;
} BamFile_t;

typedef struct BamAlignment {
    GenomicLocation *loc;
    uint8_t strand;
    uint8_t spliced_alignment;
    uint8_t is_rmsk;
    uint8_t spliced_annotation;
    LinkedList *spliced_loc;
    char *annotation;
    char *barcode;
    uint16_t flag;
} BamAlignment;

typedef struct PairedAlignmentIterator {
    HashTable *read_table;
} PairedAlignmentIterator_t;

PairedAlignmentIterator_t *PairedAlignmentIteratorInit();
void PairedAlignmentIteratorDestroy(PairedAlignmentIterator_t *it);

typedef struct PairedAlignment {
    BamAlignment *mate1;
    BamAlignment *mate2;
    uint8_t paired;
    char *qname;
    char *annotation;
    uint8_t _free_barcode;
    char *barcode;
    uint8_t is_rmsk;
    uint8_t spliced_annotation;
} PairedAlignment;

void PairedAlignmentDestroy(PairedAlignment *paired_aln);

static inline BamFile_t *BamOpen(char *path)
{
    BamFile_t *bamfile = (BamFile_t *)malloc(sizeof(BamFile_t));
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
PairedAlignment *parsePairedAlignment(
    PairedAlignmentIterator_t *peIter, 
    BamFile_t *bamfile, 
    char *bam_tag,
    char *bam_barcode);
uint8_t assignFeature(const uint64_t *rect, void *item, void *udata);
uint8_t annotateBam(BGTF *fp, BamAlignment *baln, char *attr_name);
void bamToLoc( bam_hdr_t *bamHdr, bam1_t *aln, GenomicLocation *loc);


#endif