#include "bgtf.h"
#include "bgtfIO.h"
#include "hashtable.h"
#include "endianness.h"

/**
 * compress functions derived from bgzip
 */

static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";


// get the compress level from the mode string
static inline int mode2level(const char *__restrict mode)
{
	int i, compress_level = -1;
	for (i = 0; mode[i]; ++i)
		if (mode[i] >= '0' && mode[i] <= '9') break;
	if (mode[i]) compress_level = (int)mode[i] - '0';
	if (strchr(mode, 'u')) compress_level = -2;
	return compress_level;
}


static uint8_t BGTFCompress(void *_dst, int *dlen, void *src, int slen, int level)
{
	uint32_t crc;
	z_stream zs;
	uint8_t *dst = (uint8_t*)_dst;

	// compress the body
	zs.zalloc = NULL; zs.zfree = NULL;
	zs.next_in  = (Bytef*)src;
	zs.avail_in = slen;
	zs.next_out = dst + BGTF_BLOCK_HEADER_LENGTH;
	zs.avail_out = *dlen - BGTF_BLOCK_HEADER_LENGTH - BGTF_BLOCK_FOOTER_LENGTH;
	if (deflateInit2(&zs, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) return -1; // -15 to disable zlib header/footer
	if (deflate(&zs, Z_FINISH) != Z_STREAM_END) return -1;
	if (deflateEnd(&zs) != Z_OK) return -1;
	*dlen = zs.total_out + BGTF_BLOCK_HEADER_LENGTH + BGTF_BLOCK_FOOTER_LENGTH;
	// write the header
	memcpy(dst, g_magic, BGTF_BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
	packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes
	// write the footer
	crc = crc32(crc32(0L, NULL, 0L), (Bytef*)src, slen);
	packInt32((uint8_t*)&dst[*dlen - 8], crc);
	packInt32((uint8_t*)&dst[*dlen - 4], slen);
	return 0;  
}

static uint8_t deflate_block(BGTFRWBuffer_t *buffer, int block_length)
{
	int comp_size = BGTF_MAX_BLOCK_SIZE;
	if (BGTFCompress(buffer->compressed_block, 
                    &comp_size, 
                    buffer->uncompressed_block, 
                    block_length, 
                    buffer->compress_level) != 0) {
		buffer->errcode |= BGTF_ERR_ZLIB;
		return -1;
	}
	buffer->block_offset = 0;
	return comp_size;
}

static uint8_t inflate_block(BGTFRWBuffer_t *buffer, int block_length)
{
	z_stream zs;
	zs.zalloc = NULL;
	zs.zfree = NULL;
	zs.next_in = (Bytef*)buffer->compressed_block + 18;
	zs.avail_in = block_length - 16;
	zs.next_out = (Bytef*)buffer->uncompressed_block;
	zs.avail_out = BGTF_MAX_BLOCK_SIZE;

	if (inflateInit2(&zs, -15) != Z_OK) {
		buffer->errcode |= BGTF_ERR_ZLIB;
		return -1;
	}
	if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
		inflateEnd(&zs);
		buffer->errcode |= BGTF_ERR_ZLIB;
		return -1;
	}
	if (inflateEnd(&zs) != Z_OK) {
		buffer->errcode |= BGTF_ERR_ZLIB;
		return -1;
	}
	return zs.total_out;
}

/*******
 * Read *
 ******/

BGTF *GTFRead(char *path, char *attr_name)
{
    log("Loading GTF file");
    FILE *fp;
    char buf[BGTF_IO_BUFFER_SIZE];
    char readbuf[BGTF_IO_BUFFER_SIZE];
    uint8_t binary_checked = 0, sorted_checked = 0;
    if (!(fp = fopen(path, "r")))
    {
        fprintf(stderr, "Failed to open file");
        exit(1);
    }
    BGTF *file = BGTFInit();
    BGTFRecord_t *curr = NULL, *head = NULL, *tail = NULL;
    // Initializing the array list storing the records
    ArrayList *records = ArrayListCreate(0xFF);
    BGTFRecord_t *record;
    // attach the record deallocator to the array list
    ArrayListSetDeallocationFunction(records, BGTFRecordDestroy);
    BGTFRecord_t *prev = NULL;
    while (fgets(buf, BGTF_IO_BUFFER_SIZE, fp) != NULL)
    {
        if (!binary_checked)
        {
            if (is_binary(buf, strlen(buf)))
            {
                fatal("The GTF File is binary. Try reading it as BGTF File");
            };
            binary_checked = 1;
        }
        if (!(record = parseOneRecord(readbuf, buf)))
            continue;
        if (!sorted_checked && prev)
        {
            if (RecordCompare(prev, record) > 0)
            {
                sorted_checked = 1;
                warn("The GTF file is not sorted.")
            }
        }
        prev = record;
        ArrayListPush(records, record);
    }
    file->records = records;
    if (sorted_checked)
        BGTFSortRecord(file);
    BGTFBuildIndex(file, attr_name);
    
    char *chrom;
    BGTFRecordIdx *value;
    RTree_t *rt;
    /*
    HashTableForEach(
        file->index->chrom_idx,
        chrom,
        value, {
            logf("%llu, %s, %x, %x", file->index->chrom_idx->hashFunction(chrom),
                 chrom,
                 ((BGTFRecordIdx *)value)->head,
                 ((BGTFRecordIdx *)value)->tail);
        })
    */
        // BGTFflushRecord(stdout, file);
    log("Loading Finished");
    return file;
}

void BGTFAttrfprint (FILE *__stream__, RecordAttrKV_t *attr, char *k)
{
    enum BGTF_ATTR_MAGIC magic = HashTableStringHashFunction(k);
    switch (magic)
    {
    case GENE_ID:
        fprintf(__stream__, "%s", attr->gene_id);
        break;
    case GENE_VERSION:
        fprintf(__stream__, "%d", attr->gene_version);
        break;
    case GENE_NAME:
        fprintf(__stream__, "%s", attr->gene_name);
        break;
    case TRANSCRIPT_ID:
        fprintf(__stream__, "%s", attr->transcript_id);
        break;
    case TRANSCRIPT_VERSION:
        fprintf(__stream__, "%s", attr->transcript_version);
        break;
    case TRANSCRIPT_NAME:
        fprintf(__stream__, "%s", attr->transcript_name);
        break;
    case EXON_NUMBER:
        fprintf(__stream__, "%d", attr->exon_number);
        break;
    case EXON_ID:
        fprintf(__stream__, "%s", attr->exon_id);
        break;
    case FAMILY_ID:
       fprintf(__stream__, "%s", attr->family_id);
        break;
    case CLASS_ID:
        fprintf(__stream__, "%s", attr->class_id);
        break;
    default:
        break;
    }
}

void BGTFfprint(FILE *__stream__, BGTFRecord_t *record, char *k)
{
    char c;
    if (record->strand)
        c = '+';
    else
        c = '-';
    fprintf(__stream__, "%s\t%s\t%s\t%llu\t%llu\t.\t%c\t.\t",
            record->chrom,
            record->source,
            record->type,
            record->start,
            record->end,
            c);
    if (!k) {fprintf(__stream__, "\n");return;}
    BGTFAttrfprint(__stream__, record->attrs, k);
    fprintf(__stream__, "\n");
}

void printGenomicLocation(GenomicLocation *loc)
{
    fprintf(stdout, "%s:%llu-%llu", loc->chrom, loc->start, loc->end);
}

void printBamAlignment(struct BamAlignment *t)
{
    printGenomicLocation(t->loc);
    if (t->spliced_alignment > 0) {
        printf(" [S:%d] \n", t->spliced_alignment);
    } else {
        printf(" [U] \n");
    }
    GenomicLocation *loc = NULL;
    if (t->spliced_alignment > 0) 
    {
        LinkedListForEach(GenomicLocation, t->spliced_loc, loc, {
            printf("   ");
            printGenomicLocation(loc);
            printf("\n");
        })
        printf("   ");
        printf("\n");
    }
}


void BGTFflushRecord(FILE *__stream__, BGTF *fp, char *k)
{
    BGTFRecord_t *record = NULL;
    
    ArrayListForEach(fp->records, record, {
        BGTFfprint(stdout, record, k);
    })
}

static inline uint32_t serializeString(void *dest, char *src, uint32_t magic)
{
    uint32_t c = strlen(src);
    memcpy(dest, &c, 4);
    memcpy(dest + 4, &magic, 4);
    memcpy(dest + 8, src, c);
    return c + 8;
}

static inline uint32_t deserializeString(RecordAttrKV_t *attr, char *src)
{
    uint32_t c;
    uint32_t magic;
    memcpy(&c, src, 4);
    memcpy(&magic, src + 4, 4);
    char *v = (char *)malloc(c+1);
    memset(v, 0, c+1);
    memcpy(v, src + 8, c);
    setRecordAttrKV_helper(attr, magic, v);
    return c + 8;
}


uint32_t RecordAttrSerialize(void *dest, int *dlen, RecordAttrKV_t *attr, int *slen)
{
    uint32_t p = 4;
    memcpy(dest + p, &(attr->exon_number), 2);
    memcpy(dest + p + 2, &(attr->gene_version), 2);
    p = 8;
    if (attr->gene_id) p += serializeString(dest + p, attr->gene_id, GENE_ID);
    if (attr->gene_name) p += serializeString(dest + p, attr->gene_name, GENE_NAME);
    if (attr->transcript_id) p += serializeString(dest + p, attr->transcript_id, TRANSCRIPT_ID);
    if (attr->transcript_name) p += serializeString(dest + p, attr->transcript_name, TRANSCRIPT_NAME);
    if (attr->exon_id) p += serializeString(dest + p, attr->exon_id, EXON_ID);
    if (attr->family_id) p += serializeString(dest + p, attr->family_id, FAMILY_ID);
    if (attr->class_id) p += serializeString(dest + p, attr->class_id, CLASS_ID);
    memcpy(dest, &p, 4);  // the total length of record string
    return p;
}


static inline void RecordAttrDeserialize(char *src,  RecordAttrKV_t *attr)
{
    int c;
    uint32_t slen;
    uint32_t p = 4;
    memcpy(&(attr->exon_number), src + p, 2);
    memcpy(&(attr->gene_version), src + p + 2,  2);
    memcpy(&c, src, 4); // the total length of record string
    p = 8;
    while (c > 8) {
        slen = deserializeString(attr, src + p);
        p += slen;
        c -= slen;
    }
}

static int BGTFSerializeRecord(void *dest, int *dlen, BGTFRecord_t *record, int *slen)
{
    assert(dlen > slen);
    uint32_t c;
    uint32_t p = 8;
    // write chromosome
    c = strlen(record->chrom) + 1;
    memcpy(dest + p, record->chrom, strlen(record->chrom) + 1);
    memcpy(dest + p - 4, &c, 4);
    p += c + 4;
    // write source
    c = strlen(record->source) + 1;
    memcpy(dest + p, record->source, strlen(record->source) + 1);
    memcpy(dest + p - 4, &c, 4);
    p += c + 4;
    // write type
    c = strlen(record->type) + 1;
    memcpy(dest + p, record->type, strlen(record->type) + 1);
    memcpy(dest + p - 4, &c, 4);
    p += c;
    // write start and end
    memcpy(dest + p, &(record->start), 8);
    memcpy(dest + p + 8, &(record->end), 8);
    p += 16;
    // we do not save column 6 and 8 here
    // write strand
    memcpy(dest + p, &(record->strand), 1);
    p += 1;
    // read record
    p += RecordAttrSerialize(dest + p, dlen, record->attrs, slen);
    memcpy(dest, &p, 4);
    *dlen = p;
    return 0;
}

static BGTFRecord_t *BGTFDeserializeRecord(void *src)
{
    BGTFRecord_t *record = (BGTFRecord_t *)malloc(sizeof(BGTFRecord_t));
    record->attrs = (RecordAttrKV_t *)malloc(sizeof(RecordAttrKV_t));
    memset(record->attrs, 0, sizeof(RecordAttrKV_t));
    uint32_t slen;
    uint32_t p;
    uint32_t c;

    memcpy(&slen, src, 4);
    // read chromosome
    memcpy(&c, src + 4, 4);
    record->chrom = malloc(c);
    memcpy(record->chrom, src + 8, c);
    p = 8 + c;
    // read source
    memcpy(&c, src + p, 4);
    p += 4;
    record->source = malloc(c);
    memcpy(record->source, src + p, c);
    p += c;
    // read type
    memcpy(&c, src + p, 4);
    p += 4;
    record->type = malloc(c);
    memcpy(record->type, src + p, c);
    p += c;
    // read start and end
    memcpy(&(record->start), src + p, 8);
    memcpy(&(record->end), src + p + 8, 8);
    p += 16;
    // we do not save column 6 and 8 here
    // read strand
    memcpy(&(record->strand), src + p, 1);
    p += 1;
    RecordAttrDeserialize(src + p, record->attrs);
    return record;
}

/********
 * Write *
 *******/

/**
 * @brief 
 * 
 * @param mode 1 for read and 0 for write
 * @param compress_level 
 * @return BGTFRWBuffer_t* 
 */
BGTFRWBuffer_t *BGTFRWBufferInit(uint8_t mode, int compress_level)
{
    BGTFRWBuffer_t *buffer;
    if (!(buffer = malloc(sizeof(BGTFRWBuffer_t)))) {
        fatal("BGTFRWBuffer Initialization failed");
    }
    if (mode) {
        buffer->is_write = 0;
        // TODO: Check whether the file is actually compressed
        buffer->is_compressed = 1;
        buffer->uncompressed_block = malloc(BGTF_MAX_BLOCK_SIZE);
        buffer->compressed_block = malloc(BGTF_MAX_BLOCK_SIZE);
    } else {
        buffer->is_write = 1;
        if (compress_level == -2) {
            buffer->is_compressed = 0;
            return buffer; // not compressing 
        }
        buffer->is_compressed = 1;
        buffer->compress_level = compress_level < 0 ? Z_DEFAULT_COMPRESSION : compress_level;
        buffer->uncompressed_block = malloc(BGTF_MAX_BLOCK_SIZE);
        buffer->compressed_block = malloc(BGTF_MAX_BLOCK_SIZE);
        buffer->compress_level = buffer->compress_level > 9 ? Z_DEFAULT_COMPRESSION : compress_level;
    }
    return buffer;
}

void BGTFRWBufferDestroy(BGTFRWBuffer_t *buffer)
{
    free(buffer->compressed_block);
    free(buffer->uncompressed_block);
    free(buffer);
}

BGTFRW *BGTFRWInit(char *path, char *mode)
{
    FILE *fp = NULL;
    BGTFRW *write = (BGTFRW *)malloc(sizeof(BGTFRW));
    assert(compressBound(BGTF_BLOCK_SIZE) < BGTF_MAX_BLOCK_SIZE);
    if (strchr(mode, 'r')) {
        fp = fopen(path, mode);
        if (!fp) fatal("file open failed");
        write->fp = fp;
        write->buffer = BGTFRWBufferInit(0, mode2level(mode));
    } else if (strchr(mode, 'w') || strchr(mode, 'a')) {
        fp = fopen(path, mode);
        if (!fp) fatal("file open failed");
        write->fp = fp;
        write->buffer = BGTFRWBufferInit(0, mode2level(mode));
    } else {
        fatal("Unknown file opening mode");
    }
    return write;
}

void BGTFRWDestroy(BGTFRW *write)
{
    fclose(write->fp);
    BGTFRWBufferDestroy(write->buffer);
    free(write);
}

static uint8_t writeAtPos(void *ptr, size_t sz, size_t nmemb, size_t pos, FILE *fp)
{
    size_t curpos = ftell(fp);
    if (fseek(fp, pos, SEEK_SET))
        return 1;
    if (fwrite(ptr, sz, nmemb, fp) != nmemb)
        return 2;
    if (fseek(fp, curpos, SEEK_SET))
        return 3;
    return 0;
}


uint8_t BGTFWriteChromList(FILE *__stream__, BGTF *file)
{
    HashTable *chrom_table = file->index->chrom_idx;
    char *chrom, *value;
    uint16_t chrom_len;
    HashTableForEach(file->index->chrom_idx, chrom, value, {
        chrom_len = strlen(chrom);
        if (fwrite(&chrom_len, sizeof(uint16_t), 1, __stream__) != 1)
            return 1;
        if (!(fwrite(chrom, sizeof(uint8_t), strlen(chrom), __stream__) == strlen(chrom)))
            return 1;
    }) return 0;
}

uint8_t BGTFWriteHdr(BGTF *file)
{
    FILE *__stream__ = file->write->fp;
    uint32_t magic = BGTF_MAGIC;
    // 58 bytes of zeros
    void *p = calloc(58, sizeof(uint8_t));

    uint16_t two = 4;
    if (!__stream__)
        return 1;
    if (fseek(__stream__, 0, SEEK_SET))
        return 2; // 4 bytes

    if (fwrite(&magic, sizeof(uint32_t), 1, __stream__) != 1)
        return 3; // 6 bytes
    if (fwrite(&two, sizeof(uint16_t), 1, __stream__) != 1)
        return 5; // 8 bytes
    if (fwrite(p, sizeof(uint8_t), 58, __stream__) != 58)
        return 6; // 64 bytes

    // TODO: Do we need a summary for GTF file? 
    file->hdr->summary_offset = ftell(__stream__); 
    // record the summary offset at 64 bytes
    if (fwrite(p, sizeof(uint8_t), 32, __stream__) != 32)
        return 10;                            
    // SummaryOffset + 32 bytes empty = 96 bytes

    file->hdr->cl_offset = ftell(__stream__); 

    // record the chromosome list offset at 96 bytes
    if (BGTFWriteChromList(__stream__, file))
        return 11;
    if (writeAtPos(&(file->hdr->cl_offset), sizeof(uint64_t), 1, 0x8, __stream__)) return 7;

    // record the number of records at 0x40 (64 bytes)
    writeAtPos(&(file->n_records), sizeof(size_t), 1, 0x40, __stream__);
    // record the number of chromosomes at 0x48 (64 bytes)
    writeAtPos(&(file->n_chroms), sizeof(size_t), 1, 0x48, __stream__);

    // Save space for the number of blocks
    uint64_t offset = ftell(__stream__);
    offset = ROUNDUP(ftell(__stream__), 0x200) - offset;
    if (fwrite(p, sizeof(uint8_t), offset, __stream__) != offset)
        return 9;
    file->hdr->record_offset = ftell(__stream__);

    // record the number of records at 0x10 (16 bytes)
    if (writeAtPos(&(file->hdr->record_offset), sizeof(uint64_t), 1, 0x10, __stream__)) return 8;

    // TODO: save the on-disk chromosome offset
    free(p);
    return 0;
}


uint8_t BGTFRWRecords(BGTF *file, uint8_t mode)
{
    FILE *__stream__ = file->write->fp;
    BGTFRecord_t *record = NULL;
    fseek(__stream__, file->hdr->record_offset, SEEK_SET);
    int blocksz;
    int block_length; 
    if (mode) {
        if (file->records) {
            warn("Seems that the record ArrayList has been initialized");
        }
        file->records = ArrayListCreate(0xFF);
        for (size_t i = 0; i < file->n_records; ++i) 
        {
            if (!(fread(&block_length, 4, 1, __stream__) == 1)) {
                fatalf("Read failed at %s line %s", __FILE__, __LINE__);
            }
            if (!(fread(file->write->buffer->compressed_block, block_length, 1, __stream__) == 1)) {
                fatalf("Read failed at %s line %s", __FILE__, __LINE__);
            }
            blocksz = inflate_block(file->write->buffer, block_length);

            ArrayListPush(file->records, BGTFDeserializeRecord(file->write->buffer->uncompressed_block));
        }
    } else {
        ArrayListForEach(file->records, record, {
            BGTFSerializeRecord(file->write->buffer->uncompressed_block, &blocksz, record, NULL);
            block_length = deflate_block(file->write->buffer, blocksz);
            if (!(fwrite(&block_length, 4, 1, __stream__) == 1)) {
                fatalf("Write failed at %s line %s", __FILE__, __LINE__);
            }
            if (!(fwrite(file->write->buffer->compressed_block, block_length, 1, __stream__) == 1)) {
                fatalf("Write failed at %s line %s", __FILE__, __LINE__);
            }
        })
    }

}

uint8_t BGTFSave(BGTF *file, char *path, char *mode)
{
    logf("Saving  BGTF file to %s", path);
    if (!file->is_write) {
        file->is_write = 1;
        file->write = BGTFRWInit(path, mode);
    }
    FILE *out = fopen(path, mode);
    BGTFWriteHdr(file);
    BGTFRWRecords(file, 0);
    BGTFRWDestroy(file->write);
    file->is_write = 0;
    fclose(out);
    return 0;
}


/********
 * Read *
 *******/

uint8_t BGTFReadChromList(FILE *__stream__, BGTF *file, size_t nitems)
{
    fseek(__stream__, file->hdr->cl_offset, SEEK_SET);
    HashTable *chrom_table = file->index->chrom_idx;
    char chrom[BGTF_IO_BUFFER_SIZE];
    void *value;
    uint16_t chrom_len;
    for (size_t i = 0;  i < nitems; ++i) {
        if (!fread(&chrom_len, sizeof(uint16_t), 1, __stream__) == 1) return 1;
        if (!(fread(chrom, sizeof(uint8_t), chrom_len, __stream__) == chrom_len)) return 2;
        // We don't actually need to read chromlist now. Just testing.
        // warnf("%s", chrom);
    }
    return 0;
}



uint8_t BGTFReadHdr(FILE *__stream__, BGTF *file)
{
    int errno;
    uint32_t magic;
    uint16_t two;
    void *p = calloc(58, sizeof(uint8_t));
    if (!__stream__)
        return 1;
    if (fseek(__stream__, 0, SEEK_SET))
        return 2; // 4 bytes
    if (fread(&magic, sizeof(uint32_t), 1, __stream__) != 1)
        return 3; // 8 bytes
    if (magic != BGTF_MAGIC) {
        return 4;
    }
    fseek(__stream__, 0x8, SEEK_SET);
    if (fread(&(file->hdr->cl_offset), sizeof(uint64_t), 1, __stream__)) // 0x8
    fseek(__stream__, 0x10, SEEK_SET);
    if (fread(&(file->hdr->record_offset), sizeof(uint64_t), 1, __stream__)) // 0x10

    fseek(__stream__, 0x48, SEEK_SET);
    size_t nitems; // n_records
    if (fread(&nitems, sizeof(size_t), 1, __stream__) != 1) return 5;
    if (errno = (BGTFReadChromList(__stream__, file, nitems))) {
        fatalf("BGTFReadChromList error %d", errno);
    };
    fseek(__stream__, 0x40, SEEK_SET);
    if (fread(&nitems, sizeof(size_t), 1, __stream__) != 1) return 6;
    file->n_records = nitems;
    fseek(__stream__, file->hdr->record_offset, SEEK_SET);
    free(p);
    return 0;
}

uint8_t BGTFLoad(BGTF *file, char *path, char *attr_name)
{
    logf("Loading BGTF file %s", path);
    int errno;
    if (!file->is_write) {
        file->is_write = 1;
    }
    file->write = BGTFRWInit(path, "r");
    if (errno = BGTFReadHdr(file->write->fp, file)) {
        fatalf("BGTFLoad error %d", errno);
    };
    logf("Total record %d", file->n_records);
    BGTFRWRecords(file, 1);  // Read binary records into memory
    BGTFBuildIndex(file, attr_name);    // Build in-memory index
    log("Loading finished");
}

HashTable *readBarcodeFile(char *path)
{
    FILE *fp;
    char buf[BGTF_IO_BUFFER_SIZE];
    logf("Read barcode file %s", path);
    if (!(fp = fopen(path, "r"))) {
        fatalf("Read barcode file %s failed", path);
    }
    uint8_t l;
    
    HashTable *barcode_table = HashTableCreate(0x100);
    HashTableSetHashFunction(barcode_table, HashTableStringHashFunction);
    HashTableSetKeyComparisonFunction(barcode_table, strcmp);
    ArrayList *sort_arr = ArrayListCreate(0x2);
    while (fgets(buf, BGTF_IO_BUFFER_SIZE, fp) != NULL) {
        l = strlen(buf);
        buf[l-1] = '\0';
        char *b = malloc(l);
        memcpy(b, buf, l);
        ArrayListPush(sort_arr, b);
    }
    char *k;
    ArrayListSort(sort_arr, strcmp);
    size_t idx = 1;
    ArrayListForEach(sort_arr, k, {
        HashTablePut(barcode_table, k, idx);
        idx++;
    })
    fclose(fp);
    return barcode_table;
}

HashTable *bgtf_listdir(char *path, char *prefix, char *suffix)
{
    size_t idx = 1;
    DIR *dir;
    char *fname;
    int prefix_l = 0, suffix_l = 0;
    if (prefix) prefix_l = strlen(prefix);
    if (suffix) suffix_l = strlen(suffix);
    struct dirent *ent;
    if ((dir = opendir(path)) != NULL) {
        HashTable *dir_table = HashTableCreate(0x100);
        HashTableSetHashFunction(dir_table, HashTableStringHashFunction);
        HashTableSetKeyComparisonFunction(dir_table, strcmp);
        ArrayList *sort_arr = ArrayListCreate(0x2);
        
        while ((ent = readdir(dir)) != NULL) {
            fname = ent->d_name;
            if (ent->d_type == DT_DIR) continue;
            if ((prefix) && (strncmp(fname, prefix, prefix_l) != 0)) continue;
            if ((suffix) && (strncmp(fname + strlen(fname) - suffix_l, suffix, suffix_l) != 0)) continue;
            ArrayListPush(sort_arr, bgtf_str_concat(path, fname));
        }
        closedir(dir);
        ArrayListSort(sort_arr, strcmp);
        ArrayListForEach(sort_arr, fname, {
            logf("%s", fname);
            HashTablePut(dir_table, fname, idx);
            idx++;
        })
        ArrayListDestroy(sort_arr);
        return dir_table;
    } else {
        fatalf("The directory cannot be opened: %s", path);
    }
}
void destroyBarcodeTable(HashTable *barcode_table)
{
    HashTableSetDeallocationFunctions(barcode_table, free, NULL);
    HashTableDestroy(barcode_table);
}

int view_main(int argc, char *argv[])
{
    int c;
    char *file_path = argv[1];
    BGTF *rfp = BGTFInit();
    BGTFLoad(rfp, file_path, NULL);
    BGTFflushRecord(stdout, rfp, argv[2]);
    BGTFClose(rfp);
    return 0;
}