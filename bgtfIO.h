#ifndef LIBBGTFIO_H
#define LIBBGTFIO_H

#include <zlib.h>
#include "bgtf.h"

/*! \file bgtfIO.h
 * These are (typically internal) IO functions, so there's generally no need for you to directly use them!
 */

#define BGTF_BLOCK_SIZE 0xff00
#define BGTF_MAX_BLOCK_SIZE 0x10000
#define BGTF_BLOCK_HEADER_LENGTH 18
#define BGTF_BLOCK_FOOTER_LENGTH 8

enum BGTF_COMPRESS_ERR {
    BGTF_ERR_ZLIB   = 1,
    BGTF_ERR_HEADER = 2,
    BGTF_ERR_IO     = 4,
    BGTF_ERR_MISUSE = 8
};


BGTF *GTFRead(char *path);

void BGTFfprint(FILE *__stream__, BGTFRecord_t *record, char *k);

uint8_t BGTFSave(BGTF *file, char *path, char *mode);
uint8_t BGTFLoad(BGTF *file, char *path);
#endif 