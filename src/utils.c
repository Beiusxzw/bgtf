#include <string.h>
#include "utils.h"



uint8_t is_binary(const void *buf, const size_t buf_len)
{
    size_t suspicious_bytes = 0;
    size_t total_bytes = buf_len > 512 ? 512 : buf_len;
    const unsigned char *buf_c = buf;
    size_t i;

    if (buf_len == 0) {
        /* Is an empty file binary? Is it text? */
        return 0;
    }

    if (buf_len >= 3 && buf_c[0] == 0xEF && buf_c[1] == 0xBB && buf_c[2] == 0xBF) {
        /* UTF-8 BOM. This isn't binary. */
        return 0;
    }

    if (buf_len >= 5 && strncmp(buf, "%PDF-", 5) == 0) {
        /* PDF. This is binary. */
        return 1;
    }

    for (i = 0; i < total_bytes; i++) {
        if (buf_c[i] == '\0') {
            /* NULL char. It's binary */
            return 1;
        } else if ((buf_c[i] < 7 || buf_c[i] > 14) && (buf_c[i] < 32 || buf_c[i] > 127)) {
            /* UTF-8 detection */
            if (buf_c[i] > 193 && buf_c[i] < 224 && i + 1 < total_bytes) {
                i++;
                if (buf_c[i] > 127 && buf_c[i] < 192) {
                    continue;
                }
            } else if (buf_c[i] > 223 && buf_c[i] < 240 && i + 2 < total_bytes) {
                i++;
                if (buf_c[i] > 127 && buf_c[i] < 192 && buf_c[i + 1] > 127 && buf_c[i + 1] < 192) {
                    i++;
                    continue;
                }
            }
            suspicious_bytes++;
            /* Disk IO is so slow that it's worthwhile to do this calculation after every suspicious byte. */
            /* This is true even on a 1.6Ghz Atom with an Intel 320 SSD. */
            /* Read at least 32 bytes before making a decision */
            if (i >= 32 && (suspicious_bytes * 100) / total_bytes > 10) {
                return 1;
            }
        }
    }
    if ((suspicious_bytes * 100) / total_bytes > 10) {
        return 1;
    }

    return 0;
}

char *bgtf_str_cpy(char *src)
{
    char *dest = malloc(strlen(src)+1);
    memset(dest, 0, strlen(src)+1);
    memcpy(dest, src, strlen(src));
    return dest;
}

char *bgtf_str_concat(char *s1, char *s2)
{
    char *dest = malloc(strlen(s1) + strlen(s2) + 1);
    memset(dest, 0, strlen(s1) + strlen(s2) + 1);
    memcpy(dest, s1, strlen(s1));
    memcpy(dest + strlen(s1), s2, strlen(s2));
    return dest;
}

uint8_t bgtf_str_startswith(char *s, char *prefix)
{
    int prefix_l = strlen(prefix);
    return (strncmp(s, prefix, prefix_l) == 0);
}

uint8_t bgtf_str_endswith(char *s, char *suffix)
{
    int suffix_l = strlen(suffix);
    return (strncmp(s + strlen(s) - suffix_l, suffix, suffix_l) == 0);
}