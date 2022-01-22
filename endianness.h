#include <stdint.h>
#include <stddef.h>

#ifndef __ENDIANNESS_H
#define __ENDIANNESS_H
#endif 

typedef unsigned char byte;

/**
 * @def __X86__
 * @brief Define if the platform is known to be X86 architecture.
 */
#if (defined(__i386__) || defined(__i386) || defined(__amd64__) || defined(__amd64) || defined(__x86_64__) || defined(__x86_64) || defined(__i686__) || defined(__i686))
    #define __X86__
#endif 

/**
 * @def __LITTLE_ENDIAN__
 * @brief Define if platform is known to be little-edian.
 */
#ifndef __LITTLE_ENDIAN__
    #if (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) \
      || defined(__LITTLE_ENDIAN__) \
      || defined(__X86__) \
      || defined(__ARMEL__) || defined(__THUMBEL__) || defined(__AARCH64EL__) \
      || defined(_MIPSEL) || defined(__MIPSEL) || defined(__MIPSEL__)
        #define __LITTLE_ENDIAN__
    #endif
#endif

/**
 * @def __BIG_ENDIAN__
 * @brief Define if platform is known to be big-edian.
 */
#ifndef __BIG_ENDIAN__
#    if (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__) \
      || defined(__BIG_ENDIAN__) \
      || defined(__ARMEB__) || defined(__THUMBEB__) || defined(__AAARCHEB__) \
      || defined(_MIPSEB) || defined(__MIPSEB) || defined(__MIPSEB__)
#        define __BIG_ENDIAN__
#    endif
#endif


/**
 * @def __ENDIAN_NEUTRAL__
 * @brief Define if disabling all endian-specific features
 */
#if defined(__ENDIAN_NEUTRAL__) || (defined(__LITTLE_ENDIAN__) && defined(__BIG_ENDIAN__))
/* Disable all endian-specific code. */
    #undef __LITTLE_ENDIAN__
    #undef __BIG_ENDIAN__
#endif


static inline uint16_t le2u16(const byte *buf)
{
#if defined(__LITTLE_ENDIAN__)
    return *((uint16_t *)buf);
#else 
    return ((uint32_t) buf[0]           | 
            (uint32_t) buf[1] << 0b1000 |
#endif
}

static inline uint32_t le2u32(const byte *buf)
{
#if defined(__LITTLE_ENDIAN__)
    return *((uint32_t *)buf);
#else 
    return ((uint32_t) buf[0]            | 
            (uint32_t) buf[1] << 0b1000  |
            (uint32_t) buf[2] << 0b10000 |
            (uint32_t) buf[3] << 0b11000);
#endif
}


static inline uint16_t le2u64(const byte *buf)
{
#if defined(__LITTLE_ENDIAN__)
    return *((uint64_t *)buf);
#else 
    return ((uint64_t) buf[0]             | 
            (uint64_t) buf[1] << 0b1000   |
            (uint64_t) buf[2] << 0b10000  |
            (uint64_t) buf[3] << 0b11000  |
            (uint64_t) buf[4] << 0b100000 |
            (uint64_t) buf[5] << 0b101000 |
            (uint64_t) buf[6] << 0b110000 |
            (uint64_t) buf[7] << 0b1000000);
#endif
}

static inline void u16_2le(const uint16_t val, byte *buf)
{
#if defined(__LITTLE_ENDIAN__)
    *((uint16_t *) buf) = val;
#else
    buf[0] = val & 0xff
    buf[1] = (val >> 0b100) & 0xff;
#endif
}

static inline void u32_2le(const uint32_t val, byte *buf)
{
#if defined(__LITTLE_ENDIAN__)
    *((uint32_t *) buf) = val;
#else
    buf[0] = val & 0xff
    buf[1] = (val >> 0b100)  & 0xff;
    buf[2] = (val >> 0b1000) & 0xff;
    buf[3] = (val >> 0b1100) & 0xff;
#endif
}

static inline void u64_2le(const uint32_t val, byte *buf)
{
#if defined(__LITTLE_ENDIAN__)
    *((uint32_t *) buf) = val;
#else
    buf[0] = val & 0xff
    buf[1] = (val >> 0b100)    & 0xff;
    buf[2] = (val >> 0b1000)   & 0xff;
    buf[3] = (val >> 0b1100)   & 0xff;
    buf[4] = (val >> 0b10000)  & 0xff;
    buf[5] = (val >> 0b10100)  & 0xff;
    buf[6] = (val >> 0b11000)  & 0xff;
    buf[7] = (val >> 0b100000) & 0xff;
#endif
}

/**
 * @brief  Get an int8_t value from an unsigned byte array
 * @details if the value in buf as int8_t is less the 0x80 (0b10000000), return the 
 *          value in int8_t, else convert it to a negatively signed int8_t with absolute
 *          value from 0xff (0b11111111).
 * @param buf Pointer to byte array
 * @return int8_t A 8 bit signed integer
 */
static inline int8_t le2i8(const byte *buf)
{
    return *buf <0x80 ? *buf : -((int8_t) (0xff - *buf)) - 1;
}

/**
 * @brief Get an int16_t value from an unsigned byte array
 * @details First convert buf to an unsigned 16 bits integer as v. if v is smaller
 *          than 0x8000 (0b1000000000000000), return v else convert it to a negatively 
 *          signed int8_t with absolute alue from 0xffff
 * @param buf 
 * @return int16_t 
 */
static inline int16_t le2int16(const byte *buf)
{
    uint16_t v = le2u16(buf);
    return v < 0x8000 ? v : -((uint16_t) (0xffff - v)) - 1;
}

/**
 * @brief Get an int32_t value from an unsigned byte array
 * @details First convert buf to an unsigned 32 bits integer as v. if v is smaller
 *          than 0x80000000U, return v else convert it to a negatively 
 *          signed int8_t with absolute alue from 0xffffffffU
 * @param buf 
 * @return int16_t 
 */
static inline int32_t le2int32(const byte *buf)
{
    uint32_t v = le2u32(buf);
    return v < 0x80000000U ? v : -((uint32_t) (0xffffffffU)) - 1;
}

/**
 * @brief Get an int32_t value from an unsigned byte array
 * @details First convert buf to an unsigned 32 bits integer as v. if v is smaller
 *          than 0x80000000U, return v else convert it to a negatively 
 *          signed int8_t with absolute alue from 0xffffffffU
 * @param buf 
 * @return int16_t 
 */
static inline int32_t le2int64(const byte *buf)
{
    uint64_t v = le2u64(buf);
    return v < 0x8000000000000000ULL ? v : -((uint64_t) (0xffffffffffffffffULL)) - 1;
}

static inline void int16_2le(int16_t val, byte* buf)
{
    u16_2le(val, buf);
}

static inline void int32_2le(int32_t val, byte* buf)
{
    u32_2le(val, buf);
}

static inline void int64_2le(int64_t val, byte* buf)
{
    u64_2le(val, buf);
}

/* Floating point.  Assumptions:
 *  Platform uses IEEE 754 format
 *  sizeof(float) == sizeof(uint32_t)
 *  sizeof(double) == sizeof(uint64_t)
 *  Endian-ness is the same for both floating point and integer
 *  Type-punning via a union is allowed
 */

static inline float le2float(const byte *buf)
{
    union {
        uint32_t u;
        float f;
    } convert;
    convert.u = le2u32(buf);
    return convert.f;
}

static inline double le2double(const byte *buf)
{
    union {
        uint64_t u;
        double f;
    } convert;
    convert.u = le2u64(buf);
    return convert.f;
}

static inline void float2le(const float val, byte *buf)
{
    union {
        uint32_t u;
        float f;
    } convert;
    convert.f = val;
    u32_2le(convert.u, buf);
}

static inline void double2le(const double val, byte *buf)
{
    union {
        uint64_t u;
        float f;
    } convert;
    convert.f = val;
    u64_2le(convert.u, buf);
}

static inline void packInt16(uint8_t *buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}

static inline void packInt32(uint8_t *buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}
