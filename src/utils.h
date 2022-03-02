#ifndef BGTF_UTILS_H
#define BGTF_UTILS_H

#include <errno.h>
#include <setjmp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include <dirent.h>
#include <signal.h>
#include <sys/ucontext.h>
#include <sys/param.h>
#include <execinfo.h>

#define MIN(a,b) (a) < (b) ? (a) : (b)
#define MAX(a,b) (a) < (b) ? (b) : (a)
#define ABS(a) (a) < 0 ? -(a) : (a)
#define MAX_THREADS  sysconf(_SC_NPROCESSORS_ONLN)

/* Optimization hints */
#if defined __GNUC__ || defined __llvm__
#define likely(x) __builtin_expect(!!(x), 1) 
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define likely(x) (x)
#define unlikely(x) (x)
#endif 

/* Defines a unique identifier of type const void*. */
#define concat(x,y) x##y

/* Automatic function overload */
/**
 * @example: int add1(int a) {return a}
 *           int add2(int a, int b) {return a+b}
 */
#define VA_NARGS_EXPAND(x) x
#define VA_NARGS(...) VA_NARGS_EXPAND(__VA_NARGS(0, ##__VA_ARGS__, \
10,9,9,7,6,5,4,3,2,1,0))
#define __VA_NARGS(_0,_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,N) N
#define VA_ARGS_FUNC(func, ...) VA_NARGS_EXPAND(concat(func, VA_NARGS(__VA_ARGS__))(__VA_ARGS__))

#define void_uid(name) \
    static const int concat(name, ___) = 0;\
    const void *name = & concat(name, ___);

#define ROUNDDOWN(a, n)               \
    ({                                \
        uint32_t __a = (uint32_t)(a); \
        (typeof(a))(__a - __a % (n)); \
    })
// Round up to the nearest multiple of n
#define ROUNDUP(a, n)                                         \
    ({                                                        \
        uint32_t __n = (uint32_t)(n);                         \
        (typeof(a))(ROUNDDOWN((uint32_t)(a) + __n - 1, __n)); \
    })

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))


/**
 * @brief Takes a pointer to a member variable and computes pointer to the 
 * structure that contains it.
 * @param ptr ptr to the member
 * @param type type of the structure
 * @param member member name
 */
#define parent_struct(ptr, type, member) \
    (ptr ? ((type*) (((char*) ptr) - offsetof(type, member))) : NULL)

#define warn(s)  do { \
    fprintf(stdout, "\033[1;33mWarning:\033[0m\t%s\n", s); \
} while (0);

#define warnf(fmt, args...) { \
    fprintf(stdout, "\033[1;33mWarning:\033[0m\t"); \
    fprintf(stdout, fmt, ## args); \
    fprintf(stdout, "\n"); \
} while (0);

#define fatal(s) do { \
    fprintf(stderr, "\033[1;31mError\033[0m:\t%s\t", s); \
    fprintf(stderr, "at (%s:%d)", __FILE__, __LINE__); \
    fprintf(stderr, "\n"); \
    exit(1); \
} while (0);


#define fatalf(fmt, args...) { \
    fprintf(stderr, "\033[1;31mError:\033[0m\t"); \
    fprintf(stderr, fmt, ## args); \
    fprintf(stderr, "\tat (%s:%d)", __FILE__, __LINE__); \
    fprintf(stderr, "\n"); \
    exit(1); \
} while (0);

#define log(s) do { \
    struct tm *t =localtime(&(time_t){ time(NULL) }); \
    char *ptr = asctime(t); \
    ptr[strlen(ptr)-1] = '\0'; \
    fprintf(stdout, "%s:\t%s\n", ptr, s); \
} while (0);


#define logf(fmt, args...) { \
    struct tm *t =localtime(&(time_t){ time(NULL) }); \
    char *ptr = asctime(t); \
    ptr[strlen(ptr)-1] = '\0'; \
    fprintf(stdout, "%s:\t", ptr); \
    fprintf(stdout, fmt, ## args); \
    fprintf(stdout, "\n"); \
} while (0);

/* File operations */
uint8_t is_binary(const void *buf, const size_t buf_len);

/* String operations */
char *bgtf_str_cpy(char *src);
char *bgtf_str_concat(char *s1, char *s2);
uint8_t bgtf_str_startswith(char *s, char *prefix);
uint8_t bgtf_str_endswith(char *s, char *suffix);

/* Hash operations */
/**
 * @brief __ac_X31_hash_string
 * 
 * @param key 
 * @return uint32_t 
 */
static inline uint32_t __ac_X31_hash_string(const void *key) {
    char *s = (char *)key;
	uint32_t h = (uint32_t)*s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + (uint32_t)*s;
	return h;
}

/**
 * @brief __ac_X15_hash_string
 * 
 * @param key 
 * @return uint15_t 
 */
static inline uint16_t __ac_X15_hash_string(const void *key) {
    char *s = (char *)key;
	uint16_t h = (uint16_t)*s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + (uint16_t)*s;
	return h;
}

#ifdef BGTF_DEBUG
static inline void crash_handler(int sig_num, siginfo_t *info, ucontext_t *ucontext)
{ /* when we crash, lets print out something usefull */
    void *array[50];
    void *caller_address;
    char **messages;
    int size, i;
    /* Get the address at the time the signal was raised */
    #if defined(__i386__) // gcc specific
    caller_address = (void *) ucontext->uc_mcontext->__ss.__eip; // EIP: x86 specific
    #elif defined(__x86_64__) // gcc specific
    caller_address = (void *) (ucontext->uc_mcontext->__ss.__rip); // RIP: x86_64 specific
    #else
    #error Unsupported architecture. // TODO: Add support for other arch.
    #endif

    switch (sig_num)
    {
#ifdef SIGSEGV
    case SIGSEGV:
        puts("This program has caused a \033[1;31mSegmentation fault\033[0m.");
        break;
#endif /* SIGSEGV */
#ifdef SIGFPE
    case SIGFPE:
        puts("This program has caused a \033[1;31mFloating point exception\033[0m.");
        break;
#endif /* SIGFPE */
#ifdef SIGILL
    case SIGILL:
        puts("This program has attempted an \033[1;31mIllegal Instruction\033[0m");
        break;
#endif /* SIGILL */
#ifdef SIGPIPE
    case SIGPIPE:
        puts("This program tried to write to a broken pipe");
        break;
#endif /* SIGPIPE */
#ifdef SIGBUS
    case SIGBUS:
        puts("This program had a \033[1;31mbus error\033[0m");
        break;
#endif /* SIGBUS */
    }

    size = backtrace(array, 5);
    array[1] = caller_address;
    messages = backtrace_symbols(array, size);

    /* skip first stack frame (points here) */
    for (i = 1; i < size && messages != NULL; ++i)
    {
    logf("[bt]: (%d) %s", i, messages[i]);
    }
    fatalf("exiting");
}
#endif

#endif