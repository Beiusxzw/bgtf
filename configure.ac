# Configure script for libbgtf
#
# Author: 
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

dnl Process this file with autoconf to produce a configure script
AC_INIT([LIBBGTF], [0.0.1], [xueziweisz@gmail.com])
AC_PREREQ(2.63)  dnl This version introduced 4-argument AC_CHECK_HEADER
AC_CONFIG_SRCDIR(./src/bgtf.c)
AC_CONFIG_HEADERS(./src/config.h)

dnl Copyright notice to be copied into the generated configure script
AC_COPYRIGHT([
This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

dnl Variant of AC_MSG_ERROR that ensures subsequent make(1) invocations fail
dnl until the configuration error is resolved and configure is run again.
AC_DEFUN([MSG_ERROR],
  [cat > config.mk <<'EOF'
ifneq ($(MAKECMDGOALS),distclean)
$(error Resolve configure error first)
endif
EOF
   AC_MSG_ERROR([$1], [$2])])

AC_PROG_CC
AC_PROG_RANLIB

dnl Avoid chicken-and-egg problem where pkg-config supplies the
dnl PKG_PROG_PKG_CONFIG macro, but we want to use it to check
dnl for pkg-config...
m4_ifdef([PKG_PROG_PKG_CONFIG], [PKG_PROG_PKG_CONFIG], [PKG_CONFIG=""])

need_crypto=no
pc_requires=
static_LDFLAGS=$LDFLAGS
static_LIBS='-lpthread -lz -lm'
private_LIBS=$LDFLAGS



AC_ARG_ENABLE([bz2],
  [AS_HELP_STRING([--disable-bz2],
                  [omit support for BZ2-compressed CRAM files])],
  [], [enable_bz2=yes])

AC_ARG_ENABLE([gcs],
  [AS_HELP_STRING([--enable-gcs],
                  [support Google Cloud Storage URLs])],
  [], [enable_gcs=check])

AC_SYS_LARGEFILE

AC_ARG_ENABLE([libcurl],
  [AS_HELP_STRING([--enable-libcurl],
                  [enable libcurl-based support for http/https/etc URLs])],
  [], [enable_libcurl=check])

AC_ARG_ENABLE([lzma],
  [AS_HELP_STRING([--disable-lzma],
                  [omit support for LZMA-compressed CRAM files])],
  [], [enable_lzma=yes])

AC_ARG_ENABLE([plugins],
  [AS_HELP_STRING([--enable-plugins],
                  [enable separately-compiled plugins for file access])],
  [], [enable_plugins=no])
AC_SUBST(enable_plugins)


AC_ARG_WITH([libdeflate],
  [AS_HELP_STRING([--with-libdeflate],
                  [use libdeflate for faster crc and deflate algorithms])],
  [], [with_libdeflate=check])


basic_host=${host_alias:-unknown-`uname -s`}
AC_MSG_CHECKING([shared library type for $basic_host])
case $basic_host in
  *-cygwin* | *-CYGWIN*)
    host_result="Cygwin DLL"
    PLATFORM=CYGWIN
    PLUGIN_EXT=.cygdll
    ;;
  *-darwin* | *-Darwin*)
    host_result="Darwin dylib"
    PLATFORM=Darwin
    PLUGIN_EXT=.bundle
    ;;
  *-msys* | *-MSYS* | *-mingw* | *-MINGW*)
    host_result="MSYS dll"
    PLATFORM=MSYS
    PLUGIN_EXT=.dll
    # This also sets __USE_MINGW_ANSI_STDIO which in turn makes PRId64,
    # %lld and %z printf formats work.  It also enforces the snprintf to
    # be C99 compliant so it returns the correct values (in kstring.c).
    CPPFLAGS="$CPPFLAGS -D_XOPEN_SOURCE=600"
    ;;
  *)
    host_result="plain .so"
    PLATFORM=default
    PLUGIN_EXT=.so
    ;;
esac
AC_MSG_RESULT([$host_result])
AC_SUBST([PLATFORM])


dnl FIXME This pulls in dozens of standard header checks
AC_FUNC_MMAP
AC_CHECK_FUNCS([gmtime_r fsync drand48 srand48_deterministic])

# Darwin has a dubious fdatasync() symbol, but no declaration in <unistd.h>
AC_CHECK_DECL([fdatasync(int)], [AC_CHECK_FUNCS(fdatasync)])



AC_SEARCH_LIBS([log], [m], [],
  [MSG_ERROR([log() not found])])

zlib_devel=ok
dnl Set a trivial non-empty INCLUDES to avoid excess default includes tests
AC_CHECK_HEADER([zlib.h], [], [zlib_devel=missing], [;])
AC_CHECK_LIB(z, inflate,  [], [zlib_devel=missing])

if test $zlib_devel != ok; then
  MSG_ERROR([zlib development files not found

The dependencies compression routines from the zlib library <http://zlib.net> and
requires zlib development files to be installed on the build
machine; you may need to ensure a package such as zlib1g-dev (on Debian or
Ubuntu Linux) or zlib-devel (on RPM-based Linux distributions or Cygwin)
is installed.

FAILED.  This error must be resolved in order to build successfully.])
fi

AC_OUTPUT

