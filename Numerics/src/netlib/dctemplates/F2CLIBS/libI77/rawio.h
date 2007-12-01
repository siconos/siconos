#ifdef KR_headers
extern FILE *fdopen();
#else
#ifdef MSDOS
#include "io.h"
#define close _close
#define creat _creat
#define open _open
#define read _read
#define write _write
#endif
#ifdef __cplusplus
extern "C" {
#endif
#ifndef MSDOS
#ifdef OPEN_DECL
  extern int creat(const char*, int), open(const char*, int);
#endif
  extern int close(int);
  extern int read(int, void*, size_t), write(int, void*, size_t);
  extern int unlink(const char*);
#ifndef _POSIX_SOURCE
#ifndef NON_UNIX_STDIO
  /*
   * This type declaration may not be consistent with the one in <stdio.h>.
   * extern FILE *fdopen(int, const char*);
   */
  extern FILE *fdopen(int, char*);
#endif
#endif
#endif

  extern char *mktemp(char*);

#ifdef __cplusplus
}
#endif
#endif

#include "fcntl.h"

#ifndef O_WRONLY
#define O_RDONLY 0
#define O_WRONLY 1
#endif
