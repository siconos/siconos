#ifndef NON_UNIX_STDIO
#include "sys/types.h"
#include "sys/stat.h"
#endif
#include "f2c.h"
#include "fio.h"
#include "fmt.h"  /* for struct syl */
#include "rawio.h"  /* for fcntl.h, fdopen */
#ifdef NON_UNIX_STDIO
#ifdef KR_headers
extern char *malloc();
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
#endif
#endif

/*global definitions*/
unit f__units[MXUNIT];  /*unit table*/
flag f__init; /*0 on entry, 1 after initializations*/
cilist *f__elist; /*active external io list*/
flag f__reading;  /*1 if reading, 0 if writing*/
flag f__cplus, f__cblank;
char *f__fmtbuf;
flag f__external; /*1 if external io, 0 if internal */
#ifdef KR_headers
int (*f__doed)(), (*f__doned)();
int (*f__doend)(), (*f__donewrec)(), (*f__dorevert)();
int (*f__getn)(), (*f__putn)(); /*for formatted io*/
#else
int (*f__getn)(void), (*f__putn)(int); /*for formatted io*/
int (*f__doed)(struct f__syl*, char*, ftnlen), (*f__doned)(struct f__syl*);
int (*f__dorevert)(void), (*f__donewrec)(void), (*f__doend)(void);
#endif
flag f__sequential; /*1 if sequential io, 0 if direct*/
flag f__formatted;  /*1 if formatted io, 0 if unformatted*/
FILE *f__cf;  /*current file*/
unit *f__curunit; /*current unit*/
int f__recpos;  /*place in current record*/
int f__cursor, f__scale;

/*error messages*/
char *F_err[] =
{
  "error in format",        /* 100 */
  "illegal unit number",        /* 101 */
  "formatted io not allowed",     /* 102 */
  "unformatted io not allowed",     /* 103 */
  "direct io not allowed",      /* 104 */
  "sequential io not allowed",      /* 105 */
  "can't backspace file",       /* 106 */
  "null file name",       /* 107 */
  "can't stat file",        /* 108 */
  "unit not connected",       /* 109 */
  "off end of record",        /* 110 */
  "truncation failed in endfile",     /* 111 */
  "incomprehensible list input",      /* 112 */
  "out of free space",        /* 113 */
  "unit not connected",       /* 114 */
  "read unexpected character",      /* 115 */
  "bad logical input field",      /* 116 */
  "bad variable type",        /* 117 */
  "bad namelist name",        /* 118 */
  "variable not in namelist",     /* 119 */
  "no end record",        /* 120 */
  "variable count incorrect",     /* 121 */
  "subscript for scalar variable",    /* 122 */
  "invalid array section",      /* 123 */
  "substring out of bounds",      /* 124 */
  "subscript out of bounds",      /* 125 */
  "can't read file",        /* 126 */
  "can't write file",       /* 127 */
  "'new' file exists",        /* 128 */
  "can't append to file"        /* 129 */
};
#define MAXERR (sizeof(F_err)/sizeof(char *)+100)

#ifdef KR_headers
f__canseek(f) FILE *f; /*SYSDEP*/
#else
f__canseek(FILE *f) /*SYSDEP*/
#endif
{
#ifdef NON_UNIX_STDIO
  return !isatty(fileno(f));
#else
  struct stat x;

  if(fstat(fileno(f), &x) < 0)
    return(0);
#ifdef S_IFMT
  switch(x.st_mode & S_IFMT)
  {
  case S_IFDIR:
  case S_IFREG:
    if(x.st_nlink > 0)  /* !pipe */
      return(1);
    else
      return(0);
  case S_IFCHR:
    if(isatty(fileno(f)))
      return(0);
    return(1);
#ifdef S_IFBLK
  case S_IFBLK:
    return(1);
#endif
  }
#else
#ifdef S_ISDIR
  /* POSIX version */
  if(S_ISREG(x.st_mode) || S_ISDIR(x.st_mode))
  {
    if(x.st_nlink > 0)  /* !pipe */
      return(1);
    else
      return(0);
  }
  if(S_ISCHR(x.st_mode))
  {
    if(isatty(fileno(f)))
      return(0);
    return(1);
  }
  if(S_ISBLK(x.st_mode))
    return(1);
#else
  Help! How does fstat work on this system ?
#endif
#endif
  return(0);  /* who knows what it is? */
#endif
}

void
#ifdef KR_headers
f__fatal(n, s) char *s;
#else
f__fatal(int n, char *s)
#endif
{
  if(n < 100 && n >= 0) perror(s);  /*SYSDEP*/
  else if(n >= (int)MAXERR || n < -1)
  {
    fprintf(stderr, "%s: illegal error number %d\n", s, n);
  }
  else if(n == -1) fprintf(stderr, "%s: end of file\n", s);
  else
    fprintf(stderr, "%s: %s\n", s, F_err[n - 100]);
  if(f__curunit)
  {
    fprintf(stderr, "apparent state: unit %d ", f__curunit - f__units);
    fprintf(stderr, f__curunit->ufnm ? "named %s\n" : "(unnamed)\n",
            f__curunit->ufnm);
  }
  else
    fprintf(stderr, "apparent state: internal I/O\n");
  if(f__fmtbuf)
    fprintf(stderr, "last format: %s\n", f__fmtbuf);
  fprintf(stderr, "lately %s %s %s %s", f__reading ? "reading" : "writing",
          f__sequential ? "sequential" : "direct", f__formatted ? "formatted" : "unformatted",
          f__external ? "external" : "internal");
  sig_die(" IO", 1);
}
/*initialization routine*/
VOID
f_init(Void)
{
  unit *p;

  f__init = 1;
  p = &f__units[0];
  p->ufd = stderr;
  p->useek = f__canseek(stderr);
#ifdef NON_UNIX_STDIO
  setbuf(stderr, (char *)malloc(BUFSIZ));
#else
  stderr->_flag &= ~_IONBF;
#endif
  p->ufmt = 1;
  p->uwrt = 1;
  p = &f__units[5];
  p->ufd = stdin;
  p->useek = f__canseek(stdin);
  p->ufmt = 1;
  p->uwrt = 0;
  p = &f__units[6];
  p->ufd = stdout;
  p->useek = f__canseek(stdout);
  /* IOLBUF and setvbuf only in system 5+ */
#ifdef COMMENTED_OUT
  if(isatty(fileno(stdout)))
  {
    extern char _sobuf[];
    setbuf(stdout, _sobuf);
    /* setvbuf(stdout, _IOLBF, 0, 0); /* the buf arg in setvbuf? */
    p->useek = 1; /* only within a record no bigger than BUFSIZ */
  }
#endif
  p->ufmt = 1;
  p->uwrt = 1;
}
#ifdef KR_headers
f__nowreading(x) unit *x;
#else
f__nowreading(unit *x)
#endif
{
  long loc;
  extern char *f__r_mode[];
  if(!x->ufnm)
    goto cantread;
  loc = ftell(x->ufd);
  if(freopen(x->ufnm, f__r_mode[x->ufmt], x->ufd) == NULL)
  {
cantread:
    errno = 126;
    return(1);
  }
  x->uwrt = 0;
  (void) fseek(x->ufd, loc, SEEK_SET);
  return(0);
}
#ifdef KR_headers
f__nowwriting(x) unit *x;
#else
f__nowwriting(unit *x)
#endif
{
  long loc;
  int k;
  extern char *f__w_mode[];

  if(!x->ufnm)
    goto cantwrite;
  if(x->uwrt == 3)    /* just did write, rewind */
  {
#ifdef NON_UNIX_STDIO
    if(!(f__cf = x->ufd =
                   freopen(x->ufnm, f__w_mode[x->ufmt], x->ufd)))
#else
    if(close(creat(x->ufnm, 0666)))
#endif
      goto cantwrite;
  }
  else
  {
    loc = ftell(x->ufd);
#ifdef NON_UNIX_STDIO
    if(!(f__cf = x->ufd =
                   freopen(x->ufnm, f__w_mode[x->ufmt + 2], x->ufd)))
#else
    if(fclose(x->ufd) < 0
        || (k = x->uwrt == 2 ? creat(x->ufnm, 0666)
                : open(x->ufnm, O_WRONLY)) < 0
        || (f__cf = x->ufd = fdopen(k, f__w_mode[x->ufmt])) == NULL)
#endif
    {
      x->ufd = NULL;
cantwrite:
      errno = 127;
      return(1);
    }
    (void) fseek(x->ufd, loc, SEEK_SET);
  }
  x->uwrt = 1;
  return(0);
}

int
#ifdef KR_headers
err__fl(f, m, s) int f, m;
char *s;
#else
err__fl(int f, int m, char *s)
#endif
{
  if(!f)
    f__fatal(m, s);
  if(f__doend)
    (*f__doend)();
  return errno = m;
}
