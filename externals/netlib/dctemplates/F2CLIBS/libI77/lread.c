#include "f2c.h"
#include "fio.h"
#include "fmt.h"
#include "lio.h"
#include "ctype.h"
#include "fp.h"

extern char *f__fmtbuf;
#ifdef KR_headers
extern double atof();
extern char *malloc(), *realloc();
int (*f__lioproc)(), (*l_getc)(), (*l_ungetc)();
#else
#undef abs
#undef min
#undef max
#include "stdlib.h"
int (*f__lioproc)(ftnint*, char*, ftnlen, ftnint), (*l_getc)(void),
    (*l_ungetc)(int, FILE*);
#endif
int l_eof;

#define isblnk(x) (f__ltab[x+1]&B)
#define issep(x) (f__ltab[x+1]&SX)
#define isapos(x) (f__ltab[x+1]&AX)
#define isexp(x) (f__ltab[x+1]&EX)
#define issign(x) (f__ltab[x+1]&SG)
#define iswhit(x) (f__ltab[x+1]&WH)
#define SX 1
#define B 2
#define AX 4
#define EX 8
#define SG 16
#define WH 32
char f__ltab[128 + 1] = /* offset one for EOF */
{
  0,
  0, 0, AX, 0, 0, 0, 0, 0, 0, WH | B, SX | WH, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  SX | B | WH, 0, AX, 0, 0, 0, 0, AX, 0, 0, 0, SG, SX, SG, 0, SX,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, EX, EX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  AX, 0, 0, 0, EX, EX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

#ifdef ungetc
static int
#ifdef KR_headers
un_getc(x, f__cf) int x;
FILE *f__cf;
#else
un_getc(int x, FILE *f__cf)
#endif
{
  return ungetc(x, f__cf);
}
#else
#define un_getc ungetc
#ifdef KR_headers
extern int ungetc();
#endif
#endif

t_getc(Void)
{
  int ch;
  if(f__curunit->uend) return(EOF);
  if((ch = getc(f__cf)) != EOF) return(ch);
  if(feof(f__cf))
    f__curunit->uend = l_eof = 1;
  return(EOF);
}
integer e_rsle(Void)
{
  int ch;
  if(f__curunit->uend) return(0);
  while((ch = t_getc()) != '\n' && ch != EOF);
  return(0);
}

flag f__lquit;
int f__lcount, f__ltype, nml_read;
char *f__lchar;
double f__lx, f__ly;
#define ERR(x) if(n=(x)) return(n)
#define GETC(x) (x=(*l_getc)())
#define Ungetc(x,y) (*l_ungetc)(x,y)

#ifdef KR_headers
l_R(poststar) int poststar;
#else
l_R(int poststar)
#endif
{
  char s[FMAX + EXPMAXDIGS + 4];
  register int ch;
  register char *sp, *spe, *sp1;
  long e, exp;
  int havenum, havestar, se;

  if(!poststar)
  {
    if(f__lcount > 0)
      return(0);
    f__lcount = 1;
  }
  f__ltype = 0;
  exp = 0;
  havestar = 0;
retry:
  sp1 = sp = s;
  spe = sp + FMAX;
  havenum = 0;

  switch(GETC(ch))
  {
  case '-':
    *sp++ = ch;
    sp1++;
    spe++;
  case '+':
    GETC(ch);
  }
  while(ch == '0')
  {
    ++havenum;
    GETC(ch);
  }
  while(isdigit(ch))
  {
    if(sp < spe) *sp++ = ch;
    else ++exp;
    GETC(ch);
  }
  if(ch == '*' && !poststar)
  {
    if(sp == sp1 || exp || *s == '-')
    {
      errfl(f__elist->cierr, 112, "bad repetition count");
    }
    poststar = havestar = 1;
    *sp = 0;
    f__lcount = atoi(s);
    goto retry;
  }
  if(ch == '.')
  {
    GETC(ch);
    if(sp == sp1)
      while(ch == '0')
      {
        ++havenum;
        --exp;
        GETC(ch);
      }
    while(isdigit(ch))
    {
      if(sp < spe)
      {
        *sp++ = ch;
        --exp;
      }
      GETC(ch);
    }
  }
  se = 0;
  if(issign(ch))
    goto signonly;
  if(isexp(ch))
  {
    GETC(ch);
    if(issign(ch))
    {
signonly:
      if(ch == '-') se = 1;
      GETC(ch);
    }
    if(!isdigit(ch))
    {
bad:
      errfl(f__elist->cierr, 112, "exponent field");
    }

    e = ch - '0';
    while(isdigit(GETC(ch)))
    {
      e = 10 * e + ch - '0';
      if(e > EXPMAX)
        goto bad;
    }
    if(se)
      exp -= e;
    else
      exp += e;
  }
  (void) Ungetc(ch, f__cf);
  if(sp > sp1)
  {
    ++havenum;
    while(*--sp == '0')
      ++exp;
    if(exp)
      sprintf(sp + 1, "e%ld", exp);
    else
      sp[1] = 0;
    f__lx = atof(s);
  }
  else
    f__lx = 0.;
  if(havenum)
    f__ltype = TYLONG;
  else
    switch(ch)
    {
    case ',':
    case '/':
      break;
    default:
      if(havestar && (ch == ' '
                      || ch == '\t'
                      || ch == '\n'))
        break;
      if(nml_read > 1)
      {
        f__lquit = 2;
        return 0;
      }
      errfl(f__elist->cierr, 112, "invalid number");
    }
  return 0;
}

static int
#ifdef KR_headers
rd_count(ch) register int ch;
#else
rd_count(register int ch)
#endif
{
  if(ch < '0' || ch > '9')
    return 1;
  f__lcount = ch - '0';
  while(GETC(ch) >= '0' && ch <= '9')
    f__lcount = 10 * f__lcount + ch - '0';
  Ungetc(ch, f__cf);
  return f__lcount <= 0;
}

l_C(Void)
{
  int ch, nml_save;
  double lz;
  if(f__lcount > 0) return(0);
  f__ltype = 0;
  GETC(ch);
  if(ch != '(')
  {
    if(nml_read > 1 && (ch < '0' || ch > '9'))
    {
      Ungetc(ch, f__cf);
      f__lquit = 2;
      return 0;
    }
    if(rd_count(ch))
      if(!f__cf || !feof(f__cf))
        errfl(f__elist->cierr, 112, "complex format");
      else
        err(f__elist->cierr, (EOF), "lread");
    if(GETC(ch) != '*')
    {
      if(!f__cf || !feof(f__cf))
        errfl(f__elist->cierr, 112, "no star");
      else
        err(f__elist->cierr, (EOF), "lread");
    }
    if(GETC(ch) != '(')
    {
      Ungetc(ch, f__cf);
      return(0);
    }
  }
  else
    f__lcount = 1;
  while(iswhit(GETC(ch)));
  Ungetc(ch, f__cf);
  nml_save = nml_read;
  nml_read = 0;
  if(ch = l_R(1))
    return ch;
  if(!f__ltype)
    errfl(f__elist->cierr, 112, "no real part");
  lz = f__lx;
  while(iswhit(GETC(ch)));
  if(ch != ',')
  {
    (void) Ungetc(ch, f__cf);
    errfl(f__elist->cierr, 112, "no comma");
  }
  while(iswhit(GETC(ch)));
  (void) Ungetc(ch, f__cf);
  if(ch = l_R(1))
    return ch;
  if(!f__ltype)
    errfl(f__elist->cierr, 112, "no imaginary part");
  while(iswhit(GETC(ch)));
  if(ch != ')') errfl(f__elist->cierr, 112, "no )");
  f__ly = f__lx;
  f__lx = lz;
  nml_read = nml_save;
  return(0);
}
l_L(Void)
{
  int ch;
  if(f__lcount > 0) return(0);
  f__ltype = 0;
  GETC(ch);
  if(isdigit(ch))
  {
    rd_count(ch);
    if(GETC(ch) != '*')
      if(!f__cf || !feof(f__cf))
        errfl(f__elist->cierr, 112, "no star");
      else
        err(f__elist->cierr, (EOF), "lread");
    GETC(ch);
  }
  if(ch == '.') GETC(ch);
  switch(ch)
  {
  case 't':
  case 'T':
    f__lx = 1;
    break;
  case 'f':
  case 'F':
    f__lx = 0;
    break;
  default:
    if(isblnk(ch) || issep(ch) || ch == EOF)
    {
      (void) Ungetc(ch, f__cf);
      return(0);
    }
    else  errfl(f__elist->cierr, 112, "logical");
  }
  f__ltype = TYLONG;
  f__lcount = 1;
  while(!issep(GETC(ch)) && ch != EOF);
  (void) Ungetc(ch, f__cf);
  return(0);
}
#define BUFSIZE 128
l_CHAR(Void)
{
  int ch, size, i;
  char quote, *p;
  if(f__lcount > 0) return(0);
  f__ltype = 0;
  if(f__lchar != NULL) free(f__lchar);
  size = BUFSIZE;
  p = f__lchar = (char *)malloc((unsigned int)size);
  if(f__lchar == NULL)
    errfl(f__elist->cierr, 113, "no space");

  GETC(ch);
  if(isdigit(ch))
  {
    /* allow Fortran 8x-style unquoted string...  */
    /* either find a repetition count or the string */
    f__lcount = ch - '0';
    *p++ = ch;
    for(i = 1;;)
    {
      switch(GETC(ch))
      {
      case '*':
        if(f__lcount == 0)
        {
          f__lcount = 1;
          goto noquote;
        }
        p = f__lchar;
        goto have_lcount;
      case ',':
      case ' ':
      case '\t':
      case '\n':
      case '/':
        Ungetc(ch, f__cf);
      /* no break */
      case EOF:
        f__lcount = 1;
        f__ltype = TYCHAR;
        return *p = 0;
      }
      if(!isdigit(ch))
      {
        f__lcount = 1;
        goto noquote;
      }
      *p++ = ch;
      f__lcount = 10 * f__lcount + ch - '0';
      if(++i == size)
      {
        f__lchar = (char *)realloc(f__lchar,
                                   (unsigned int)(size += BUFSIZE));
        p = f__lchar + i;
      }
    }
  }
  else(void) Ungetc(ch, f__cf);
have_lcount:
  if(GETC(ch) == '\'' || ch == '"') quote = ch;
  else if(isblnk(ch) || (issep(ch) && ch != '\n') || ch == EOF)
  {
    (void) Ungetc(ch, f__cf);
    return(0);
  }
  else
  {
    /* Fortran 8x-style unquoted string */
    *p++ = ch;
    for(i = 1;;)
    {
      switch(GETC(ch))
      {
      case ',':
      case ' ':
      case '\t':
      case '\n':
      case '/':
        Ungetc(ch, f__cf);
      /* no break */
      case EOF:
        f__ltype = TYCHAR;
        return *p = 0;
      }
noquote:
      *p++ = ch;
      if(++i == size)
      {
        f__lchar = (char *)realloc(f__lchar,
                                   (unsigned int)(size += BUFSIZE));
        p = f__lchar + i;
      }
    }
  }
  f__ltype = TYCHAR;
  for(i = 0;;)
  {
    while(GETC(ch) != quote && ch != '\n'
          && ch != EOF && ++i < size) *p++ = ch;
    if(i == size)
    {
newone:
      f__lchar = (char *)realloc(f__lchar,
                                 (unsigned int)(size += BUFSIZE));
      p = f__lchar + i - 1;
      *p++ = ch;
    }
    else if(ch == EOF) return(EOF);
    else if(ch == '\n')
    {
      if(*(p - 1) != '\\') continue;
      i--;
      p--;
      if(++i < size) *p++ = ch;
      else goto newone;
    }
    else if(GETC(ch) == quote)
    {
      if(++i < size) *p++ = ch;
      else goto newone;
    }
    else
    {
      (void) Ungetc(ch, f__cf);
      *p = 0;
      return(0);
    }
  }
}
#ifdef KR_headers
c_le(a) cilist *a;
#else
c_le(cilist *a)
#endif
{
  f__fmtbuf = "list io";
  if(a->ciunit >= MXUNIT || a->ciunit < 0)
    err(a->cierr, 101, "stler");
  f__scale = f__recpos = 0;
  f__elist = a;
  f__curunit = &f__units[a->ciunit];
  if(f__curunit->ufd == NULL && fk_open(SEQ, FMT, a->ciunit))
    err(a->cierr, 102, "lio");
  f__cf = f__curunit->ufd;
  if(!f__curunit->ufmt) err(a->cierr, 103, "lio")
    return(0);
}
#ifdef KR_headers
l_read(number, ptr, len, type) ftnint *number, type;
char *ptr;
ftnlen len;
#else
l_read(ftnint *number, char *ptr, ftnlen len, ftnint type)
#endif
{
#define Ptr ((flex *)ptr)
  int i, n, ch;
  doublereal *yy;
  real *xx;
  for(i = 0; i < *number; i++)
  {
    if(f__lquit) return(0);
    if(l_eof)
      err(f__elist->ciend, EOF, "list in")
      if(f__lcount == 0)
      {
        f__ltype = 0;
        for(;;)
        {
          GETC(ch);
          switch(ch)
          {
          case EOF:
            goto loopend;
          case ' ':
          case '\t':
          case '\n':
            continue;
          case '/':
            f__lquit = 1;
            goto loopend;
          case ',':
            f__lcount = 1;
            goto loopend;
          default:
            (void) Ungetc(ch, f__cf);
            goto rddata;
          }
        }
      }
rddata:
    switch((int)type)
    {
    case TYINT1:
    case TYSHORT:
    case TYLONG:
#ifdef TYQUAD
    case TYQUAD:
#endif
    case TYREAL:
    case TYDREAL:
      ERR(l_R(0));
      break;
    case TYCOMPLEX:
    case TYDCOMPLEX:
      ERR(l_C());
      break;
    case TYLOGICAL1:
    case TYLOGICAL2:
    case TYLOGICAL:
      ERR(l_L());
      break;
    case TYCHAR:
      ERR(l_CHAR());
      break;
    }
    while(GETC(ch) == ' ' || ch == '\t');
    if(ch != ',' || f__lcount > 1)
      Ungetc(ch, f__cf);
loopend:
    if(f__lquit) return(0);
    if(f__cf)
    {
      if(feof(f__cf))
        err(f__elist->ciend, (EOF), "list in")
        else if(ferror(f__cf))
        {
          clearerr(f__cf);
          errfl(f__elist->cierr, errno, "list in");
        }
    }
    if(f__ltype == 0) goto bump;
    switch((int)type)
    {
    case TYINT1:
    case TYLOGICAL1:
      Ptr->flchar = (char)f__lx;
      break;
    case TYLOGICAL2:
    case TYSHORT:
      Ptr->flshort = (short)f__lx;
      break;
    case TYLOGICAL:
    case TYLONG:
      Ptr->flint = f__lx;
      break;
#ifdef TYQUAD
    case TYQUAD:
      Ptr->fllongint = f__lx;
      break;
#endif
    case TYREAL:
      Ptr->flreal = f__lx;
      break;
    case TYDREAL:
      Ptr->fldouble = f__lx;
      break;
    case TYCOMPLEX:
      xx = (real *)ptr;
      *xx++ = f__lx;
      *xx = f__ly;
      break;
    case TYDCOMPLEX:
      yy = (doublereal *)ptr;
      *yy++ = f__lx;
      *yy = f__ly;
      break;
    case TYCHAR:
      b_char(f__lchar, ptr, len);
      break;
    }
bump:
    if(f__lcount > 0) f__lcount--;
    ptr += len;
    if(nml_read)
      nml_read++;
  }
  return(0);
#undef Ptr
}
#ifdef KR_headers
integer s_rsle(a) cilist *a;
#else
integer s_rsle(cilist *a)
#endif
{
  int n;

  if(!f__init) f_init();
  if(n = c_le(a)) return(n);
  f__reading = 1;
  f__external = 1;
  f__formatted = 1;
  f__lioproc = l_read;
  f__lquit = 0;
  f__lcount = 0;
  l_eof = 0;
  if(f__curunit->uwrt && f__nowreading(f__curunit))
    err(a->cierr, errno, "read start");
  l_getc = t_getc;
  l_ungetc = un_getc;
  f__doend = xrd_SL;
  return(0);
}
