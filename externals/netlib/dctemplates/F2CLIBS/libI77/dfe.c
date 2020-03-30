#include "f2c.h"
#include "fio.h"
#include "fmt.h"

y_rsk(Void)
{
  if(f__curunit->uend || f__curunit->url <= f__recpos
      || f__curunit->url == 1) return 0;
  do
  {
    getc(f__cf);
  }
  while(++f__recpos < f__curunit->url);
  return 0;
}
y_getc(Void)
{
  int ch;
  if(f__curunit->uend) return(-1);
  if((ch = getc(f__cf)) != EOF)
  {
    f__recpos++;
    if(f__curunit->url >= f__recpos ||
        f__curunit->url == 1)
      return(ch);
    else  return(' ');
  }
  if(feof(f__cf))
  {
    f__curunit->uend = 1;
    errno = 0;
    return(-1);
  }
  err(f__elist->cierr, errno, "readingd");
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
  return 0;
#endif
}
#ifdef KR_headers
y_putc(c)
#else
y_putc(int c)
#endif
{
  f__recpos++;
  if(f__recpos <= f__curunit->url || f__curunit->url == 1)
    putc(c, f__cf);
  else
    err(f__elist->cierr, 110, "dout");
  return(0);
}
y_rev(Void)
{
  /*what about work done?*/
  if(f__curunit->url == 1 || f__recpos == f__curunit->url)
    return(0);
  while(f__recpos < f__curunit->url)
    (*f__putn)(' ');
  f__recpos = 0;
  return(0);
}
y_err(Void)
{
  err(f__elist->cierr, 110, "dfe");
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
  return 0;
#endif
}

y_newrec(Void)
{
  if(f__curunit->url == 1 || f__recpos == f__curunit->url)
  {
    f__hiwater = f__recpos = f__cursor = 0;
    return(1);
  }
  if(f__hiwater > f__recpos)
    f__recpos = f__hiwater;
  y_rev();
  f__hiwater = f__cursor = 0;
  return(1);
}

#ifdef KR_headers
c_dfe(a) cilist *a;
#else
c_dfe(cilist *a)
#endif
{
  f__sequential = 0;
  f__formatted = f__external = 1;
  f__elist = a;
  f__cursor = f__scale = f__recpos = 0;
  if(a->ciunit > MXUNIT || a->ciunit < 0)
    err(a->cierr, 101, "startchk");
  f__curunit = &f__units[a->ciunit];
  if(f__curunit->ufd == NULL && fk_open(DIR, FMT, a->ciunit))
    err(a->cierr, 104, "dfe");
  f__cf = f__curunit->ufd;
  if(!f__curunit->ufmt) err(a->cierr, 102, "dfe")
    if(!f__curunit->useek) err(a->cierr, 104, "dfe")
      f__fmtbuf = a->cifmt;
  (void) fseek(f__cf, (long)f__curunit->url * (a->cirec - 1), SEEK_SET);
  f__curunit->uend = 0;
  return(0);
}
#ifdef KR_headers
integer s_rdfe(a) cilist *a;
#else
integer s_rdfe(cilist *a)
#endif
{
  int n;
  if(!f__init) f_init();
  if(n = c_dfe(a))return(n);
  f__reading = 1;
  if(f__curunit->uwrt && f__nowreading(f__curunit))
    err(a->cierr, errno, "read start");
  f__getn = y_getc;
  f__doed = rd_ed;
  f__doned = rd_ned;
  f__dorevert = f__donewrec = y_err;
  f__doend = y_rsk;
  if(pars_f(f__fmtbuf) < 0)
    err(a->cierr, 100, "read start");
  fmt_bg();
  return(0);
}
#ifdef KR_headers
integer s_wdfe(a) cilist *a;
#else
integer s_wdfe(cilist *a)
#endif
{
  int n;
  if(!f__init) f_init();
  if(n = c_dfe(a)) return(n);
  f__reading = 0;
  if(f__curunit->uwrt != 1 && f__nowwriting(f__curunit))
    err(a->cierr, errno, "startwrt");
  f__putn = y_putc;
  f__doed = w_ed;
  f__doned = w_ned;
  f__dorevert = y_err;
  f__donewrec = y_newrec;
  f__doend = y_rev;
  if(pars_f(f__fmtbuf) < 0)
    err(a->cierr, 100, "startwrt");
  fmt_bg();
  return(0);
}
integer e_rdfe(Void)
{
  (void) en_fio();
  return(0);
}
integer e_wdfe(Void)
{
  (void) en_fio();
  return(0);
}
