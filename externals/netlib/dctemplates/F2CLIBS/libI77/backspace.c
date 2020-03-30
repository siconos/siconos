#include "f2c.h"
#include "fio.h"
#ifdef KR_headers
integer f_back(a) alist *a;
#else
integer f_back(alist *a)
#endif
{
  unit *b;
  int i, n, ndec;
  long x;
  char buf[32];
  if(a->aunit >= MXUNIT || a->aunit < 0)
    err(a->aerr, 101, "backspace")
    b = &f__units[a->aunit];
  if(b->useek == 0) err(a->aerr, 106, "backspace")
    if(b->ufd == NULL)
    {
      fk_open(1, 1, a->aunit);
      return(0);
    }
  if(b->uend == 1)
  {
    b->uend = 0;
    return(0);
  }
  if(b->uwrt)
  {
    (void) t_runc(a);
    if(f__nowreading(b))
      err(a->aerr, errno, "backspace")
    }
  if(b->url > 0)
  {
    long y;
    x = ftell(b->ufd);
    y = x % b->url;
    if(y == 0) x--;
    x /= b->url;
    x *= b->url;
    (void) fseek(b->ufd, x, SEEK_SET);
    return(0);
  }

  if(b->ufmt == 0)
  {
    (void) fseek(b->ufd, -(long)sizeof(int), SEEK_CUR);
    (void) fread((char *)&n, sizeof(int), 1, b->ufd);
    (void) fseek(b->ufd, -(long)n - 2 * sizeof(int), SEEK_CUR);
    return(0);
  }
  for(ndec = 2;; ndec = 1)
  {
    long y;
    y = x = ftell(b->ufd);
    if(x < sizeof(buf)) x = 0;
    else x -= sizeof(buf);
    (void) fseek(b->ufd, x, SEEK_SET);
    n = fread(buf, 1, (int)(y - x), b->ufd);
    for(i = n - ndec; i >= 0; i--)
    {
      if(buf[i] != '\n') continue;
      (void) fseek(b->ufd, (long)(i + 1 - n), SEEK_CUR);
      return(0);
    }
    if(x == 0)
    {
      (void) fseek(b->ufd, 0L, SEEK_SET);
      return(0);
    }
    else if(n <= 0) err(a->aerr, (EOF), "backspace")
      (void) fseek(b->ufd, x, SEEK_SET);
  }
}
