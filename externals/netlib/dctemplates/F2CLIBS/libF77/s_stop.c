#include "stdio.h"
#include "f2c.h"

#ifdef KR_headers
extern void f_exit();
VOID s_stop(s, n) char *s;
ftnlen n;
#else
#undef abs
#include "stdlib.h"
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif
void f_exit(void);

int s_stop(char *s, ftnlen n)
#endif
{
int i;

if(n > 0)
{
  fprintf(stderr, "STOP ");
  for(i = 0; i < n ; ++i)
    putc(*s++, stderr);
  fprintf(stderr, " statement executed\n");
}
f_exit();
exit(0);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
return 0; /* NOT REACHED */
}
#endif
}
