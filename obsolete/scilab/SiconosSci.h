#include "stack-c.h"
#include <string.h>
#include <stdio.h>
#include "machine.h"
//#include "sun/link.h"
//#include "siconos.h"

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif


/***************************************************
 * Declarations for siconos gateway                *
 ***************************************************/
typedef int (*gate_function)(char *);

extern int C2F(SiconosGateway)();

extern int sci_gateway(char *name, gate_function f);
