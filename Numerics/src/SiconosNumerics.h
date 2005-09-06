#ifndef SICONOSNUMERICS_H
#define SICONOSNUMERICS_H
/*!\file SiconosNumerics.h
 *   \author Vincent Acary
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#if  defined(RIOS) && !defined(CLAPACK)
#define F77NAME(x) x
#else
#define F77NAME(x) x##_
#endif

#include "blaslapack.h"


#include "solverpack.h"

//#include "odepack.h"

#endif // SICONOSNUMERICS_H

