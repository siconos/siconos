#ifndef LICENSE_H
#define LICENSE_H

#include "Types.h"

FUN_DECL(Int) License_SetString(const Char *str);
FUN_DECL(Int) License_SetMagic(Unsigned Int magic);

FUN_DECL(Int) License_GetAlgorithms(Char *algs, Int len);
FUN_DECL(Int) License_GetArchitectures(Char *archs, Int len);
FUN_DECL(Int) License_GetExpiration(Int *day, Int *month, Int *year);
FUN_DECL(Int) License_GetMaxVersion(Int *ver);
FUN_DECL(Int) License_GetMaxBuild(Int *day, Int *month, Int *year);
FUN_DECL(Int) License_GetRestrictions(Int *n, Int *nnz);

FUN_DECL(Int) License_GetDemoRestrictions(Int *n, Int *nnz);

FUN_DECL(License_Termination) License_GetTermination(Void);
FUN_DECL(License_Termination) License_GetDemoTermination(Void);

#endif

