/*
  Dummy MA57 interface to be used when MA57 isn't available.
  10 Sept 08: First version.
*/

#include "lbl.h"
#include "ma57.h"

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "MA57_Initialize"
  Ma57_Data *Ma57_Initialize( int nz, int n, FILE *logfile ) {
      return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Analyze"
  int Ma57_Analyze( Ma57_Data *ma57 ) {
      return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Factorize"
  int Ma57_Factorize( Ma57_Data *ma57, double A[] ) {
      return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Solve"
  int Ma57_Solve( Ma57_Data *ma57, double x[] ) {
      return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Refine"
  int Ma57_Refine( Ma57_Data *ma57, double x[], double rhs[],
                   double A[], int maxitref, int job ) {
      return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Finalize"
  void Ma57_Finalize( Ma57_Data *ma57 ) {
      return;
  }

  /* ================================================================= */

#ifdef __cplusplus
}              /* Closing brace for  extern "C"  block */
#endif
