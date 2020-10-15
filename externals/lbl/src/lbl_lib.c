#include "cblas.h"
#include "lbl.h"

#define USE_27  lbl->lblsolver == 0
#define USE_57  lbl->lblsolver == 1

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Initialize"
  LBL_Data *LBL_Initialize( int nz, int n, FILE *outfile, int lblsolver ) {

      LBL_Data *lbl = LBL_Calloc( 1, sizeof( LBL_Data ) );
      lbl->lblsolver = lblsolver;

      if (USE_27) {
          lbl->ma27 = MA27_Initialize(nz, n, outfile);
          lbl->irn  = lbl->ma27->irn;
          lbl->jcn  = lbl->ma27->jcn;
      }
      else if (USE_57) {
          lbl->ma57 = Ma57_Initialize(nz, n, outfile);
          lbl->irn  = lbl->ma57->irn;
          lbl->jcn  = lbl->ma57->jcn;
      }

      return lbl;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Analyze"
  int LBL_Analyze( LBL_Data *lbl, int iflag ) {
      if (USE_27)
          return MA27_Analyze(lbl->ma27, iflag);
      else if (USE_57)
          return Ma57_Analyze(lbl->ma57);
      else
          return -1;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Factorize"
  int LBL_Factorize( LBL_Data *lbl, double *A ) {
      if (USE_27)
          return MA27_Factorize(lbl->ma27, A);
      else if (USE_57)
          return Ma57_Factorize(lbl->ma57, A);
      else
          return -1;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Solve"
  int LBL_Solve( LBL_Data *lbl, double x[] ) {
      if (USE_27)
          return MA27_Solve(lbl->ma27, x);
      else if (USE_57)
          return Ma57_Solve(lbl->ma57, x);
      else
          return -1;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Refine"
  int LBL_Refine( LBL_Data *lbl, double x[], double rhs[], double A[],
                  double tol, int maxitref, int job ) {
      if (USE_27)
          return MA27_Refine(lbl->ma27, x, rhs, A, tol, maxitref);
      else if (USE_57)
          return Ma57_Refine(lbl->ma57, x, rhs, A, maxitref, job);
      else
          return -1;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Finalize"
  void LBL_Finalize( LBL_Data *lbl ) {
      if (USE_27)
          MA27_Finalize(lbl->ma27);
      if (USE_57)
          Ma57_Finalize(lbl->ma57);
      LBL_Free( lbl );
      return;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_set_int_parm"
  void LBL_set_int_parm(LBL_Data *lbl, int parm, int val) {
      if (USE_27) {
          /* Not yet implemented */
      }
      if (USE_57)
          switch (parm) {
          case LBL_I_PIV_SELECTION:
              Ma57_set_int_parm(lbl->ma57, MA57_I_PIV_SELECTION, val);
              break;
          case LBL_I_PIV_NUMERICAL:
              Ma57_set_int_parm(lbl->ma57, MA57_I_PIV_NUMERICAL, val);
              break;
          case LBL_I_SCALING:
              Ma57_set_int_parm(lbl->ma57, MA57_I_SCALING, val);
              break;
          }
      return;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_set_real_parm"
  void LBL_set_real_parm(LBL_Data *lbl, int parm, double val) {
      if (USE_27) {
          /* Not yet implemented */
      }
      if (USE_57)
          switch (parm) {
          case LBL_D_PIV_THRESH:
              Ma57_set_real_parm(lbl->ma57, MA57_D_PIV_THRESH, val);
              break;
          case LBL_D_PIV_NONZERO:
              Ma57_set_real_parm(lbl->ma57, MA57_D_PIV_NONZERO, val);
              break;
          }
      return;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_get_int_parm"
  int LBL_get_int_parm(LBL_Data *lbl, int parm) {
      if (USE_27) {
          /* Not yet implemented */
      }
      if (USE_57)
          switch (parm) {
          case LBL_I_PIV_SELECTION:
              return Ma57_get_int_parm(lbl->ma57, MA57_I_PIV_SELECTION);
          case LBL_I_PIV_NUMERICAL:
              return Ma57_get_int_parm(lbl->ma57, MA57_I_PIV_NUMERICAL);
          case LBL_I_SCALING:
              return Ma57_get_int_parm(lbl->ma57, MA57_I_SCALING);
          }
      return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_get_real_parm"
  double LBL_get_real_parm(LBL_Data *lbl, int parm) {
      if (USE_27) {
          /* Not yet implemented */
      }
      if (USE_57)
          switch (parm) {
          case LBL_D_PIV_THRESH:
              return Ma57_get_real_parm(lbl->ma57, MA57_D_PIV_THRESH);
          case LBL_D_PIV_NONZERO:
              return Ma57_get_real_parm(lbl->ma57, MA57_D_PIV_NONZERO);
          }
      return 0.0;
  }

  /* ================================================================= */

#ifdef __cplusplus
}              /* Closing brace for  extern "C"  block */
#endif
