/*
  10 Sept 08: Process_Error_Code is now static and thus private to
              this library.
*/

#include <cblas.h>
#include "lbl.h"
#include "ma57.h"

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Process_Error_Code"
static int Process_Error_Code( Ma57_Data *ma57, int nerror ) {

    LOGMSG( " [%3d]", nerror );
    /* Take appropriate action according to error code */
    /* Inflation factors are from the MA57 spec sheet  */
    switch( nerror ) {

    case 0:
      break;

      /* Warnings */
    case 1:
      LOGMSG( " Found and ignored %d indices out of range.\n", ma57->info[2] );
      break;

    case 2:
      LOGMSG( " Found and summed %d duplicate entries.\n", ma57->info[3] );
      break;

    case 3:
      LOGMSG( " Found duplicate and out-of-range indices.\n" );
      break;

    case 4:
      LOGMSG( " Matrix has rank %d, deficient.\n", ma57->info[24] );
      break;

    case 5:
      LOGMSG( " Found %d pivot sign changes in definite matrix.\n",
              ma57->info[25] );
      break;

    case 8:
      LOGMSG( " Infinity norm of solution found to be zero.\n" );
      break;

    case 10:
      LOGMSG( " Insufficient real space. Increased LFACT to %d.\n",
              ma57->lfact );
      break;

    case 11:
      LOGMSG( " Insufficient integer space. Increased LIFACT to %d.\n",
              ma57->lifact );
      break;

      /* Errors */
    case -1:
      LOGMSG( " Value of n is out of range: %d.\n", ma57->info[1] );
      break;

    case -2:
      LOGMSG( " Value of nz is out of range: %d.\n", ma57->info[1] );
      break;

    case -3:
      LOGMSG( " Adjusted size of array FACT to %d.\n", ma57->lfact );
      break;

    case -4:
      LOGMSG( " Adjusting size of array IFACT to %d.\n", ma57->lifact );
      break;

    case -5:
      LOGMSG( " Small pivot encountered at pivot step %d. Threshold = %lf.\n",
              ma57->info[1], ma57->cntl[1] );

    case -6:
      LOGMSG( " Change in pivot sign detected at pivot step %d.\n",
              ma57->info[1] );
      break;

    case -7:
      LOGMSG( " Erroneous sizing of array FACT or IFACT.\n" );
      break;

    case -8:
      LOGMSG( " Iterative refinement failed to converge.\n" );
      break;

    case -9:
      LOGMSG( " Error in user-supplied permutation array in component %d.\n",
              ma57->info[1] );
      break;

    case -10:
      LOGMSG( " Unknown pivoting strategy: %d\n.", ma57->info[1] );
      break;

    case -11:
      LOGMSG( " Size of RHS must be %d. Received %d.\n",
              ma57->n, ma57->info[1] );
      break;

    case -12:
      LOGMSG( " Invalid value of JOB (%d).\n", ma57->info[1] );
      break;

    case -13:
      LOGMSG( " Invalid number of iterative refinement steps (%d).\n",
              ma57->info[1] );
      break;

    case -14:
      LOGMSG( " Failed to estimate condition number.\n" );
      break;

    case -15:
      LOGMSG( " LKEEP has value %d, less than minimum allowed.\n",
              ma57->info[1] );
      break;

    case -16:
      LOGMSG( " Invalid number of RHS (%d).\n", ma57->info[1] );
      break;

    case -17:
      LOGMSG( " Increasing size of LWORK to %d.\n", ma57->lwork );
      break;

    case -18:
      LOGMSG( " MeTiS library not available or not found.\n" );
      break;

    default:
      LOGMSG( " Unrecognized flag from Factorize()." );
      nerror = -30;
    }    
    return nerror;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "MA57_Initialize"
  Ma57_Data *Ma57_Initialize( int nz, int n, FILE *logfile ) {

    /* Call initialize subroutine MA57ID and set defaults */
    Ma57_Data *ma57 = (Ma57_Data *)LBL_Calloc( 1, sizeof(Ma57_Data) );
    ma57->logfile   = logfile;

    LOGMSG( " MA57 :: Initializing..." );

    ma57->n         = n;
    ma57->nz        = nz;
    ma57->fetched   = 0;
    ma57->irn       = (int *)LBL_Calloc( nz, sizeof(int) );
    ma57->jcn       = (int *)LBL_Calloc( nz, sizeof(int) );
    ma57->lkeep     = 5*n + nz + imax(n,nz) + 42 + n; // Add n to suggested val.
    ma57->keep      = (int *)LBL_Calloc( ma57->lkeep, sizeof(int) );
    ma57->iwork     = (int *)LBL_Calloc( 5*n, sizeof(int) );
    ma57->work      = NULL; // Will be initialized in Ma57_Solve()

    LOGMSG( " Calling ma57id..." );
    MA57ID( ma57->cntl, ma57->icntl ); // Initialize all parameters

    // Ensure some default parameters are appropriate
    ma57->icntl[0] = -1;  // Stream for error messages.
    ma57->icntl[1] = -1;  // Stream for warning messages.
    ma57->icntl[2] = -1;  // Stream for monitoring printing.
    ma57->icntl[3] = -1;  // Stream for printing of statistics.
    ma57->icntl[4] =  0;  // Verbosity: 0=none, 1=errors, 2=1+warnings,
                          //            3=2+monitor, 4=3+input,output
    ma57->icntl[5] =  5;  // Pivot selection strategy:
                          //  0: AMD using MC47
                          //  1: User-supplied pivot sequence
                          //  2: AMD with dense row strategy
                          //  3: MD as in MA27
                          //  4: MeTiS
                          //  5: Automatic (4 with fallback on 2)
    ma57->icntl[6] =  1;  // Numerical pivoting strategy
                          //  1: Do pivoting using value in cntl(1)
                          //  2: No pivoting; exit on sign change or zero
                          //  3: No pivoting; exit if pivot modulus < cntl(2)
                          //  4: No pivoting; alter pivots to have same sign
    ma57->icntl[7] =  1;  // Memory will be re-allocated if necessary
    ma57->icntl[15] = 0;  // No need to scale system before factorizing    
    ma57->fetched = 0;

    LOGMSG( " done\n");
    return ma57;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Analyze"
  int Ma57_Analyze( Ma57_Data *ma57 ) {

    int finished = 0, error;
    LOGMSG( " MA57 :: Analyzing..." );

    /* Unpack data structure and call MA57AD */
    while( !finished ) {
      LOGMSG( "\n         Calling ma57ad... " );

      MA57AD( &(ma57->n), &(ma57->nz), ma57->irn, ma57->jcn,
              &(ma57->lkeep), ma57->keep, ma57->iwork, ma57->icntl,
              ma57->info, ma57->rinfo );

      error = ma57->info[0];

      if( !error )
        finished = 1;
      else {
        error = Process_Error_Code( ma57, error );
        if( error != -3 && error != -4 && error != -8 && error != -14 )
          return error;
      }
    }

    // Allocate data for Factorize()
    ma57->lfact = ceil( LFACT_GROW * ma57->info[8] );
    ma57->fact = (double *)LBL_Calloc( ma57->lfact, sizeof(double) );
    ma57->lifact = ceil( LIFACT_GROW * ma57->info[9] );
    ma57->ifact = (int *)LBL_Calloc( ma57->lifact, sizeof(int) );
    LBL_Free( ma57->iwork );
    ma57->iwork = (int *)LBL_Calloc( ma57->n, sizeof(int) );
    ma57->work = NULL;

    LOGMSG( " done\n");
    return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Factorize"
  int Ma57_Factorize( Ma57_Data *ma57, double A[] ) {

    // To ensure consistency with the MA27 interface, the user must pass
    // the array A, *including* its zero-th element. Internally, A+1 is
    // passed to MA57.

    double *newFact;
    int    *newIfact;
    int     newSize, one = 1, zero = 0;
    int     finished = 0, error;
    LOGMSG( " MA57 :: Factorizing...\n" );
        
    /* Unpack data structure and call MA57BD */
    while( !finished ) {
      LOGMSG( "         Calling ma57bd... " );

      MA57BD( &(ma57->n), &(ma57->nz), A, ma57->fact, &(ma57->lfact),
              ma57->ifact, &(ma57->lifact), &(ma57->lkeep), ma57->keep,
              ma57->iwork, ma57->icntl, ma57->cntl, ma57->info, ma57->rinfo );

      error = ma57->info[0];

      if( !error ) {

        finished = 1;

      } else {

        error = Process_Error_Code( ma57, error );

        if( error == -3 || error == 10 ) {

          /* Resize real workspace */
          newSize = (error == 10) ? ma57->lfact : ma57->info[16];
          newSize = ceil( LFACT_GROW * newSize );
          LOGMSG("         Resizing real workspace for factors to %d", newSize);
          newFact = (double *)LBL_Calloc( newSize, sizeof(double) );
          MA57ED( &(ma57->n), &zero, ma57->keep, ma57->fact, &(ma57->lfact),
                  newFact, &newSize, ma57->ifact, &(ma57->lifact), NULL,
                  &zero, ma57->info );
          LBL_Free( ma57->fact );
          ma57->fact = newFact; newFact = NULL;
          ma57->lfact = newSize;

        } else if( error == -4 || error == 11 ) {

          /* Resize integer workspace */
          newSize = (error == 11) ? ma57->lifact : ma57->info[17];
          newSize = ceil( LIFACT_GROW * newSize );
          LOGMSG("         Resizing int workspace for factors to %d", newSize);
          newIfact = (int *)LBL_Calloc( newSize, sizeof(int) );
          MA57ED( &(ma57->n), &one, ma57->keep, ma57->fact, &(ma57->lfact),
                  NULL, &zero, ma57->ifact, &(ma57->lifact), newIfact,
                  &newSize, ma57->info );
          LBL_Free( ma57->ifact );
          ma57->ifact = newIfact; newIfact = NULL;
          ma57->lifact = newSize;

        } else {
          finished = 1;
          if( error < 0 ) return error;
        }
      }
    }

    LOGMSG("\n");
    LOGMSG("         %-23s: %6i   %-28s: %6i\n",
           "No. of 2x2 pivots"            , ma57->info[21],
           "Pivot step for modification"  , ma57->info[26]);
    LOGMSG("         %-28s: %6i   %-28s: %6i\n",
           "No. of negative e-vals"       , ma57->info[23],
           "No. of entries in factor"     , ma57->info[13]);
    LOGMSG("         %-28s: %6i   %-28s: %6i\n",
           "Rank of factorization"        , ma57->info[24],
           "No. of pivot sign changes"    , ma57->info[25]);
    LOGMSG("         done\n");
    return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Solve"
  int Ma57_Solve( Ma57_Data *ma57, double x[] ) {

    int one = 1, finished = 0, error;

    LOGMSG( " MA57 :: Solving..." );

    ma57->job = 1;
    ma57->lrhs = ma57->n;
    ma57->nrhs = 1;
    ma57->lwork = ma57->n * ma57->nrhs;
    if( !ma57->work ) {
      LOGMSG("\n         Sizing work array to size %d ", ma57->lwork);
      ma57->work = (double *)LBL_Calloc( ma57->lwork, sizeof(double) );
    }

    while( !finished ) {
      LOGMSG( "\n         Calling ma57cd... " );

      /* Unpack data structure and call MA57CD */
      MA57CD( &(ma57->job), &(ma57->n), ma57->fact, &(ma57->lfact),
              ma57->ifact, &(ma57->lifact), &ma57->nrhs, x, &(ma57->lrhs),
              ma57->work, &(ma57->lwork), ma57->iwork, ma57->icntl,
              ma57->info );

      error = ma57->info[0];
      if( error == -17 ) {
        LBL_Free( ma57->work );
        ma57->lwork = ceil( 1.2 * ma57->lwork );
        LOGMSG("\n         Resizing work array to size %d ", ma57->lwork);
        ma57->work = (double *)LBL_Calloc( ma57->lwork, sizeof(double) );
      } else
        finished = 1;
      if( error ) error = Process_Error_Code( ma57, error );
    }

    LOGMSG( " done\n" );
    return error;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Refine"
  int Ma57_Refine( Ma57_Data *ma57, double x[], double rhs[],
                   double A[], int maxitref, int job ) {
    
    int error;

    LOGMSG( " MA57 :: Performing iterative refinement..." );

    ma57->job = job;
    /* For values of 'job' that demand the presence of the residual,
     * it should be placed by the user in ma57->residual prior to calling
     * this function.
     */

    /* Allocate work space */
    ma57->icntl[8] = imax( 1, maxitref );  // Number of refinement iterations
    ma57->icntl[9] = 1; // Return estimates of condition number
    LBL_Free( ma57->iwork );
    ma57->iwork = (int *)LBL_Calloc( ma57->n, sizeof(int) );

    LBL_Free( ma57->work );
    ma57->lwork = ma57->n;
    if( ma57->icntl[8] > 1 ) {
      ma57->lwork += 2 * ma57->n;
      if( ma57->icntl[9] > 0 )
        ma57->lwork += 2 * ma57->n;
    }
    ma57->work = (double *)LBL_Calloc( ma57->lwork, sizeof(double) );
    ma57->residual = (double *)LBL_Calloc( ma57->n, sizeof(double) );

    /* Perform iterative refinement */
    MA57DD( &(ma57->job), &(ma57->n), &(ma57->nz), A, ma57->irn,
            ma57->jcn, ma57->fact, &(ma57->lfact), ma57->ifact,
            &(ma57->lifact), rhs, x, ma57->residual, ma57->work,
            ma57->iwork, ma57->icntl, ma57->cntl, ma57->info,
            ma57->rinfo );

    error = ma57->info[0];
    if( error ) error = Process_Error_Code(ma57, ma57->info[0]);
    LOGMSG( " done\n" );
    return error;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_Finalize"
  void Ma57_Finalize( Ma57_Data *ma57 ) {

    /* Free allocated memory */
    LOGMSG( " MA57 :: Deallocating data arrays..." );
    LBL_Free( ma57->irn );
    LBL_Free( ma57->jcn );
    LBL_Free( ma57->keep );
    LBL_Free( ma57->iwork );
    LBL_Free( ma57->fact );
    LBL_Free( ma57->ifact );
    LBL_Free( ma57->work );
    if( ma57->residual ) LBL_Free(ma57->residual);
    LOGMSG( " done.\n" );
    free(ma57);
    return;
  }

  /* ================================================================= */
    /* Currently, this routine doesn't do error checking on the
       allowable parameter values. */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_set_int_parm"
  void Ma57_set_int_parm( Ma57_Data *ma57, int parm, int val ) {

      switch (parm) {
      case MA57_I_PIV_SELECTION:
          ma57->icntl[MA57_I_PIV_SELECTION] = val;
          LOGMSG(" Set %-20s = %10d\n","I_PIV_SELECTION",val);
          break;
      case MA57_I_PIV_NUMERICAL:
          ma57->icntl[MA57_I_PIV_NUMERICAL] = val;
          LOGMSG(" Set %-20s = %10d\n","I_PIV_NUMERICAL",val);
          break;
      case MA57_I_SCALING:
          ma57->icntl[MA57_I_SCALING] = val;
          LOGMSG(" Set %-20s = %10d\n","I_SCALING",val);
          break;                 
      }
      return;
  }

  /* ================================================================= */
    /* Currently, this routine doesn't do error checking on the
       allowable parameter values. */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_set_real_parm"
  void Ma57_set_real_parm( Ma57_Data *ma57, int parm, double val ) {

      switch (parm) {
      case MA57_D_PIV_THRESH:
          ma57->cntl[MA57_D_PIV_THRESH] = val;
          LOGMSG(" Set %-20s = %10.2e\n","D_PIV_THRESH",val);
          break;
      case MA57_D_PIV_NONZERO:
          ma57->cntl[MA57_D_PIV_NONZERO] = val;
          LOGMSG(" Set %-20s = %10.2e\n","D_PIV_NONZERO",val);
          break;
      }
      return;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_get_int_parm"
  int Ma57_get_int_parm( Ma57_Data *ma57, int parm ) {

      switch (parm) {
      case MA57_I_PIV_SELECTION:
          return ma57->icntl[MA57_I_PIV_SELECTION];
      case MA57_I_PIV_NUMERICAL:
          return ma57->icntl[MA57_I_PIV_NUMERICAL];
      case MA57_I_SCALING:
          return ma57->icntl[MA57_I_SCALING];
      }
      return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Ma57_get_real_parm"
  double Ma57_get_real_parm( Ma57_Data *ma57, int parm ) {

      switch (parm) {
      case MA57_D_PIV_THRESH:
          return ma57->cntl[MA57_D_PIV_THRESH];
      case MA57_D_PIV_NONZERO:
          return ma57->cntl[MA57_D_PIV_NONZERO];
      }
      return 0.0;
  }

  /* ================================================================= */

#ifdef __cplusplus
}              /* Closing brace for  extern "C"  block */
#endif
