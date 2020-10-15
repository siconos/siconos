/*
  10 Sept 08: Process_Error_Code is now static and thus private to
              this library.
*/

#include "cblas.h"
#include "lbl.h"
#include "ma27.h"

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

  /*----------------------------------------------------------------------
    --
    -- Our very own dnrm_infty.
    --*/

  static double
  cblas_dnrm_infty( const int N, const double *X, const int incX ) {
    return fabs( X[ cblas_idamax( N, X, incX ) ] );
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "Process_Error_Code"
static int Process_Error_Code( Ma27_Data *ma27, int nerror ) {

    int newsize;

    PRINT( "    [nerror %-2d]\n", nerror );

    /* Take appropriate action according to error code */
    /* Inflation factors are from the MA27 spec sheet  */
    PRINT( "Message: " );
    switch( nerror ) {
      /* Warnings */
    case 1:
      PRINT( "Found %d indices out of range\n", ma27->info[1] );
      break;

    case 2:
      PRINT( "Found %d sign changes in supposedly definite matrix\n",
             ma27->info[1] );
      break;

    case 3:
      PRINT( "Matrix has rank %d, deficient\n", ma27->info[1] );
      break;

      /* Errors */
    case -1:
      PRINT( "Value of n is out of range: %d\n", ma27->n );
      break;

    case -2:
      PRINT( "Value of nz is out of range: %d\n", ma27->nz );
      break;

    case -3:
      ma27->liw = ceil( 1.2 * ma27->info[1] );
      PRINT( "Adjusting size of array IW to %d\n", ma27->liw );
      LBL_Free( ma27->iw ); ma27->iw = NULL;
      ma27->iw = LBL_Calloc( ma27->liw, sizeof( int ) );
      break;

    case -4:
      newsize = ceil( 1.2 * ma27->info[1] );
      PRINT( "Adjusting size of array FACTORS to %d\n", newsize );
      LBL_Free( ma27->factors );
      ma27->factors = LBL_Calloc( newsize, sizeof( double ) );
      ma27->la = newsize;
      break;

    case -5:
      PRINT( "Matrix singularity detected at pivot step %d\n",
             ma27->info[1] );
      break;

    case -6:
      PRINT( "Change in pivot sign detected at pivot step %d\n",
             ma27->info[1] );
      break;

    case -7:
      PRINT( "Value of nsteps out of range: %d\n", ma27->nsteps );
      break;

    default:
      PRINT( "Unrecognized flag from Factorize()" );
      nerror = -30;
    }    
    return nerror;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "MA27_Initialize"
  Ma27_Data *MA27_Initialize( int nz, int n, FILE *outfile ) {

    /* Call initialize subroutine MA27ID and set defaults */
    Ma27_Data *ma27 = LBL_Calloc( 1, sizeof( Ma27_Data ) );
    ma27->outfile   = outfile; // (Re)Set output file.
    PRINT( " MA27 :: Initializing..." );

    ma27->n         = n;
    ma27->nz        = nz;
    ma27->fetched   = 0;
    ma27->la        = ceil( 1.2 * nz );
    ma27->liw       = imax( ceil( 1.2 * ( 2*nz + 3*n )), LIW_MIN );
    ma27->irn       = LBL_Calloc( nz, sizeof( int ) );
    ma27->jcn       = LBL_Calloc( nz, sizeof( int ) );
    ma27->iw        = LBL_Calloc( ma27->liw, sizeof( int ) );
    ma27->ikeep     = LBL_Calloc( 3*n, sizeof( int ) );
    ma27->iw1       = LBL_Calloc( 2*n, sizeof( int ) );
    ma27->factors   = LBL_Calloc( ma27->la, sizeof( double ) );
    ma27->residual  = LBL_Calloc( ma27->n, sizeof( double ) );

    MA27ID( ma27->icntl, ma27->cntl );
    ma27->icntl[0]  = 0;  // Stream for error messages.
    ma27->icntl[1]  = 0;  // Stream for diagnotic messages.
    ma27->icntl[2]  = 0;  // Verbosity: 0=none, 1=partial, 2=full
    ma27->cntl[0]   = 1.0e-15; // Pivot-for-stability threshold

    PRINT( " done\n" );

    return ma27;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "MA27_Analyze"
  int MA27_Analyze( Ma27_Data *ma27, int iflag ) {

    int finished = 0, error;

    PRINT( " MA27 :: Analyzing..." );

    ma27->iflag = iflag;

    /* Unpack data structure and call MA27AD */
    while( !finished ) {

      MA27AD(&(ma27->n), &(ma27->nz), ma27->irn, ma27->jcn,
             ma27->iw,
             &(ma27->liw), ma27->ikeep, ma27->iw1, &(ma27->nsteps),
             &(ma27->iflag), ma27->icntl, ma27->cntl, ma27->info,
             &(ma27->ops) );

      error = ma27->info[0];

      if (ma27->icntl[2]) {
        fflush(stdout);
        fflush(stderr);
      }

      if( !error )
        finished = 1;
      else {
        error = Process_Error_Code( ma27, error );
        if( error != -3 && error != -4 ) return error;
      }
    }

    /* Adjust size of factors (if necessary) */
    if( ma27->info[4] > ma27->nz ) {
      ma27->la = ceil( 1.2 * ma27->info[4] );
      LBL_Free( ma27->factors );
      ma27->factors = LBL_Calloc( ma27->la, sizeof( double ) );
    }

    /* Adjust size of w1 (if necessary) */
    if( ma27->nsteps > 2 * ma27->n ) {
      LBL_Free( ma27->iw1 );
      ma27->iw1 = LBL_Calloc( ma27->nsteps, sizeof( int ) );
    }

    /* For now we assume the front size is maximal. */
    ma27->w = LBL_Calloc( ma27->n, sizeof( double ) );

    PRINT( " done\n" );

    return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "MA27_Factorize"
  int MA27_Factorize( Ma27_Data *ma27, double *A ) {

    int finished = 0, error;
    double u, new_u;

    PRINT( " MA27 :: Factorizing..." );

    /* Copy A into factors. */
    cblas_dcopy( ma27->nz, A, 1, ma27->factors, 1 );

    /* Unpack data structure and call MA27BD */
    while( !finished ) {

      MA27BD( &(ma27->n), &(ma27->nz), ma27->irn, ma27->jcn,
              ma27->factors, &(ma27->la), ma27->iw, &(ma27->liw),
              ma27->ikeep, &(ma27->nsteps), &(ma27->maxfrt),
              ma27->iw1, ma27->icntl, ma27->cntl, ma27->info );

      error = ma27->info[0];
      u = ma27->cntl[0];     // Threshold controling stability by pivoting

      if( !error ) {
        finished = 1;
      } else if( (error == 3 || error == -5) && (u <= PIV_MAX) ) {
        // Matrix is singular and u is less than max
        /* Adjust pivot tolerance */
        if( u == 0.0 )
          u =  1.0e-6;
        else {
          new_u = fmin( PIV_MAX, fmax( PIV_MIN, 100 * u ) );
          if( new_u == u ) {
            /* We have tried all allowed pivot tolerances */
            PRINT("MA27: Failed to factorize matrix\n" );
            PRINT("      Order = %-d, ", ma27->n );
            if( error == 3 ) PRINT("Rank =");
            else PRINT("Singularity at pivot step" );
            PRINT(" %-d\n", ma27->info[1] );
            return -2;
          }
        }
        PRINT("Adjusting pivot-for-stability threshold to %g\n", new_u);
        ma27->cntl[0] = new_u;

        //if( error == 3 ) finished = 1;  // A factorization was produced!
      }
      else {
        error = Process_Error_Code( ma27, error );
        if( error == -4 ) {
          /* Must re-initialize factors */
          cblas_dcopy( ma27->nz, A, 1, ma27->factors, 1 );
        }
        /* If error > 0, Ma27 only issued a warning; we keep on going */
        if( error < 0 && error != -3 && error != -4 ) return error;
      }
    }

    PRINT( "Used %d 2x2 pivots in factorization\n", ma27->info[13] );
    PRINT( "Factorization: %-d nonzeros in factors (density: %6.3f%%)\n",
           ma27->info[7],
           100.0*(double)ma27->info[7]/(double)((ma27->n)*(ma27->n)) );

    PRINT( " done\n" );

    return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "MA27_Solve"
  int MA27_Solve( Ma27_Data *ma27, double x[] ) {

    PRINT( " MA27 :: Solving..." );

    /* Unpack data structure and call MA27CD */
    MA27CD( &(ma27->n), ma27->factors, &(ma27->la), ma27->iw,
            &(ma27->liw), ma27->w, &(ma27->maxfrt), x,
            ma27->iw1, &(ma27->nsteps), ma27->icntl,
            ma27->info );

    PRINT( " done\n" );
    return ma27->info[0];
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "MA27_Refine"
  int MA27_Refine( Ma27_Data *ma27, double x[], double rhs[], double A[],
                   double tol, int maxitref ) {

    int    n = ma27->n, error, i, j, k, nitref;
    double b_norm, resid_norm;

    PRINT( " MA27 :: Performing iterative refinement...\n" );

    /* Compute initial residual */
    b_norm = cblas_dnrm_infty( n, rhs, 1 );        
    cblas_dcopy( n, rhs, 1, ma27->residual, 1 );
    for( k = 0; k < ma27->nz; k++ ) {
      i = ma27->irn[k];
      j = ma27->jcn[k];
      ma27->residual[i] -= A[k] * x[j];
      if( i != j ) ma27->residual[j] -= A[k] * x[i];
    }
    resid_norm = cblas_dnrm_infty( n, ma27->residual, 1 );

    PRINT( " Norm of residual: %8.2e", resid_norm );
    PRINT( " Tolerance: %8.2e\n", tol*(1+b_norm) );

    /* Perform iterative refinements, if required */
    /* Work with rtmp, and copy result into residual at the end */
    nitref = 0;
    while( nitref < maxitref && resid_norm > tol * (1+b_norm) ) {

      nitref++;

      /* Solve system again with residual as rhs */
      cblas_dcopy( n, ma27->residual, 1, rhs, 1 );
      MA27CD( &n, ma27->factors, &(ma27->la), ma27->iw,
              &(ma27->liw), ma27->w, &(ma27->maxfrt),
              rhs, ma27->iw1, &(ma27->nsteps),
              ma27->icntl, ma27->info );

      error = ma27->info[0];
      if( error ) return error;

      /* Update solution: x <- x + rhs */
      cblas_daxpy( n, 1.0, rhs, 1, x, 1 );
          
      /* Update residual: residual <- residual - A rhs */
      for( k = 0; k < ma27->nz; k++ ) {
        i = ma27->irn[k];
        j = ma27->jcn[k];
        ma27->residual[i] -= A[k] * rhs[j];
        if( i != j )
          ma27->residual[j] -= A[k] * rhs[i];
      }
      resid_norm = cblas_dnrm_infty( n, ma27->residual, 1 );
            
      PRINT("   Ref %-d: Norm of residual: %8.2e\n", nitref, resid_norm);
    }

    PRINT( " done\n" );
    return 0;
  }

  /* ================================================================= */

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "MA27_Finalize"
  void MA27_Finalize( Ma27_Data *ma27 ) {

    /* Free allocated memory */

    PRINT( " MA27 :: Deallocating data arrays..." );

    LBL_Free( ma27->irn );
    LBL_Free( ma27->jcn );
    LBL_Free( ma27->iw  );
    LBL_Free( ma27->ikeep );
    LBL_Free( ma27->iw1 );
    LBL_Free( ma27->factors );
    LBL_Free( ma27->residual );
    LBL_Free( ma27->w   );
    free( ma27 );
    PRINT( " done\n" );

    return;
  }

  /* ================================================================= */

#ifdef __cplusplus
}              /* Closing brace for  extern "C"  block */
#endif
