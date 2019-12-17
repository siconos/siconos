/*
 * ex4xc_cb.c
 *
 * C implementation of the external equations for example 4
 * (DEA, parameter estimation) in the file ex4x.gms.
 * Here, we use the callback capability for the messages.
 * The special aspect of this implementation is that we tell the
 * solvers that some of the derivatives are constant and we declare
 * the derivate w.r.t cv to be constant = +1.
 *
 */

#if 0

#include <math.h>
#include <assert.h>

#include "geheader.h"

#include "SiconosBlas.h"

#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Gams.h"


#define BOGUS_EXTEQ 2

static int nchars;
static int mmode;
/* use the do-while loop so the macro acts like one statement */
#define MSGCB(mode,buf) do { mmode = mode; nchars=strlen(buf); msgcb(&mmode,&nchars,buf,nchars);} while (0)


/* Goal of this plugin: we want to compute v from the problem data and r */

/* Some general comments:
 * - We need some data that are not going to be provided by the solvers: the
 *   matrix M, H and the vector f -> we have to resort using a static pointer
 *   to the problem data
 *
 * - gams will call this plugin for each component of v.
 */


/* TODO use thread local storage here ?  */
static GlobalFrictionContactProblem* GFCP;

void set_pointer_data(GlobalFrictionContactProblem* gfcp)
{
  GFCP = gfcp;
}

static inline double solve_lu(NumericsMatrix* M, unsigned i, rhs)
{
  nb_matrix = div(i, 3);
}

GE_API int GE_CALLCONV
gefunc(int *icntr, double *x, double *func, double *d, msgcb_t msgcb)
{

  assert(GFCP && "No GlobalFrictionContactProblem has been set!");

#if defined(SOME_DEBUG_STUFF)
  char msgBuf[256];
#endif
  double t,  dtdh;
  double t1, dt1dh;
  double t2, dt2dh;
  double h, cv, dfdh, dfdcv, f;

  if(icntr[I_Mode] == DOINIT)
  {
    /*
     * Initialization Mode:
     * Write a "finger print" to the status file so errors in the DLL
     * can be detected more easily. This should be done before anything
     * can go wrong. Also write a line to the log just to show it.
     */
    MSGCB(LOGFILE | STAFILE,"");
    MSGCB(LOGFILE | STAFILE,"--- GEFUNC in ex4xc_cb.c is being initialized.");

    /*  Test the equation count and return 2 if bogus */
    if(icntr[I_Neq] != 1)
    {
      MSGCB(LOGFILE | STAFILE,
            "--- Model has the wrong number of external equations.");
      return BOGUS_EXTEQ;
    }
    if(2 != icntr[I_Nz])
    {
      MSGCB(LOGFILE | STAFILE,
            "--- The external equation should be fully dense.");
      return BOGUS_EXTEQ;
    }
    if(2 != icntr[I_Nvar])
    {
      MSGCB(LOGFILE | STAFILE,
            "--- The external equation should have 2 variables.");
      return BOGUS_EXTEQ;
    }
    /* Define number of constant derivatives    */
    icntr[I_ConstDeriv] = 1;

    /* Form   */
    NM_setup(GFCP->H);
    return 0;
  } /* initialization mode */
  else if(icntr[I_Mode] == DOCONSTDERIV)
  {
    assert(0 && "Why are we here? icntr[I_Mode] == DOCONSTDERIV ");
  }
  else if(icntr[I_Mode] == DOTERM)
  {
    /* Termination mode: free allocated storage */

    return 0;
  } /* termination mode */
  else if(icntr[I_Mode] == DOEVAL)
  {
    /*
     * Function index: there is only one equation here,
     * but we check the equation number just to show the principle.
     */
    if(icntr[I_Eqno] != 1)
    {
      MSGCB(STAFILE | LOGFILE," ** Eqno has unexpected value.");
      return BOGUS_EXTEQ;
    }

    NumericsMatrix Mlu = ((GFC3D_Gams*) GFCP->env)->Mlu;
    double* rhs = ((GFC3D_Gams*) GFCP->env)->rhs;
    /* set rhs = Hr + f */
    cblas_dcopy(n, func, 1, rhs, 1);



    /* solve Mv = Hr + f = rhs */

    /* get our values from the array passed in, just to be neat */
    h = x[0];
    cv = x[1];

#if defined(SOME_DEBUG_STUFF)
    sprintf(msgBuf, "              dh = %g, dcv = %g, f = %f",
            dfdh, dfdcv, f);
    MSGCB(STAFILE | LOGFILE, msgBuf);
#endif

    if(icntr[I_Dofunc])
    {
      *func = f;
    }

    if(icntr[I_Dodrv])
    {
      assert(0 && "Computing derivative is not implemented yet ...!");
    }
    return 0;
  } /* Function and Derivative Evaluation Mode */
  else
  {
    MSGCB(STAFILE | LOGFILE, " ** Mode not defined.");
    return 2;
  }
} /* gefunc */

#endif
