/*!\file pfc_2D_projc.c
 *
 * \fn  pfc_2D_projc( int n , double mu , double *z , double *p , int *status )
 *
 * pfc_2D_projc is a specific projection operator related to CPG (conjugated projected gradient) algorithm
 *              for primal contact problem with friction.\n
 *
 * Ref: Renouf, M. and Alart, P. "" Comp. Method Appl. Mech. Engrg. (2004).
 *
 * \param n       Unchanged parameter which represents the half dimension of the system.
 * \param mu      Unchanged parameter which represents the friction coefficient
 * \param z       Modified parameter which retruns the corrected iterate.
 * \param p       Unchanged parameter which contains the components of the descent direction.
 * \param status  Unchanged parameter which contains the vector status
 *
 * \author Mathieu Renouf.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void pfc_2D_projc(int nc , double mu , double *z , double *p , int *status)
{

  int i;

  for (i = 0 ; i < nc ; ++i)
  {

    /* No contact case */

    if (z[2 * i] < 0.0)
    {

      z[2 * i  ]  = 0.0;
      z[2 * i + 1]  = 0.0;
      status[i] = 0;
    }
    else
    {
      /* contact case */
      if (p[2 * i + 1] == 0.0)
      {
        /* sliding contact */
        if (z[2 * i + 1] > 0.0)
        {
          z[2 * i + 1] = mu * z[2 * i];
          status[i] = 2;
        }
        else
        {
          z[2 * i + 1] = -mu * z[2 * i];
          status[i] = 3;
        }
      }
      else
      {
        /* slide forward */
        if (z[2 * i + 1] < -mu * z[2 * i])
        {
          z[2 * i + 1]  = -mu * z[2 * i];
          status[i] = 3;
        }
        /* slide backward */
        else if (z[2 * i + 1] > mu * z[2 * i])
        {
          z[2 * i + 1] = mu * z[2 * i];
          status[i] = 2;
        }
        /* sticking contact */
        else
        {
          status[i] = 1;
        }
      }
    }
  }
}

