#include "NaturalMapGenerated.h"
#include "assert.h"                       // for assert
#include "fc3d_NaturalMapABGenerated.h"   // for fc3d_NaturalMapABGenerated
#include "fc3d_NaturalMapFABGenerated.h"  // for fc3d_NaturalMapFABGenerated
#include "fc3d_NaturalMapFGenerated.h"    // for fc3d_NaturalMapFGenerated
#include "op3x3.h"                        // for cpy3x3, cpy3, SET3

void fc3d_NaturalMapFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B)
{
  double result[21];

  assert(reaction);
  assert(velocity);
  assert(rho);

  SET3(reaction);
  SET3(velocity);
  SET3(rho);


  if(f && A && B)
  {

    fc3d_NaturalMapFABGenerated(
      *reaction0, *reaction1, *reaction2,
      *velocity0, *velocity1, *velocity2,
      mu,
      *rho0, *rho1, *rho2,
      result);
    cpy3(result, f);
    cpy3x3(result + 3, A);
    cpy3x3(result + 12, B);
  }

  else
  {
    if(f)
    {
      fc3d_NaturalMapFGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3(result, f);
    }

    if(A && B)
    {
      fc3d_NaturalMapABGenerated(
        *reaction0, *reaction1, *reaction2,
        *velocity0, *velocity1, *velocity2,
        mu,
        *rho0, *rho1, *rho2,
        result);
      cpy3x3(result, A);
      cpy3x3(result + 9, B);
    }
  }
}
