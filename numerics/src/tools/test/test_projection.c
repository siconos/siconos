#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <string.h>
#include <projectionOnRollingCone.h>

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif


static double orthogonality(double * reaction, double * reaction_ori)
{
  double ortho =0.0;
  for (int i =0; i<5; i++)
  {
    ortho += (reaction_ori[i] - reaction[i])*(reaction[i]);
  }
  return ortho;
}

static double compute_diff(double * reaction, double * reaction_again)
{
  double diff =0.0;
  for (int i =0; i<5; i++)
  {
    diff += (reaction_again[i] - reaction[i])*(reaction_again[i] - reaction[i]);
  }
  return sqrt(diff);
}


static int test_projection(double * reaction, double mu, double mur)
{
  printf("################################ \n");
  printf("start test \n");
  printf("################################ \n");
  unsigned int status;
  unsigned int status_again;
  double reaction_ori[5];
  double reaction_again[5];
  double ortho=0.0, diff =0.0;
  DEBUG_EXPR(NV_display(reaction,5););

  for (int i =0; i<5; i++)  reaction_ori[i] = reaction[i];
  status = projectionOnRollingCone(reaction,mu,mur);
  DEBUG_EXPR(NV_display(reaction,5););
  display_status_rolling_cone(status);
  ortho = orthogonality(reaction,reaction_ori);
  printf("ortho = %f\n", ortho);
  assert(ortho < 1e-14);

  for (int i =0; i<5; i++)  reaction_again[i] = reaction[i];
  status_again = projectionOnRollingCone(reaction_again,mu,mur);
  DEBUG_EXPR(NV_display(reaction_again,5););
  display_status_rolling_cone(status_again);
  diff = compute_diff(reaction,reaction_again);
  printf("diff = %f\n", diff);
  assert(diff < 1e-14);

  return status;

}
int main(void)
{

  int info=0;
  unsigned int status;
  double reaction[5];
  double mu=1.0, mur=1.0;


  reaction[0] = 1.0;
  reaction[1] = 0.0;
  reaction[2] = 0.0;
  reaction[3] = 0.0;
  reaction[4] = 0.0;

  status = test_projection(reaction, mu, mur);
  if (status != PROJRCONE_INSIDE)
  {
    info+=1;
  }
  reaction[0] = -1.0;
  reaction[1] = 0.5;
  reaction[2] = 0.0;
  reaction[3] = 0.5;
  reaction[4] = 0.0;
  status = test_projection(reaction, mu, mur);
  if (status != PROJRCONE_DUAL)
  {
    info+=1;
  }

  reaction[0] = 1.0;
  reaction[1] = 2.0;
  reaction[2] = 0.0;
  reaction[3] = 0.0;
  reaction[4] = 0.0;
  status = test_projection(reaction, mu, mur);
  if (status != PROJRCONE_BOUNDARY_FRICTION)
  {
    info+=1;
  }

  reaction[0] = 1.0;
  reaction[1] = 0.0;
  reaction[2] = 0.0;
  reaction[3] = 2.0;
  reaction[4] = 0.0;
  status = test_projection(reaction, mu, mur);
  if (status != PROJRCONE_BOUNDARY_ROLLING)
  {
    info+=1;
  }

  reaction[0] = 1.0;
  reaction[1] = 2.0;
  reaction[2] = 0.0;
  reaction[3] = 2.0;
  reaction[4] = 0.0;
  status = test_projection(reaction, mu, mur);
  if (status != PROJRCONE_BOUNDARY_FRICTION_ROLLING)
  {
    info+=1;
  }

  reaction[0] = 1.0;
  reaction[1] = 2.0;
  reaction[2] = 1.0;
  reaction[3] = 2.0;
  reaction[4] = 1.0;
  status = test_projection(reaction, mu, mur);
  if (status != PROJRCONE_BOUNDARY_FRICTION_ROLLING)
  {
    info+=1;
  }

  return info;

}
