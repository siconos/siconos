#include "pinv.h"
#include "cond.h"
#include <stdlib.h>
#include <stdio.h>

#include <string.h>
#include <projectionOnRollingCone.h>

static void display(unsigned int status)
{
  printf("status = %i\n", status);
  if (status == PROJRCONE_INSIDE)
  {
    printf("PROJRCONE_INSIDE reaction was inside the cone\n");
  }
  else if (status == PROJRCONE_DUAL)
  {
    printf("PROJRCONE_DUAL reaction was inside the dual cone\n");
  }
  else if (status == PROJRCONE_BOUNDARY_FRICTION)
  {
    printf("PROJRCONE_BOUNDARY_FRICTION reaction is projected on the friction cone\n");
  }
  else if (status == PROJRCONE_BOUNDARY_ROLLING)
  {
    printf("PROJRCONE_BOUNDARY_ROLLING reaction is projected on the rolling cone\n");
  }
  else if (status == PROJRCONE_BOUNDARY_FRICTION_ROLLING)
  {
    printf("PROJRCONE_BOUNDARY_FRICTION_ROLLING reaction is projected on the both cones\n");
  }
  

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
  status = projectionOnRollingCone(reaction,mu,mur);
  display(status);
  if (status != PROJRCONE_INSIDE)
  {
    info+=1;
  }

  reaction[0] = 1.0;
  reaction[1] = 2.0;
  reaction[2] = 0.0;
  reaction[3] = 0.0;
  reaction[4] = 0.0;
  status = projectionOnRollingCone(reaction,mu,mur);
  display(status);
  
  reaction[0] = 1.0;
  reaction[1] = 0.0;
  reaction[2] = 0.0;
  reaction[3] = 2.0;
  reaction[4] = 0.0;
  status = projectionOnRollingCone(reaction,mu,mur);
  display(status);
  
  reaction[0] = 1.0;
  reaction[1] = 2.0;
  reaction[2] = 0.0;
  reaction[3] = 2.0;
  reaction[4] = 0.0;
  status = projectionOnRollingCone(reaction,mu,mur);
  display(status);
  
  reaction[0] = -1.0;
  reaction[1] = 0.5;
  reaction[2] = 0.0;
  reaction[3] = 0.5;
  reaction[4] = 0.0;
  status = projectionOnRollingCone(reaction,mu,mur);
  display(status); 
  return info;

}
