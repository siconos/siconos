//#include "robotmodel.h"
/* gcc test.c NLEffects.c -o tt -lm */
#define N_DOF 3

main()
{
  double q[N_DOF] = {0, 0};
  double qdot[N_DOF] = {0, 0};
  double NL[N_DOF] = {0, 0};

  NLEffects(NL, q, qdot);

  printf("-> %lf %lf \n", NL[0], NL[1]);

}
