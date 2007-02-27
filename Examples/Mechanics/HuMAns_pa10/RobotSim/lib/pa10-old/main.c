//#include "robotmodel.h"
/* gcc main.c inertie.c -o mm -lm */
#define N_DOF 3
main()
{
  double q[N_DOF] = {0, 0, 0};
  double M[N_DOF * N_DOF] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  //modele_inertie(q,M);

  printf("-> %lf %lf %lf \n", M[0], M[1], M[2]);
  printf("-> %lf %lf %lf \n", M[3], M[4], M[5]);
  printf("-> %lf %lf %lf \n\n", M[6], M[7], M[8]);

  printf("-> %lf %lf %lf \n", M[0, 0], M[0, 1], M[0, 2]);
  printf("-> %lf %lf %lf \n", M[1, 0], M[1, 1], M[1, 2]);
  printf("-> %lf %lf %lf \n\n", M[2, 0], M[2, 1], M[2, 2]);
  // (3*nligne)*i +j

}
