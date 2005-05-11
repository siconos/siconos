#include <stdio.h>


extern "C" void vectorField(int *sizeOfX, double *time, double *x, double *xdot)
{
  /* input parameter : sizeOfX (size of the vector X); time ; x (pointer to X vector);
   * output parameter : xdot (pointer to Xdot vector)
   */
  printf("Call of the function 'vectorField' of the basic plugin.\nYou have to implement this function.\n");

  // address of one element (double) of X
  //double **xElem;
  //xElem = &x;

  // address of one element (double) of Xdot
  //double **xdotElem;
  //xdotElem = &xdot;

  // to access to the next element, you must do :
  // xElem++;
  // xdotElem++;

}

extern "C" void computeMass(int *sizeOfq, double *time, double *q, double *mass)
{
  /* input parameter : sizeOfq (size of the vector q); time ; q (pointer to q vector);
   * output parameter : mass (pointer to mass matrix)
   */
  printf("Call of the function 'computeMass' of the basic plugin.\nYou have to implement this function.\n");

  // address of one element (double) of q
  //double **qElem;
  //qElem = &q;

  // address of one element (double) of mass
  //double **massElem;
  //massElem = &mass;

  // to access to the next element, you must do :
  // qElem++;
  // massElem++;

}

extern "C" void computeQNLInertia(int *sizeOfq, double *q, double *velocity, double *Q)
{
  /* input parameter : sizeOfq (size of the vector q); q (pointer to q vector); velocity (pointer to velocity vector);
   * output parameter : Q (pointer to Q vector)
   */

  printf("Call of the function 'computeQNLInertia' of the basic plugin.\nYou have to implement this function.\n");

}

extern "C" void computeFInt(int *sizeOfq, double *time, double *q, double *velocity, double *fInt)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector); velocity (pointer to velocity vector);
   * output parameter : fInt (pointer to Fint *vector)
   */

  printf("Call of the function 'computeFInt' of the basic plugin.\nYou have to implement this function.\n");

}

extern "C" void computeFExt(int *sizeOfq, double *time, double *fExt)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector);
   * output parameter : fExt (pointer to Fext vector)
   */

  printf("Call of the function 'computeFExt' of the basic plugin.\nYou have to implement this function.\n");

}

extern "C" void computeJacobianX(double *time, double *x, double *jacob)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output parameter : jacob (pointer to JacobianX matrix)
   */

  printf("Call of the function 'computeJacobianX' of the basic plugin.\nYou have to implement this function.\n");
}


extern "C" void computeJacobianQFInt(int *sizeOfq, double *time, double *q, double *velocity, double *jacob)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector); velocity (pointer to velocity vector);
   * output parameter : jacob (pointer to JacobianCoordFint *matrix)
   */

  printf("Call of the function 'computeJacobianQFInt' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeJacobianVelocityFInt(int *sizeOfq, double *time, double *q, double *velocity, double *jacob)
{
  printf("Call of the function 'computeJacobianVelocityFInt' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeJacobianQQNLInertia(int *sizeOfq, double *q, double *velocity, double *jacob)
{
  printf("Call of the function 'computeJacobianQQNLInertia' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeJacobianVelocityQNLInertia(int *sizeOfq, double *q, double *velocity, double *jacob)
{
  printf("Call of the function 'computeJacobianVelocityQNLInertia' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeOutput(int *sizeOfX, double* x, double *time, double* lambda, double* y)
{
  /* input parameter : sizeOfX (size of the vector X); x (pointer to x vector); time; lambda (pointer to lambda vector)
   * output parameter : y (pointer to vector y )
   */
  printf("Call of the function 'computeOutput' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeInput(int *sizeOfX, double* x, double *time, double* lambda, double* r)
{
  /* input parameter : sizeOfX (size of the vector X); x (pointer to x vector); time; lambda (pointer to lambda vector)
   * output parameter : r (pointer to vector r )
   */
  printf("Call of the function 'computeInput' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeU(int *sizeOfU, double* u, double *time/* time */)
{
  /* input parameter : time
   * output patarmeter : sizeOfU (size of the vector u); u (pointer to u vector)
   */
  printf("Call of the function 'computeU' of the basic plugin.\nYou have to implement this function.\n");

}

extern "C" void computeF(int *sizeOfF, double* f, double *time/* time */)
{
  /* input parameter : time
   * output patarmeter : sizeOfF (size of the vector f); f (pointer to f vector)
   */
  printf("Call of the function 'computeF' of the basic plugin.\nYou have to implement this function.\n");

}
