#include <stdio.h>

// ==== Dynamical System ====

// function to compute u
extern "C" void computeU(unsigned int* sizeOfU, unsigned int* sizeOfX, const double* time, double* xPtr, double* xDotPtr, double* UPtr)
{
  /* input parameter :
   *  sizeOfU (size of the vector u)
   *  sizeOfX (size of the vector x)
   *  time
   *  pointer to x
   *  pointer to xDot
   *  pointer to u (in-out parameter)
   */

  unsigned int size = *sizeOfU;
  for (unsigned int i = 0; i < size; i++)
  {
    if (i < *sizeOfX)
      UPtr[i] = *time * xPtr[i];
    else
      UPtr[i] = 0;
  }
}

// function to compute T
extern "C" void computeT(unsigned int* sizeOfU, unsigned int* sizeOfX, double* xPtr, double* TPtr)
{
  /* input parameter :
   *  sizeOfU (size of the vector u)
   *  sizeOfX (size of the vector X)
   *  pointer to x
   *  pointer to T (in-out parameter) - T is a (sizeOfX x sizeOfU) matrix
   */
  unsigned int sizeU = *sizeOfU;
  unsigned int sizeX = *sizeOfX;
  for (unsigned int j = 0; j < sizeU; j++)
  {
    for (unsigned int i = 0; i < sizeX; i++)
      TPtr[i + j * sizeX] =  xPtr[i];
  }
}

// vectorField
extern "C" void vectorField(int *sizeOfX, double *time, double *x, double *xdot)
{
  /* input parameter : sizeOfX (size of the vector X); time ; x (pointer to X vector);
   * output parameter : xdot (pointer to Xdot vector)
   */
  unsigned int size = *sizeOfX;
  for (unsigned int i = 0; i < size; i++)
    xdot[i] = *time * x[i];
}



extern "C" void computeJacobianX(double *time, double *x, double *jacob)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output parameter : jacob (pointer to JacobianX matrix)
   */

  printf("Call of the function 'computeJacobianX' of the basic plugin.\nYou have to implement this function.\n");
}


// ===== Lagrangian DS  =====

// Plugins for Fext, Fint, NNL (vectors), Mass, JacobianQNNL, JacobianVelocityNNL,
// JacobianQFint and JacobianVelocityFint (matrices)

extern "C" void computeFInt(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *fInt)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector); velocity (pointer to velocity vector);
   * output parameter : fInt (pointer to Fint *vector)
   */

  printf("Call of the function 'computeFInt' of the basic plugin.\nYou have to implement this function.\n");

}

extern "C" void computeFExt(unsigned int *sizeOfq, const double *time, double *fExt)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector);
   * output parameter : fExt (pointer to Fext vector)
   */

  printf("Call of the function 'computeFExt' of the basic plugin.\nYou have to implement this function.\n");

}

extern "C" void computeNNL(unsigned int *sizeOfq, double *q, double *velocity, double *Q)
{
  /* input parameter : sizeOfq (size of the vector q); q (pointer to q vector); velocity (pointer to velocity vector);
   * output parameter : Q (pointer to Q vector)
   */

  printf("Call of the function 'computeNNL' of the basic plugin.\nYou have to implement this function.\n");

}


extern "C" void computeMass(unsigned int *sizeOfq, const double *time, double *q, double *mass)
{
  /* input parameter : sizeOfq (size of the vector q); time ; q (pointer to q vector);
   * output parameter : mass (pointer to mass matrix)
   */
  printf("Call of the function 'computeMass' of the basic plugin.\nYou have to implement this function.\n");
}


extern "C" void computeJacobianQFInt(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *jacob)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector); velocity (pointer to velocity vector);
   * output parameter : jacob (pointer to JacobianCoordFint *matrix)
   */

  printf("Call of the function 'computeJacobianQFInt' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeJacobianVelocityFInt(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *jacob)
{
  printf("Call of the function 'computeJacobianVelocityFInt' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeJacobianQNNL(unsigned int *sizeOfq, double *q, double *velocity, double *jacob)
{
  printf("Call of the function 'computeJacobianQNNL' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void computeJacobianVelocityNNL(unsigned int *sizeOfq, double *q, double *velocity, double *jacob)
{
  printf("Call of the function 'computeJacobianVelocityNNL' of the basic plugin.\nYou have to implement this function.\n");
}

// plugin for relations
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

// ===== Linear DS  ====

// Plugins for A, B (matrices), u and f. See LinearDS.h

extern "C" void computeA(unsigned int *sizeOfA, double* APtr, const double *time)
{
  /* input parameter : time, sizeOfA (size of the matrix A)
   * output parameter : APtr (pointer to SiconosMatrix)
   */
  printf("Call of the function 'computeA' of the basic plugin.\nYou have to implement this function.\n");

}

extern "C" void computeB(unsigned int *sizeOfB, double* b, const double *time)
{
  /* input parameter : time, sizeOfB (size of the vector b);
   * output parameter : b (pointer to b vector)
   */
  printf("Call of the function 'computeB' of the basic plugin.\nYou have to implement this function.\n");

}
