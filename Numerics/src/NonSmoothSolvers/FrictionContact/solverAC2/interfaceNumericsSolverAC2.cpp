#include "SiconosNumerics.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "solverAC2.h"

USING_PART_OF_NAMESPACE_EIGEN

void frictionContact3D_solverAC2(FrictionContactProblem* problem, double *reaction, double *velocity)
{
  // on reçoit le problem, vérifions d'abord qu'il est complet:
  if (!problem->isComplete)
    std::cout << "warning: the FrictionContactProblem has member isComplete equal to 0" << std::endl;
  // et de type SparseBlockStructuredMatrix:
  if (problem->M->storageType != 1)
    std::cout << "error: the FrictionContactProblem should have a matix of type 1 (SparseBlockStructuredMatrix)" << std::endl;
  // renommons, plus courtement
  NumericsMatrix * _MM = problem->M;
  SparseBlockStructuredMatrix * _M = _MM->matrix1;

  const int d = 3;
  const int n = problem->numberOfContacts;

  // =======================================
  // ========== construction de W ==========
  // =======================================
  Eigen::SparseMatrix<double, Eigen::RowMajor> W(_MM->size0, _MM->size1);
  // 1st pass: on compte le nb de non-zeros
  int number_of_nonzero = 0;
  for (int ibloc = 0; ibloc < _M->filled1 - 1; ibloc++) // ibloc : indice de bloc-ligne
  {
    int firstblock = _M->index1_data[ibloc];
    int lastblock  = _M->index1_data[ibloc + 1] - 1;
    // nrows: nb de lignes dans la bloc-ligne courante
    int nrows = (ibloc == 0) ? _M->blocksize0[0] : (_M->blocksize0[ibloc] - _M->blocksize0[ibloc - 1]);
    for (int k = firstblock; k <= lastblock; k++)   // k : indice de bloc
    {
      int jbloc = _M->index2_data[k];
      // ncols: nb de colonnes dans le bloc courant
      int ncols = (jbloc == 0) ? _M->blocksize1[0] : (_M->blocksize1[jbloc] - _M->blocksize1[jbloc - 1]);
      number_of_nonzero += nrows * ncols;
    }
  }
  //std::cout << "Number of nonzeros = " << number_of_nonzero << std::endl;

  // 2nd pass: fill the matrix
  W.startFill(number_of_nonzero); //nnz is optional
  for (int ibloc = 0; ibloc < _M->filled1 - 1; ibloc++) // ibloc : indice de bloc-ligne
  {
    int firstblock = _M->index1_data[ibloc];
    int lastblock  = _M->index1_data[ibloc + 1] - 1;
    // nrows: nb de lignes dans la bloc-ligne courante
    int nrows = (ibloc == 0) ? _M->blocksize0[0] : (_M->blocksize0[ibloc] - _M->blocksize0[ibloc - 1]);
    //std::cout << "la bloc-ligne courante " << ibloc << " possède " << nrows << " lignes" << std::endl;
    for (int i_in = 0; i_in < nrows; i_in++)   // i_in : indice de ligne _dans_le_bloc_courant
    {
      for (int k = firstblock; k <= lastblock; k++)   // k : indice de bloc
      {
        int jbloc = _M->index2_data[k];
        // ncols: nb de colonnes dans le bloc courant
        int ncols = (jbloc == 0) ? _M->blocksize1[0] : (_M->blocksize1[jbloc] - _M->blocksize1[jbloc - 1]);
        //std::cout << "la bloc-colonne courante " << jbloc << " possède " << ncols << " colonnes" << std::endl;
        for (int j_in = 0; j_in < ncols; j_in++)
        {
          int i = (ibloc == 0 ? 0 : _M->blocksize0[ibloc - 1]) + i_in;
          int j = (jbloc == 0 ? 0 : _M->blocksize1[jbloc - 1]) + j_in;
          double val = _M->block[k][i_in + j_in * nrows];
          //std::cout << "i    = " << i    << ", j    = " << j << std::endl;
          //std::cout << i << " " << j << " " << val << std::endl;
          W.fill(i, j) = val;
        }
      }
    }
  }
  W.endFill();

  Eigen::SparseMatrix<double, Eigen::ColMajor> Wcolmaj = Eigen::SparseMatrix<double, Eigen::ColMajor>(W); // duplicate, since we want col-major :(

  /*
  //DEBUG
  Eigen::SparseMatrix<double, Eigen::ColMajor> Wcolmaj(W.rows(), W.cols());
  Wcolmaj.startFill(W.nonZeros());
  for (int j=0; j<W.cols(); j++) {
    for (int i=0; i<W.rows(); i++) {
      double val = W.coeff(i,j);
      if (val != 0.0)
        Wcolmaj.fill(i,j) = val;
    }
  }
  Wcolmaj.endFill();
  */

  // =======================================
  // ========== construction de E ==========
  // =======================================
  double *E = new double[n * d * d];
  for (int i = 0; i < n; i++)
  {
    E[d * i  ] = 1;
    E[d * i  + n * d] = 0;
    E[d * i  + 2 * n * d] = 0;
    E[d * i + 1] = 0;
    E[d * i + 1 + n * d] = 1;
    E[d * i + 1 + 2 * n * d] = 0;
    E[d * i + 2] = 0;
    E[d * i + 2 + n * d] = 0;
    E[d * i + 2 + 2 * n * d] = 1;
  }

  // =========================================
  // ========== appel de solverAC2 ===========
  // =========================================
  solverAC2 mySolver(d, n, &Wcolmaj, problem->q, problem->mu, E, reaction);
  mySolver.solveDampedNewton();  // Newton step + Goldstein-Price line-search
  //mySolver.solveCylqp();

  Eigen::Map<VectorXd> u(velocity, n * d);
  Eigen::Map<VectorXd> r(reaction, n * d);
  Eigen::Map<VectorXd> q(problem->q, n * d);
  u = Wcolmaj * r + q;

  // clean up
  delete [] E;
}
