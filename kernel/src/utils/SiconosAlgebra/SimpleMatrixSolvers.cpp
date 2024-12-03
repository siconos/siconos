/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithm>

#include "BlockMatrix.hpp"
#include "NumericsToolsNamespace.h"
#include "SiconosConfig.h"
#include "SiconosLapack.h"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "cholesky.hpp"

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

#ifdef DEBUG_MESSAGES
#include <cs.h>
#endif

// namespace lapack = boost::numeric::bindings::lapack;

void siconos::algebra::solveInPlace(SiconosMatrix &A, SiconosVector &B) {
  assert(A.cols() == B.rows());
  if (A.lu_siconos == nullptr) {
    A.lu_siconos.reset(new Eigen::FullPivLU<SiconosMatrix>(A));
    // A.isFactorized = true;
    B = A.lu_siconos->solve(B);  // TODO : avoid temp copy
    // std::cout << "solve results " << std::endl << tmp << std::endl;
  } else {
    B = A.lu_siconos->solve(B);
  }
}

void siconos::algebra::solveInPlace(SiconosMatrix &A, SiconosMatrix &B) {
  assert(A.rows() == B.rows());
  // std::cout << "A " << std::endl << A << std::endl;

  if (A.lu_siconos == nullptr) {
    A.lu_siconos.reset(new Eigen::FullPivLU<SiconosMatrix>(A));
    // A.isFactorized = true;
    B = A.lu_siconos->solve(B);  // TODO : avoid temp copy
    // std::cout << "solve results " << std::endl << tmp << std::endl;
  } else {
    SiconosMatrix tmp = A.lu().solve(B);
    B = A.lu_siconos->solve(B);
  }
  // std::cout << "A  apres" << std::endl << A << std::endl;
}

void siconos::algebra::solveInPlace(SiconosMatrix &A, MapType &B) {
  assert(A.rows() == B.rows());
  // std::cout << "A " << std::endl << A << std::endl;

  if (A.lu_siconos == nullptr) {
    A.lu_siconos.reset(new Eigen::FullPivLU<SiconosMatrix>(A));
    // A.isFactorized = true;
    B = A.lu_siconos->solve(B);  // TODO : avoid temp copy
    // std::cout << "solve results " << std::endl << tmp << std::endl;
  } else {
    SiconosMatrix tmp = A.lu().solve(B);
    B = A.lu_siconos->solve(B);
  }
  // std::cout << "A  apres" << std::endl << A << std::endl;
}

void siconos::algebra::solveByLeastSquares(SiconosMatrix &A, SiconosMatrix &B) {
  int info = 0;
  // #ifdef USE_OPTIMAL_WORKSPACE
  //   info += lapack::gels(*mat.Dense, *(B.dense()), lapack::optimal_workspace());
  // #endif
  // #ifdef USE_MINIMAL_WORKSPACE
  //   info += lapack::gels(*mat.Dense, *(B.dense()), lapack::minimal_workspace());
  // #endif
  std::cout << B << std::endl;
  int M_ = A.rows();
  int N_ = B.cols();
  int NRHS_ = B.cols();
  int LDA_ = std::max(1, M_);
  int LDB_ = std::max(LDA_, N_);

  DGELS('T', M_, N_, NRHS_, A.data(), LDA_, B.data(), LDB_, &info);
  if (info != 0) THROW_EXCEPTION("SiconosMatrix::SolveByLeastSquares failed.");
  std::cout << B << std::endl << "************************" << std::endl;
}

void siconos::algebra::solveByLeastSquares(SiconosMatrix &A, SiconosVector &B) {
  // std::cout << B << std::endl;
  SiconosMatrix tmpB(B.size(), 1);
  tmpB.col(0) = B;  // Conversion of vector to matrix. Temporary solution.
  std::cout << "tpmB" << std::endl << tmpB << std::endl;
  int info = 0;

  // #ifdef USE_OPTIMAL_WORKSPACE
  //   info += lapack::gels(*mat.Dense, tmpB, lapack::optimal_workspace());
  // #endif
  // #ifdef USE_MINIMAL_WORKSPACE
  //   info += lapack::gels(*mat.Dense, tmpB, lapack::minimal_workspace());
  // #endif

  int M_ = A.rows();
  int N_ = B.cols();
  int NRHS_ = B.cols();
  int LDA_ = std::max(1, M_);
  int LDB_ = std::max(LDA_, N_);

  DGELS('T', M_, N_, NRHS_, A.data(), LDA_, tmpB.data(), LDB_, &info);
  if (info != 0) {
    std::cout << "info = " << info << std::endl;
    THROW_EXCEPTION("SiconosMatrix::SolveByLeastSquares failed.");
  } else {
    B = tmpB.col(0);
    std::cout << "B out :" << std::endl << B << "************************" << std::endl;
  }
}
