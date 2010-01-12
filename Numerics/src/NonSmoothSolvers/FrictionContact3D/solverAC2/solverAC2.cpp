/* solverAC2
A class to solve the linear incremental friction problem.
F. Cadoux Apr. 20th, 2009 (INRIA).
*/

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/Sparse>

#include "solverAC2.h"

#define BAVARD

solverAC2::solverAC2(
  unsigned int d_in,
  unsigned int NObj,
  const unsigned int * ndof,
  const double *const * MassMat,
  const double * f_in,
  unsigned int n_in,
  const double * mu_in,
  const double * E_in,
  const double * w_in,
  const int *ObjA,
  const int *ObjB,
  const double *const HA[],
  const double *const HB[],
  double r0[],
  double v0[]) : _LUFailed(false), printFormat(6, Eigen::AlignCols, ", ", "\n", "[", "]"), d(d_in), m(_m), n(n_in), primalAvailable(true), v(v0), r(r0)
{

  // compute the total number of degrees of freedom and the first index of each object
  _m = ndof[0];
  unsigned int nnzM_sparse = ndof[0] * ndof[0];
  std::vector<unsigned int> startObj(NObj);
  startObj[0] = 0;
  for (unsigned int i = 1; i < NObj; i++)
  {
    _m += ndof[i];
    nnzM_sparse += ndof[i] * ndof[i];
    startObj[i] = startObj[i - 1] + ndof[i - 1];
  }

  // alloc for sparse mats
  M_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(m, m);
  //H_sparse = new Eigen::SparseMatrix<double,Eigen::ColMajor>(n*d, m);
  W_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);
  dfAC_du_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);
  dfAC_dr_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);
  dG_dr_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);

  // assemble M_sparse
  M_sparse->startFill(nnzM_sparse);
  for (unsigned int k = 0; k < NObj; k++)
    for (int j = 0; j < ndof[k]; j++)
      for (int i = 0; i < ndof[k]; i++)
        M_sparse->fill(startObj[k] + i, startObj[k] + j) = MassMat[k][i + ndof[k] * j];
  M_sparse->endFill();
  // check
  // printSparse(*M_sparse);

  // map of f
  f = Eigen::Map<VectorXd>(f_in, m);

  // assemble H
  H = MatrixXd::Zero(n * d, m);
  for (unsigned int i = 0; i < n; i++)
  {
    int jA = ObjA[i];
    int jB = ObjB[i];
    H.block(i * d, startObj[jA], d, ndof[jA]) += Eigen::Map<MatrixXd>(HA[i], d, ndof[jA]);
    if (jB != -1) // TODO: throw exception if jB makes no sense (<0 or >NObj)
      H.block(i * d, startObj[jB], d, ndof[jB]) -= Eigen::Map<MatrixXd>(HB[i], d, ndof[jB]);
  }

  /*
    // assemble H_sparse (to be used only if LU decomp accepts a sparse rhs)
    H_sparse->startFill();
    for (int k = 0; k<NObj; k++) { // k: index of object
      for(int j=0; j<ndof[k]; j++) // startObj[k]+j iterates over columns
        for (int ic=0; ic<n; ic++) // ic: index of contact
          for (int i=0; i<d; i++) { // d*ic+i iterates over rows
            if (ObjA[ic]==k && ObjB[ic]!=k) { // copy HA[i]
              H_sparse->fill(ic*d+i,startObj[k]+j) = HA[ic][i+d*j];
            }
            if (ObjA[ic]!=k && ObjB[ic]==k) { // copy -HB[i]
              H_sparse->fill(ic*d+i,startObj[k]+j) = -HB[ic][i+d*j];
            }
            if (ObjA[ic]==k && ObjB[ic]==k) { // copy -HB[i]
              H_sparse->fill(ic*d+i,startObj[k]+j) = HA[ic][i+d*j]-HB[ic][i+d*j];
            }
          }
    }
    H_sparse->endFill();
    printSparse(*H_sparse);
  */

  // map E, w, mu
  E = Eigen::Map<MatrixXd>(E_in, n * d, d);
  w = Eigen::Map<VectorXd>(w_in, n * d);
  mu = Eigen::Map<VectorXd>(mu_in, n);

  // W sparse
  // M should be factored by blocks and it should be a Cholesky, not LU, factorization
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatXd;
  luOfM_sparse = new Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::UmfPack>;
  luOfM_sparse->compute(*M_sparse);

  Htrans = MatrixXd(H.transpose());
  //Htrans_sparse = SparseMatXd(H_sparse->transpose());
  Eigen::MatrixXd temp2 = Eigen::MatrixXd::Zero(m, n * d);
  Eigen::MatrixXd W2;

  if (luOfM_sparse->succeeded())
  {
    if (luOfM_sparse->solve(Htrans, &temp2))
    {
      W2 = H * temp2;
    }
    else
    {
      std::cout << "solve failed" << std::endl;
    }
  }
  else
  {
    std::cout << "sparse LU of M failed" << std::endl;
  }
  q.resize(n * d);
  Eigen::VectorXd temp5 = Eigen::VectorXd::Zero(m);
  luOfM_sparse->solve(f, &temp5); // temp5 = M\f
  VectorXd temp3 = H * temp5;
  q = w - temp3;

  //compute W_sparse and build Id at the same time
  W_sparse->startFill();
  Id = Eigen::SparseMatrix<double>(n * d, n * d);
  Id.startFill();
  for (int j = 0; j < n * d; j++)
  {
    for (int i = 0; i < n * d; i++)
    {
      if (W2(i, j) > 1e-6 || W2(i, j) < -1e-6)
        W_sparse->fill(i, j) = W2(i, j);
    }
    Id.fill(j, j) = 1;
  }
  W_sparse->endFill();
  Id.endFill();

  // init rhoN, rhoT, G, dr
  rhoN = std::vector<double>(n, 1.0);
  rhoT = std::vector<double>(n, 1.0);
  G.resize(n * d);
  u.resize(n * d);
  dr.resize(n * d);
  rhs_for_lusolve.resize(m);

  // init for fAC and its derivatives
  fAC.resize(n * d);
  dfAC_du_packed.resize(n * d, d);
  dfAC_dr_packed.resize(n * d, d);

  // default param value
  setNIterMax(25);
  setNIterMax_LineSearch(20);
  setTol(1e-5);
}

solverAC2::solverAC2(unsigned int d_in, unsigned int n_in, Eigen::SparseMatrix<double> * W_sparse_in,
                     const double * q_in, const double * mu_in, const double * E_in, double r0[]) : d(d_in), n(n_in), m(_m), primalAvailable(false)
{


  primalAvailable = false;
  r = r0;
  v = NULL;
  W_sparse = W_sparse_in;

  q.resize(n * d);
  q = Eigen::VectorXd::Map(q_in, n * d);


  // alloc for sparse mats
  dfAC_du_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);
  dfAC_dr_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);
  dG_dr_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);

  // map E, mu
  E = Eigen::Map<MatrixXd>(E_in, n * d, d);
  mu = Eigen::Map<VectorXd>(mu_in, n);

  //compute W_sparse and build Id at the same time
  Id = Eigen::SparseMatrix<double>(n * d, n * d);
  Id.startFill(n * d);
  for (int j = 0; j < n * d; j++)
    Id.fill(j, j) = 1;
  Id.endFill();

  // init rhoN, rhoT, G, dr
  rhoN = std::vector<double>(n, 1.0);
  rhoT = std::vector<double>(n, 1.0);
  G.resize(n * d);
  u.resize(n * d);
  dr.resize(n * d);

  // init for fAC and its derivatives
  fAC.resize(n * d);
  dfAC_du_packed.resize(n * d, d);
  dfAC_dr_packed.resize(n * d, d);

  // default param value
  setNIterMax(200);
  setNIterMax_LineSearch(20);
  setTol(1e-5);

}

// constructor which does not fill the data (readFromFile is expected to be used afterwards)
solverAC2::solverAC2(unsigned int d_in, unsigned int n_in, unsigned int m_in, const char* filename, double r0[], double v0[])
  : d(d_in), n(n_in), m(_m), primalAvailable(true), v(v0), r(r0)
{

  _m = m_in;
  // alloc for sparse mats
  M_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(m, m);
  //H_sparse = new Eigen::SparseMatrix<double,Eigen::ColMajor>(n*d, m);
  W_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);
  dfAC_du_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);
  dfAC_dr_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);
  dG_dr_sparse = new Eigen::SparseMatrix<double, Eigen::ColMajor>(n * d, n * d);

  // init rhoN, rhoT, G, dr
  rhoN = std::vector<double>(n, 1.0);
  rhoT = std::vector<double>(n, 1.0);
  G.resize(n * d);
  u.resize(n * d);
  dr.resize(n * d);
  rhs_for_lusolve.resize(m);

  // init for fAC and its derivatives
  fAC.resize(n * d);
  dfAC_du_packed.resize(n * d, d);
  dfAC_dr_packed.resize(n * d, d);

  // set correct sizes
  f.resize(m);
  H.resize(n * d, m);
  w.resize(n * d);
  E.resize(n * d, d);
  mu.resize(n);

  // default param value
  setNIterMax(50);
  setNIterMax_LineSearch(20);
  setTol(1e-5);

  // read data
  readFromFile(filename);

  // W sparse
  // M should be factored by blocks and it should be a Cholesky, not LU, factorization
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatXd;
  luOfM_sparse = new Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::UmfPack>;
  luOfM_sparse->compute(*M_sparse);

  Htrans = MatrixXd(H.transpose());
  //Htrans_sparse = SparseMatXd(H_sparse->transpose());
  Eigen::MatrixXd temp2 = Eigen::MatrixXd::Zero(m, n * d);
  Eigen::MatrixXd W2;

  if (luOfM_sparse->succeeded())
  {
    if (luOfM_sparse->solve(Htrans, &temp2))
    {
      W2 = H * temp2;
    }
    else
    {
      std::cout << "solve failed" << std::endl;
    }
  }
  else
  {
    std::cout << "sparse LU of M failed" << std::endl;
  }
  q.resize(n * d);
  Eigen::VectorXd temp5 = Eigen::VectorXd::Zero(m);
  luOfM_sparse->solve(f, &temp5); // temp5 = M\f
  VectorXd temp3 = H * temp5;
  q = w - temp3;

  //compute W_sparse and build Id at the same time
  W_sparse->startFill();
  Id = Eigen::SparseMatrix<double>(n * d, n * d);
  Id.startFill();
  for (int j = 0; j < n * d; j++)
  {
    for (int i = 0; i < n * d; i++)
    {
      if (W2(i, j) > 1e-6 || W2(i, j) < -1e-6)
        W_sparse->fill(i, j) = W2(i, j);
    }
    Id.fill(j, j) = 1;
  }
  W_sparse->endFill();
  Id.endFill();

}

solverAC2::~solverAC2()
{

  //delete H_sparse;
  if (primalAvailable)
  {
    delete W_sparse;
    delete M_sparse;
    delete luOfM_sparse;
  }

  delete dfAC_du_sparse;
  delete dfAC_dr_sparse;
  delete dG_dr_sparse;

}

void solverAC2::compute_fAC(bool computeJacobianAsWell)
{

  // update u
  u = (*W_sparse) * VectorXd::Map(r, n * d) + q;

  VectorXd rN(1); // kind of stupid to encapsulate a double into a VectorXd,
  // but it seems necessary since no implicit cast from VectorXd to double is available
  VectorXd rT(d - 1); // d is known at compile-time, we could avoid a malloc here
  VectorXd uN(1);
  VectorXd uT(d - 1); // idem

  MatrixXd PN(1, d);
  MatrixXd PT((int)d - 1, (int)d); // explicit cast needed to avoid ambiguous call to overloaded constructor

  // main loop
  for (unsigned int i = 0; i < n; i++)
  {

    // pass to local coordinates
    PN = E.block(i * d, 0, d, 1).transpose();
    PT = E.block(i * d, 1, d, d - 1).transpose();

    rN = PN * VectorXd::Map(r, n * d).segment(i * d, d);
    rT = PT * VectorXd::Map(r, n * d).segment(i * d, d);
    uN = PN * u.segment(i * d, d);
    uT = PT * u.segment(i * d, d);

    double lambda, ny;
    VectorXd y(d - 1);

    // normal part of fAC
    double x = rN[0] - rhoN[i] * uN[0];

    if (x >= 0)
    {
      fAC[i * d] = x;
      if (computeJacobianAsWell)
      {
        dfAC_du_packed.block(i * d, 0, 1, d) = -rhoN[i] * PN;
        dfAC_dr_packed.block(i * d, 0, 1, d) = PN;
      }
    }
    else
    {
      fAC[i * d] = 0;
      if (computeJacobianAsWell)
      {
        dfAC_du_packed.block(i * d, 0, 1, d) = MatrixXd::Zero(1, d);
        dfAC_dr_packed.block(i * d, 0, 1, d) = MatrixXd::Zero(1, d);
      }
    }

    // tangent part of fAC
    lambda = mu[i] * rN[0];
    y = rT - rhoT[i] * uT;
    ny = y.norm();

    if (ny < lambda)
    {
      fAC.segment(i * d + 1, d - 1) = y;
      if (computeJacobianAsWell)
      {
        dfAC_du_packed.block(i * d + 1, 0, d - 1, d) = -rhoT[i] * PT;
        dfAC_dr_packed.block(i * d + 1, 0, d - 1, d) = PT;
      }
    }
    else if (ny > lambda && lambda > 0)
    {
      fAC.segment(i * d + 1, d - 1) = (lambda / ny) * y;
      if (computeJacobianAsWell)
      {
        dfAC_du_packed.block(i * d + 1, 0, d - 1, d) = // line break
          -rhoT[i] * (lambda / ny) * (MatrixXd::Identity(d - 1, d - 1) - (1 / (ny * ny)) * y * y.transpose()) * PT;
        dfAC_dr_packed.block(i * d + 1, 0, d - 1, d) = // line break
          mu[i] * y / ny * PN + lambda / ny * (MatrixXd::Identity(d - 1, d - 1) - (1 / (ny * ny)) * y * y.transpose()) * PT;
      }
    }
    else
    {
      fAC.segment(i * d + 1, d - 1) = VectorXd::Zero(d - 1);
      if (computeJacobianAsWell)
      {
        dfAC_du_packed.block(i * d + 1, 0, d - 1, d) = MatrixXd::Zero(d - 1, d);
        dfAC_dr_packed.block(i * d + 1, 0, d - 1, d) = MatrixXd::Zero(d - 1, d);
      }
    }

    // return to world coordinates
    // (is there a '*=' operator to shorten the lines below ?)
    fAC.segment(i * d, d) = E.block(i * d, 0, d, d) * fAC.segment(i * d, d);
    if (computeJacobianAsWell)
    {
      dfAC_du_packed.block(i * d, 0, d, d) = E.block(i * d, 0, d, d) * dfAC_du_packed.block(i * d, 0, d, d);
      dfAC_dr_packed.block(i * d, 0, d, d) = E.block(i * d, 0, d, d) * dfAC_dr_packed.block(i * d, 0, d, d);
    }
  }

  // second pass, to assemble dfAC_du and dfAC_dr
  if (computeJacobianAsWell)
  {

    // assemble dfAC_du_sparse and dfAC_dr_sparse
    dfAC_du_sparse->startFill(n * d * d);
    dfAC_dr_sparse->startFill(n * d * d);
    for (unsigned int k = 0; k < n; k++)
      for (int j = 0; j < d; j++)
        for (int i = 0; i < d; i++)
        {
          dfAC_du_sparse->fill(k * d + i, k * d + j) = dfAC_du_packed[k * d + i + (n * d) * j];
          dfAC_dr_sparse->fill(k * d + i, k * d + j) = dfAC_dr_packed[k * d + i + (n * d) * j];
        }
    dfAC_du_sparse->endFill();
    dfAC_dr_sparse->endFill();
  }

  // evaluation of function G
  G = fAC - VectorXd::Map(r, n * d);

  // evaluation of the Jacobian of G
  if (computeJacobianAsWell)
  {
    (*dG_dr_sparse) = (*dfAC_du_sparse) * (*W_sparse) + (*dfAC_dr_sparse) - Id;
  }
}

bool solverAC2::compute_dr()
{
  Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::UmfPack> luOfdG_dr_sparse(*dG_dr_sparse);
  bool flag_sparse = false;
  if (luOfdG_dr_sparse.succeeded())
  {
    G = -G; // change sign of G...
    flag_sparse = luOfdG_dr_sparse.solve(G, &dr);
    G = -G; // ...and set it back
  }
  else
  {
    _LUFailed = true ;
    //std::cout << "sparse LU of dG_dr failed" << std::endl;
  }
  return flag_sparse;
}

double solverAC2::solve()
{
  //return solveDampedNewton();
  //return solveCylqp();
  return solveSoqp();
}

// Newton step + Goldstein-Price line-search
double solverAC2::solveDampedNewton()
{
  compute_fAC(true);
  double q0 = 0.5 * G.dot(G);
  unsigned int iter = 0;
  double t; // to store the step length
  while (iter++ < NIterMax && q0 > Tol)
  {
    // compute Newton step
    bool flag_dr = compute_dr();
    if (!flag_dr)
    {
#ifdef BAVARD
      std::cout << "WARNING: computation of Newton step went wrong";
#endif
    }
    // perform line-search
    bool flag_gp = lineSearchGP(t);
    if (!flag_gp)
    {
      //std::cout << "WARNING: line-search went wrong";
    }
    // update r
    VectorXd::Map(r, n * d) += t * dr;
    // update fAC
    compute_fAC(true);
    q0 = 0.5 * G.dot(G);
  }

  // compute v
  if (primalAvailable)
  {
    VectorXd temp4 = Htrans * VectorXd::Map(r, n * d) - f;
    luOfM_sparse->solve(temp4, &rhs_for_lusolve);
    VectorXd::Map(v, m) = rhs_for_lusolve;
  }
  if (q0 > Tol)
  {
    std::cout << "solverAC2::solveDampedNewton: reached precision " << q0 << " instead of " << Tol << std::endl;
    std::cout << "solverAC2::solveDampedNewton: " << iter << " iterations with NIterMax = " << NIterMax << std::endl;
  }
  std::cout << "solverAC2::solveDampedNewton: reached precision " << q0 << std::endl;
  //std::cout << "r   = " << VectorXd::Map(r,n*d).transpose() << std::endl;
  //std::cout << "u   = " << u.transpose() << std::endl;
  //std::cout << "fAC = " << fAC.transpose() << std::endl;
  //std::cout << "norm(fAC-r) = " << (fAC - VectorXd::Map(r,n*d)).norm() << std::endl;

  return q0;
}

// Newton step
double solverAC2::solvePureNewton()
{
  compute_fAC(true);
  double q0 = 0.5 * G.dot(G);
  unsigned int iter = 0;
  double t; // to store the step length
  while (iter++ < NIterMax && q0 > Tol)
  {
    // compute Newton step
    bool flag_dr = compute_dr();
    if (!flag_dr)
      std::cout << "WARNING: computation of Newton step went wrong";
    // update r (no line-search)
    VectorXd::Map(r, n * d) += dr;
    // update fAC
    compute_fAC(true);
    q0 = 0.5 * G.dot(G);
  }

  // compute v
  if (primalAvailable)
  {
    VectorXd temp4 = Htrans * VectorXd::Map(r, n * d) - f;
    luOfM_sparse->solve(temp4, &rhs_for_lusolve);
    VectorXd::Map(v, m) = rhs_for_lusolve;
  }

  if (q0 > Tol)
    std::cout << "solverAC2::solvePureNewton: reached precision " << q0 << " instead of " << Tol << std::endl;
  return q0;
}

// fixed-point iteration
double solverAC2::solveFixedPoint()
{
  compute_fAC(false);
  double q0 = 0.5 * G.dot(G);
  unsigned int iter = 0;
  double t; // to store the step length
  while (iter++ < NIterMax && q0 > Tol)
  {
    // fixed-point step
    VectorXd::Map(r, n * d) += G;
    // update fAC
    compute_fAC(false);
    q0 = 0.5 * G.dot(G);
  }

  // compute v
  if (primalAvailable)
  {
    VectorXd temp4 = Htrans * VectorXd::Map(r, n * d) - f;
    luOfM_sparse->solve(temp4, &rhs_for_lusolve);
    VectorXd::Map(v, m) = rhs_for_lusolve;
  }
  return q0;
}

// Newton step + GP line-search, unless Newton step computation fails, then fixed-point iter
double solverAC2::solveDampedHybrid()
{
  compute_fAC(true);
  double q0 = 0.5 * G.dot(G);
  unsigned int iter = 0;
  double t; // to store the step length
  while (iter++ < NIterMax && q0 > Tol)
  {
    // compute Newton step
    bool flag_dr = compute_dr();
    if (!flag_dr)
    {
      dr = G;
      t = 1.0;
    }
    else
    {
      // perform line-search
      bool flag_gp = lineSearchGP(t);
      if (!flag_gp)
        std::cout << "WARNING: line-search went wrong";
    }
    // update r
    VectorXd::Map(r, n * d) += t * dr;
    // update fAC
    compute_fAC(true);
    q0 = 0.5 * G.dot(G);
  }

  // compute v
  if (primalAvailable)
  {
    VectorXd temp4 = Htrans * VectorXd::Map(r, n * d) - f;
    luOfM_sparse->solve(temp4, &rhs_for_lusolve);
    VectorXd::Map(v, m) = rhs_for_lusolve;
  }
  if (q0 > Tol)
    std::cout << "solverAC2::solveDampedHybrid: reached precision " << q0 << " instead of " << Tol << std::endl;
  return q0;
}

// Newton step (without LS), unless Newton step computation fails, then fixed-point iter
double solverAC2::solvePureHybrid()
{
  compute_fAC(true);
  double q0 = 0.5 * G.dot(G);
  unsigned int iter = 0;
  double t; // to store the step length
  while (iter++ < NIterMax && q0 > Tol)
  {
    // compute Newton step
    bool flag_dr = compute_dr();
    if (!flag_dr)
    {
      dr = G;
    }
    // update r
    VectorXd::Map(r, n * d) += dr;
    // update fAC
    compute_fAC(true);
    q0 = 0.5 * G.dot(G);
  }

  // compute v
  if (primalAvailable)
  {
    VectorXd temp4 = Htrans * VectorXd::Map(r, n * d) - f;
    luOfM_sparse->solve(temp4, &rhs_for_lusolve);
    VectorXd::Map(v, m) = rhs_for_lusolve;
  }
  if (q0 > Tol)
    std::cout << "solverAC2::solvePureHybrid: reached precision " << q0 << " instead of " << Tol << std::endl;
  return q0;
}


/* Goldstein-Price line-search
   returns a step size t so that q(t) is lower than q(0) for the objective function
   q(t) := 1/2 |f(x+t*dx)|^2
   q(0) and qp(0) are generally known in advance, so we expect them in the input arguments
   t_init provides an initialization for t
*/
bool solverAC2::lineSearchGP(double &t, double t_in)
{

  const double big = 1e5;
  t = t_in;
  // initial "t left" and " right"
  double tl = 0;
  double tr = big;

  // parameters (hardcoded, is it worth letting the user choose ?)
  double m1 = 0.1;
  double m2 = 0.9;

  unsigned int iter = 0;

  // value of q(t) and (dq/dt)(t) at t=0
  double q0 = 0.5 * G.dot(G); // q0 = 1/2 |G(r)|^2
  Eigen::Matrix<double, 1, 1> qp0_mat = G.transpose() * (*dG_dr_sparse) * dr; // qp0 = G'*dG_dr*dr
  double qp0 = qp0_mat[0];

  // line search
  for (iter = 0; iter < NIterMax_LineSearch; iter++)
  {
    // move r to r+t*dr, evaluate fAC and bring r back to its value
    VectorXd::Map(r, n * d) += t * dr;
    compute_fAC(false);
    VectorXd::Map(r, n * d) -= t * dr;
    // evaluate q(t)
    double qt = 0.5 * G.dot(G);
    double slope = (qt - q0) / t;
    //std::cout << "qt = " << qt << " and slope = " << slope << std::endl;
    bool C1 = (slope >= m2 * qp0);
    bool C2 = (slope <= m1 * qp0);
    // three cases
    if (C1 && C2)
    {
      //std::cout << "line-search returned t = " << t << std::endl;
      return true; // ok, success
    }
    else if (!C1)
    {
      //std::cout << "t = " << t << " is too small : slope = " << slope << ", m2*qp0 = " << m2*qp0 << std::endl;
      tl = t;
    }
    else     // not(C2)
    {
      //std::cout << "t = " << t << " is too big : slope = " << slope << ", m1*qp0 = " << m1*qp0 << std::endl;
      tr = t;
    }
    if (tr < big)
    {
      t = 0.5 * (tl + tr);
    }
    else
    {
      t = 10 * t;
    }

  }
#ifdef BAVARD
  std::cout << "GP line-search: reached NITERMAX with t = " << t << std::endl;
#endif
  return false; // failure
}

/* ================================
   "QP over a cylinder" technique (called "method of Haslinger" in my thesis)
   ================================ */

// ===== projCylinder =====
// projection of x onto the cylinder of radius s (in dimension d, n times)
// s should be >= 0 (otherwise behaviour is meaningless and div_by_zero may occur)
void solverAC2::projCylinder(VectorXd &x, const VectorXd &s)
{
  if (x.size() != n * d)
  {
    std::cout << "projK: x should be of size n*d" << std::endl;
    return;
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      // pass to local coordinates
      MatrixXd PN = E.block(i * d, 0, d, 1).transpose();
      MatrixXd PT = E.block(i * d, 1, d, d - 1).transpose();
      VectorXd xN = PN * x.segment(i * d, d);
      VectorXd xT = PT * x.segment(i * d, d);
      VectorXd p(d);

      double nxT = xT.norm();
      if (nxT <= s[i] && xN[0] >= 0)   // x lies inside the cylinder
      {
        p[0] = xN[0];
        p.segment(1, d - 1) = xT;
      }
      else if (nxT > s[i] && xN[0] >= 0)     // x lies inside normal cone
      {
        p[0] = xN[0];
        p.segment(1, d - 1) = s[i] * xT / nxT;
      }
      else     // xN < 0
      {
        p[0] = 0.;
        p.segment(1, d - 1) = s[i] * xT / nxT;
      }
      x.segment(i * d, d) = E.block(i * d, 0, d, d) * p;
    }
  }
}

// ===== cylqp =====
// cylqp (quadratic programming over a cylinder): min 1/2 r'*W*r + q'*r   s.t.   r \in Cyl(s) (cf projCylinder)
// "x" is both init. and return value
void solverAC2::cylqp(VectorXd &x, VectorXd &s)
{

  //TODO critere d'arret indep de step

  double old_val, new_val, ndx, npu;
  static double step; // step length; static may not be clever (store step as a member of the class ?)
  int i, j;
  VectorXd Wx, xtest;
  VectorXd best_iterate = x;

  // params
  int NOuterMax = 50;
  int NInnerMax = 12; //around 10-20 ?
  double initStep = 1.0;
  double npuMax = 1e-6; // caution, this tolerance should be chosen coherently with
  // the outer tolerance (in the fixed point alg.)
  // TODO faire un ndx _relatif_ (genre norm(dx)/n)

  // init
  projCylinder(x, s);
  Wx = (*W_sparse) * x;
  old_val = x.dot(0.5 * Wx + q);
  //std::cout << std::scientific << std::setprecision(8);
  //std::cout << "cylqp: start PG w/ val = " << old_val << std::endl;
  step = initStep;
  for (i = 0; i < NOuterMax; i++)   // outer loop
  {
    u = Wx + q;
    xtest = x - u; // step=1 for stopping criterion
    projCylinder(xtest, s);
    npu = (xtest - x).norm();
    if (npu < npuMax) break;
    j = 0; // counter for inner loop
    if (step <= 50 * initStep)
      step = 5 * step; // try a larger step (optimistic phase)
    //std::cout << "old_val = " << old_val << "; ";
    do
    {
      step = 0.5 * step; // reduce step (realistic phase)
      xtest = x - step * u;
      projCylinder(xtest, s);
      Wx = (*W_sparse) * xtest;
      new_val = xtest.dot(0.5 * Wx + q);
    }
    while (new_val > old_val && ++j <= NInnerMax);
    ndx = (xtest - x).norm();
    //std::cout << "new_val = " << new_val << "; step = " << step << " ;" << j << " loop(s); norm(dx) = " << ndx << std::endl;
    if (new_val < old_val) best_iterate = xtest;
    x = xtest;
    old_val = new_val;
  }
  //std::cout << "cylqp: stop PG w/ val = " << new_val << " (" << i << " iter., npu = " << npu << ", ndx = " << ndx << ", last step = " << step << ")" << std::endl;
  x = best_iterate;
}

// ===== solveCylqp =====
double solverAC2::solveCylqp()
{
  // TODO: setters/getters pour les parametres
  // TODO: last call to cylqp: passer tol << 1 et nitermax grand pour assurer uN>=0
  VectorXd s = VectorXd::Zero(n);
  VectorXd new_s = VectorXd::Zero(n);
  VectorXd rr = VectorXd::Zero(n * d);
  rr = Eigen::Map<VectorXd>(r, n * d);
  double nds;

  // params
  int fpIterMax = 100;
  double ndsMax = 1e-6; //TODO faire un nds relatif (genre norm(ds)/n)

  // eval s=mu*rN
  for (int j = 0; j < n; j++)
  {
    Eigen::Matrix<double, 1, 1> rN = E.block(j * d, 0, d, 1).transpose() * rr.segment(j * d, d);
    s[j] = mu[j] * rN[0];
  }

  // fixed-point loop
  int fpIter = 0;
  std::cout << "solveCylqp: start fp iter." << std::endl;
  do
  {
    cylqp(rr, s);
    for (int j = 0; j < n; j++)
    {
      Eigen::Matrix<double, 1, 1> rN = E.block(j * d, 0, d, 1).transpose() * rr.segment(j * d, d);
      new_s[j] = mu[j] * rN[0];
    }
    //std::cout << "s = " << s.transpose() << std::endl;
    //std::cout << "new_s = " << new_s.transpose() << std::endl;
    nds = (new_s - s).norm();
    s = new_s;
  }
  while (fpIter++ < fpIterMax && nds > ndsMax);
  std::cout << "solveCylqp: stop after " << fpIter << " fp iter. with nds = " << nds << std::endl;
  //std::cout << "rr = " << rr.transpose() << std::endl;


  u = (*W_sparse) * rr + q;
  /*
  std::cout << "uN = ";
  for (int j=0; j<n; j++) {
    Eigen::Matrix<double, 1, 1> uN = E.block(j*d, 0, d, 1).transpose() * u.segment(j*d, d);
    std::cout << uN[0] << "; ";
  }
  std::cout << std::endl;
  */
  // return r
  Eigen::Map<VectorXd>(r, n * d) = rr;

  // compute v
  if (primalAvailable)
  {
    VectorXd temp4 = Htrans * VectorXd::Map(r, n * d) - f;
    luOfM_sparse->solve(temp4, &rhs_for_lusolve);
    VectorXd::Map(v, m) = rhs_for_lusolve;
  }

  compute_fAC(false);

  //std::cout << "r   = " << rr.transpose() << std::endl;
  //std::cout << "u   = " << u.transpose() << std::endl;
  //std::cout << "fAC = " << fAC.transpose() << std::endl;
  //std::cout << "norm(fAC-r) = " << (fAC - VectorXd::Map(r,n*d)).norm() << std::endl;
  //dumpToFile(fAC.norm(), "dumpcyl.txt");

  double sq0 = (fAC - VectorXd::Map(r, n * d)).norm();
  std::cout << "solverAC2::solveCylQP: reached precision " << 0.5 * sq0*sq0 << std::endl;
  return 0.5 * sq0 * sq0;
}

/* ======================================
   "QP over a 2nd order cone" technique (called "our method" in my thesis)
   ====================================== */

// ===== projK =====
// project x onto K_mu (and return the result in x itself)
void solverAC2::projK(VectorXd& x)
{
  if (x.size() != n * d)
  {
    std::cout << "projK: x should be of size n*d" << std::endl;
    return;
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      // pass to local coordinates
      MatrixXd PN = E.block(i * d, 0, d, 1).transpose();
      MatrixXd PT = E.block(i * d, 1, d, d - 1).transpose();
      VectorXd xN = PN * x.segment(i * d, d);
      VectorXd xT = PT * x.segment(i * d, d);
      VectorXd p(d);

      double nxT = xT.norm();
      if (nxT <= mu[i] * xN[0])   // x lies inside the cone
      {
        p[0] = xN[0];
        p.segment(1, d - 1) = xT;
      }
      else if (xN[0] + mu[i] * nxT <= 0)     // x lies inside normal cone
      {
        p.setZero();
      }
      else     // "hard" case
      {
        double alpha = 1 / (1 + mu[i] * mu[i]);
        p[0] = alpha * (xN[0] + mu[i] * nxT);
        p.segment(1, d - 1) = alpha * mu[i] * (mu[i] + xN[0] / nxT) * xT;
      }
      x.segment(i * d, d) = E.block(i * d, 0, d, d) * p;
    }
  }
}

// ===== projKstar ===== (UNUSED)
// project x onto K_mu^star (and return the result in x itself)
void solverAC2::projKstar(VectorXd& x)
{
  if (x.size() != n * d)
  {
    std::cout << "projKstar: x should be of size n*d" << std::endl;
    return;
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      // pass to local coordinates
      MatrixXd PN = E.block(i * d, 0, d, 1).transpose();
      MatrixXd PT = E.block(i * d, 1, d, d - 1).transpose();
      VectorXd xN = PN * x.segment(i * d, d);
      VectorXd xT = PT * x.segment(i * d, d);
      VectorXd p(d);

      double nxT = xT.norm();
      if (mu[i] * nxT <= xN[0])   // x lies inside the cone
      {
        p[0] = xN[0];
        p.segment(1, d - 1) = xT;
      }
      else if (mu[i] * xN[0] + nxT <= 0)     // x lies inside normal cone
      {
        p.setZero();
      }
      else     // "hard" case
      {
        double alpha = 1 / (1 + mu[i] * mu[i]);
        p[0] = mu[i] * alpha * (mu[i] * xN[0] + nxT);
        p.segment(1, d - 1) = alpha * (1 + mu[i] * xN[0] / nxT) * xT;
      }
      x.segment(i * d, d) = E.block(i * d, 0, d, d) * p;
    }
  }
}

// ===== soqp =====
// soqp (2nd order cone quadratic programming): min 1/2 r'*W*r + qs'*r   s.t.   r \in K
// "x" is both init. and return value
void solverAC2::soqp(VectorXd &x, VectorXd &qs)
{
  //TODO critere d'arret indep de step

  double old_val, new_val, ndx, npu;
  static double step; // step length; static may not be clever (store step as a member of the class ?)
  int i, j;
  VectorXd Wx, xtest;
  VectorXd best_iterate = x;

  // params
  int NOuterMax = 50;
  int NInnerMax = 12; //around 10-20 ?
  double initStep = 1.0;
  double npuMax = 1e-6; // caution, this tolerance should be chosen coherently with
  // the outer tolerance (in the fixed point alg.)
  // TODO faire un ndx _relatif_ (genre norm(dx)/n)

  // init
  projK(x);
  Wx = (*W_sparse) * x;
  old_val = x.dot(0.5 * Wx + qs);
  //std::cout << std::scientific << std::setprecision(8);
  //std::cout << "soqp: start PG w/ val = " << old_val << std::endl;

  step = initStep;
  for (i = 0; i < NOuterMax; i++)   // outer loop
  {
    u = Wx + qs;
    xtest = x - u; // step=1 for stopping criterion
    projK(xtest);
    npu = (xtest - x).norm();
    if (npu < npuMax) break;
    j = 0; // counter for inner loop
    if (step <= 50 * initStep)
      step = 5 * step; // try a larger step (optimistic phase)
    //std::cout << "old_val = " << old_val << "; ";
    do
    {
      step = 0.5 * step; // reduce step (realistic phase)
      xtest = x - step * u;
      projK(xtest);
      Wx = (*W_sparse) * xtest;
      new_val = xtest.dot(0.5 * Wx + qs);
    }
    while (new_val > old_val && ++j <= NInnerMax);
    ndx = (xtest - x).norm();

    //std::cout << "new_val = " << new_val << "; step = " << step << " ;" << j << " loop(s); norm(dx) = " << ndx << std::endl;

    if (new_val < old_val) best_iterate = xtest;
    x = xtest;
    old_val = new_val;
  }
  //std::cout << "soqp: stop PG w/ val = " << new_val << " (" << i << " iter., npu = " << npu << ", ndx = " << ndx << ", last step = " << step << ")" << std::endl;

  x = best_iterate;
  //std::cout << "soqp: stop with x = " << x.transpose() << std::endl;
}

// ===== solveSoqp =====
double solverAC2::solveSoqp()
{
  // TODO: setters/getters pour les parametres
  // TODO: last call to cylqp: passer tol << 1 et nitermax grand pour assurer uN>=0
  VectorXd s = VectorXd::Zero(n);
  VectorXd new_s = VectorXd::Zero(n);
  VectorXd qs = VectorXd::Zero(n * d);
  VectorXd rr = VectorXd::Zero(n * d);
  rr = Eigen::Map<VectorXd>(r, n * d);
  double nds;

  // params
  int fpIterMax = 100;
  double ndsMax = 1e-6; //TODO faire un nds relatif (genre norm(ds)/n)

  // init u
  u = (*W_sparse) * rr + q;

  // eval s=mu*|uT| and qs^i = q^i + eN^i*s^i
  qs = q;
  for (int j = 0; j < n; j++)
  {
    Eigen::Matrix<double, 2, 1> uT = E.block(j * d, 1, d, 2).transpose() * u.segment(j * d, d);
    s[j] = mu[j] * uT.norm();
    qs.segment(j * d, d) += E.block(j * d, 0, d, 1) * s[j];
  }

  // fixed-point loop
  int fpIter = 0;
  std::cout << "solveSoqp: start fp iter." << std::endl;
  do
  {
    soqp(rr, qs);
    u = (*W_sparse) * rr + q;
    qs = q;
    for (int j = 0; j < n; j++)
    {
      Eigen::Matrix<double, 2, 1> uT = E.block(j * d, 1, d, 2).transpose() * u.segment(j * d, d);
      new_s[j] = mu[j] * uT.norm();
      qs.segment(j * d, d) += E.block(j * d, 0, d, 1) * new_s[j];
    }

    //std::cout << "s = " << s.transpose() << std::endl;
    //std::cout << "new_s = " << new_s.transpose() << std::endl;
    nds = (new_s - s).norm();
    s = new_s;
  }
  while (++fpIter < fpIterMax && nds > ndsMax);
  std::cout << "solveSoqp: stop after " << fpIter << " fp iter. with nds = " << nds << std::endl;
  //std::cout << "rr = " << rr.transpose() << std::endl;

  u = (*W_sparse) * rr + q;

  /*
  std::cout << "uN = ";
  for (int j=0; j<n; j++) {
    Eigen::Matrix<double, 1, 1> uN = E.block(j*d, 0, d, 1).transpose() * u.segment(j*d, d);
    std::cout << uN[0] << "; ";
  }
  std::cout << std::endl;
  */
  // return r
  Eigen::Map<VectorXd>(r, n * d) = rr;

  // compute v
  if (primalAvailable)
  {
    VectorXd temp4 = Htrans * VectorXd::Map(r, n * d) - f;
    luOfM_sparse->solve(temp4, &rhs_for_lusolve);
    VectorXd::Map(v, m) = rhs_for_lusolve;
  }

  compute_fAC(false);

  //std::cout << "r   = " << rr.transpose() << std::endl;
  //std::cout << "u   = " << u.transpose() << std::endl;
  //std::cout << "fAC = " << fAC.transpose() << std::endl;
  //std::cout << "norm(fAC-r) = " << (fAC - VectorXd::Map(r,n*d)).norm() << std::endl;
  //dumpToFile(fAC.norm(), "dumpcyl.txt");

  double sq0 = (fAC - VectorXd::Map(r, n * d)).norm();
  std::cout << "solverAC2::solveSoqp: reached precision " << 0.5 * sq0*sq0 << std::endl;
  return 0.5 * sq0 * sq0;
}

// ===== file input / output =====
int solverAC2::dumpToFile(double res, const char *filename)
{
  std::ofstream dumpfile;
  dumpfile.open(filename);
  if (!(dumpfile.is_open()))   // failed to open file
  {
    std::cout << "dumpToFile: failed to open file.\n";
    return 1;
  }
  dumpfile << std::setprecision(9);
  dumpfile << "tol = " << Tol << ";\n";
  dumpfile << "res = " << res << ";\n";
  dumpfile << "data.d = " << d << ";\n";
  dumpfile << "data.m = " << m << ";\n";
  dumpfile << "data.n = " << n << ";\n";

  printSparse(*M_sparse, dumpfile);

  dumpfile << "data.f = [\n";
  dumpfile << f;
  dumpfile << "\n];\n";

  dumpfile << "data.H = [\n";
  dumpfile << H;
  dumpfile << "\n];\n";

  dumpfile << "data.w = [\n";
  dumpfile << w;
  dumpfile << "\n];\n";

  dumpfile << "data.E = [\n";
  dumpfile << E;
  dumpfile << "\n];\n";

  dumpfile << "data.mu = [\n";
  dumpfile << mu;
  dumpfile << "\n];\n";

  dumpfile.close();
  return 0;
}

int solverAC2::readFromFile(const char *filename)
{

  std::ifstream dumpfile;
  dumpfile.open(filename);
  if (!(dumpfile.is_open()))   // failed to open file
  {
    std::cout << "readFromFile: failed to open file.\n";
    return 1;
  }
  char c[12]; // pas genial, voir le mecanisme avec std:getline(stream, string) ci-apres

  // tol
  dumpfile.read(c, 6);
  c[6] = '\0';
  std::string stol("tol = ");
  if (stol.compare(c))
  {
    std::cout << "warning: readFromFile: could not read tol" << std::endl;
    return 1;
  }
  else
  {
    dumpfile >> Tol;
  }
  //std::cout << "tol = " << Tol << std::endl; //dbg
  dumpfile.ignore(2); // jump the ';' and the newline

  // res
  double res;
  dumpfile.read(c, 6);
  c[6] = '\0';
  std::string sres("res = ");
  if (sres.compare(c))
  {
    std::cout << "warning: readFromFile: could not read res" << std::endl;
    return 1;
  }
  else
  {
    dumpfile >> res;
  }
  //std::cout << "res = " << res << std::endl; //dbg
  dumpfile.ignore(2); // jump the ';\n'

  // d
  unsigned int d_in;
  dumpfile.read(c, 9);
  c[9] = '\0';
  std::string sd("data.d = ");
  if (sd.compare(c))
  {
    std::cout << "warning: readFromFile: could not read data.d" << std::endl;
    return 1;
  }
  else
  {
    dumpfile >> d_in;
    if (d_in != d)
    {
      std::cout << "warning: readFromFile: d should be equal to " << d << " (not " << d_in << ")" << std::endl;
    }
  }
  //std::cout << "d = " << d_in << std::endl; //dbg
  dumpfile.ignore(2); // jump the ';\n'

  // m
  unsigned int m_in;
  dumpfile.read(c, 9);
  c[9] = '\0';
  std::string sm("data.m = ");
  if (sm.compare(c))
  {
    std::cout << "warning: readFromFile: could not read data.m" << std::endl;
    return 1;
  }
  else
  {
    dumpfile >> m_in;
    if (m_in != m)
    {
      std::cout << "warning: readFromFile: m should be equal to " << m << " (not " << m_in << ")" << std::endl;
    }
  }
  //std::cout << "m = " << m_in << std::endl; //dbg
  dumpfile.ignore(2); // jump the ';\n'

  // n
  unsigned int n_in;
  dumpfile.read(c, 9);
  c[9] = '\0';
  std::string sn("data.n = ");
  if (sn.compare(c))
  {
    std::cout << "warning: readFromFile: could not read data.n" << std::endl;
    return 1;
  }
  else
  {
    dumpfile >> n_in;
    if (n_in != n)
    {
      std::cout << "warning: readFromFile: n should be equal to " << n << " (not " << n_in << ")" << std::endl;
    }
  }
  //std::cout << "n = " << n_in << std::endl; //dbg
  dumpfile.ignore(2); // jump the ';\n'

  // M_sparse
  std::string str; //TODO: utiliser str partout
  unsigned int nnz;
  std::streampos beg_ij;
  std::getline(dumpfile, str);
  if (str.compare("ij = ["))
  {
    std::cout << "warning: readFromFile: could not read ij" << std::endl;
    std::cout << "found:" << std::endl;
    std::cout << str << std::endl;
    return 1;
  }
  else
  {
    beg_ij = dumpfile.tellg();
    nnz = 0; // count the number of nonzeros
    do
    {
      if (dumpfile.good())
      {
        std::getline(dumpfile, str);
        //std::cout << str << std::endl;
        nnz++;
      }
      else
      {
        std::cout << "warning: readFromFile: failure while reading ij" << std::endl;
        return 1;
      }
    }
    while (str.compare("];"));
  }
  nnz--; // do not count the last line ("];")

  // go back to ij = ...
  dumpfile.seekg(beg_ij);
  Eigen::MatrixXi ij(nnz, 2);
  int k = 0;
  for (int i = 0; i < nnz; i++)
  {
    if (dumpfile.good())
    {
      dumpfile >> ij(k, 0);
      dumpfile.ignore(2); // jump ", "
      dumpfile >> ij(k, 1);
      dumpfile.ignore(1); // jump '\n'
      //std::cout << i_in << ", " << j_in << std::endl;
      k++;
    }
    else
    {
      std::cout << "warning: readFromFile: failure while reading ij" << std::endl;
      return 1;
    }
  }
  dumpfile.ignore(3); // "];\n"
  //std::cout << ij << std::endl;

  // val
  std::getline(dumpfile, str);
  if (str.compare("val = ["))
  {
    std::cout << "warning: readFromFile: could not read val" << std::endl;
    std::cout << "found:" << std::endl;
    std::cout << str << std::endl;
    return 1;
  }
  else
  {
    M_sparse->startFill(nnz);
    double val;
    for (int k = 0; k < nnz; k++)
    {
      dumpfile >> val;
      dumpfile.ignore(1); // '\n'
      M_sparse->fill(ij(k, 0) - 1, ij(k, 1) - 1) = val; // caution: "-1" since we count from 0
    }
    M_sparse->endFill();
  }
  dumpfile.ignore(3); // "];\n"
  //printSparse(*M_sparse, std::cout);

  // mn
  std::getline(dumpfile, str);
  std::stringstream s_mn;
  s_mn << "mn = [" << m << ", " << m << "];";
  if (str.compare(s_mn.str()))
  {
    std::cout << "warning: readFromFile: could not read mn" << std::endl;
    std::cout << "found:" << std::endl;
    std::cout << str << std::endl;
    std::cout << "instead of:" << std::endl;
    std::cout << s_mn.str() << std::endl;
    return 1;
  }

  // data.f
  getline(dumpfile, str);
  if (str.compare("data.f = ["))
  {
    std::cout << "warning: readFromFile: could not read data.f" << std::endl;
    std::cout << "found:" << std::endl;
    std::cout << str << std::endl;
    return 1;
  }
  else
  {
    for (int i = 0; i < m; i++)
    {
      if (dumpfile.good())
      {
        dumpfile >> f(i);
        dumpfile.ignore(1); // '\n'
      }
      else
      {
        std::cout << "warning: readFromFile: failure while reading data.f" << std::endl;
        return 1;
      }
    }
  }
  //std::cout << f << std::endl; //dbg
  dumpfile.ignore(3); // jump the '];\n'

  // data.H
  getline(dumpfile, str);
  //std::cout << str << std::endl;
  if (str.compare("data.H = ["))
  {
    std::cout << "warning: readFromFile: could not read data.H" << std::endl;
    std::cout << "found:" << std::endl;
    std::cout << str << std::endl;
    return 1;
  }
  else
  {
    for (int i = 0; i < n * d; i++)
    {
      for (int j = 0; j < m; j++)
        dumpfile >> H(i, j);
      dumpfile.ignore(1); // '\n'
    }
  }
  //std::cout << H << std::endl;
  dumpfile.ignore(3); // jump the '];\n'

  // data.w
  getline(dumpfile, str);
  //std::cout << str << std::endl;
  if (str.compare("data.w = ["))
  {
    std::cout << "warning: readFromFile: could not read data.w" << std::endl;
    std::cout << "found:" << std::endl;
    std::cout << str << std::endl;
    return 1;
  }
  else
  {
    for (int i = 0; i < n * d; i++)
    {
      if (dumpfile.good())
      {
        dumpfile >> w(i);
        dumpfile.ignore(1); // '\n'
      }
      else
      {
        std::cout << "warning: readFromFile: failure while reading data.w" << std::endl;
        return 1;
      }
    }
  }
  //std::cout << w << std::endl; //dbg
  dumpfile.ignore(3); // jump the '];\n'

  getline(dumpfile, str);
  //std::cout << str << std::endl;
  if (str.compare("data.E = ["))
  {
    std::cout << "warning: readFromFile: could not read data.E" << std::endl;
    std::cout << "found:" << std::endl;
    std::cout << str << std::endl;
    return 1;
  }
  else
  {
    for (int i = 0; i < n * d; i++)
    {
      for (int j = 0; j < 3; j++)
        dumpfile >> E(i, j);
      dumpfile.ignore(1); // '\n'
    }
  }
  //std::cout << E << std::endl;
  dumpfile.ignore(3); // jump the '];\n'

  // data.mu
  getline(dumpfile, str);
  if (str.compare("data.mu = ["))
  {
    std::cout << "warning: readFromFile: could not read data.mu" << std::endl;
    std::cout << "found:" << std::endl;
    std::cout << str << std::endl;
    return 1;
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      if (dumpfile.good())
      {
        dumpfile >> mu(i);
        dumpfile.ignore(1); // '\n'
      }
      else
      {
        std::cout << "warning: readFromFile: failure while reading data.mu" << std::endl;
        return 1;
      }
    }
  }
  //std::cout << mu << std::endl; //dbg
  dumpfile.ignore(3); // jump the '];\n'

  dumpfile.close();
}

// modify r so that all normal velocities are nonnegative
// returns the total violation of this constraint afterwards
// this is reaaaally sloppy (and should probably NOT be used)
double solverAC2::getRidOfPenetration()
{

  MatrixXd E2 = MatrixXd::Zero(n, n * d); // sooooo sparse :(
  for (int i = 0; i < n; i++)
  {
    E2.block(i, i * d, 1, d) = E.block(i * d, 0, d, 1).transpose();
  }
  //std::cout << "E=\n" << E << std::endl;
  //std::cout << "E2=\n" << E2 << std::endl;

  MatrixXd A = MatrixXd::Zero(n, n * d);
  VectorXd b = VectorXd::Zero(n);
  A = -E2 * (*W_sparse);
  b = E2 * q;

  // normalize components of A,b
  for (int i = 0; i < n; i++)
  {
    double a = A.block(i, 0, 1, n * d).norm();
    A.block(i, 0, 1, n * d) = A.block(i, 0, 1, n * d) / a;
    b(i) = b(i) / a;
  }

  // compute modified r
  int nmax = 500;
  int iter;
  double tol = 1e-8;
  double p;
  VectorXd gp(n * d);
  VectorXd slack(n);
  for (iter = 0; iter < nmax; iter++)
  {
    // clean p and gp
    p = 0;
    gp = VectorXd::Zero(n * d);
    slack = A * VectorXd::Map(r, n * d) - b;
    //std::cout << "slack=\n" << slack << std::endl;

    // compute p and gp
    for (int i = 0; i < n; i++)
    {
      if (slack(i) >= 0)
      {
        p += 0.5 * slack(i) * slack(i);
        gp += slack(i) * A.block(i, 0, 1, n * d).transpose();
      }
    }
    //std::cout << "p=" << p << std::endl;
    if (p < tol) break;
    VectorXd::Map(r, n * d) -= 0.25 * gp; // step length -> adjust or choose by line search
  }

  //std::cout << "bidouille: tol=" << tol << ", reached p=" << p << " in " << iter << " iterations (max " << nmax << ")\n";

  // update u
  u = (*W_sparse) * VectorXd::Map(r, n * d) + q;
  VectorXd uNs = E2 * u;
  //std::cout << "uNs=\n" << uNs << std::endl;

  // compute v
  if (primalAvailable)
  {
    Eigen::MatrixXd temp6 = Htrans * VectorXd::Map(r, n * d) - f;
    luOfM_sparse->solve(temp6, &rhs_for_lusolve);
    VectorXd::Map(v, m) = rhs_for_lusolve;
  }
  return p;
}

// trick of the geek: could use std::copy into cout to avoid the loop
template <typename T>
void solverAC2::printVector(const std::vector<T>& vec, const char * name)
{
  int l = vec.size();
  std::cout << name << " = ";
  for (int i = 0; i < l; i++)
    std::cout << vec[i] << "; ";
  std::cout << std::endl;
}

void solverAC2::printSparse(Eigen::SparseMatrix<double> & M, std::ostream& out)
{

  // iterate over a sparse matrix : output in scilab format
  // ij
  out << "ij = [" << std::endl;
  for (int k = 0; k < M.outerSize(); ++k)
  {
    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(M, k); it; ++it)
    {
      //cout << "val = " << it.value() << ", row = " <<  it.row() << ", col = " << it.col() << ", idx = " << it.index() << endl;
      out <<  it.row() + 1 << ", " << it.col() + 1 << std::endl;
    }
  }
  out << "];" << std::endl;
  // val
  out << "val = [" << std::endl;
  for (int k = 0; k < M.outerSize(); ++k)
  {
    for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(M, k); it; ++it)
    {
      //cout << "val = " << it.value() << ", row = " <<  it.row() << ", col = " << it.col() << ", idx = " << it.index() << endl;
      out <<  it.value() << std::endl;
    }
  }
  out << "];" << std::endl;
  // mn
  out << "mn = [" << M.rows() << ", " << M.cols() << "];" << std::endl;

}

