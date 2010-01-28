#ifndef _SOLVERAC2_H
/* solverAC2
A class to solve the linear incremental friction problem.
F. Cadoux Apr. 20th, 2009 (INRIA).
*/

/* Usage:
we want to solve the problem
/
| M * v + f = H^T * r
| u = H * v + w
| (u, v) satisfy Coulomb's law
\
(cf F. Cadoux, "An optimization-based algorithm for Coulomb frictional contact", Esaim-Procs)
where
- (v,u,r) are the unknown
 * v : generalized velocities
 * u : relative velocities at contact points
 * r : contact impulses
- M, f, H, w as well as well as the friction coefficients mu and the transition matrix to the local base E are input data. They are specified by the input args of the constructor of solverAC. Note that the convention is that the normal (given in E) goes from B to A, and that if one of the colliding objects does not belong to the simulation, it must be body B. The matrix HA[i] (resp. HB[i]) gives the velocity of A (resp. B) at the i-th contact point.
*/

// preciser : (faire une doc)
// - que objA et ObjB comptent à partir de 0
// - que la normale va de B vers A
// - que HA et HB sont donnés sans le signe "-" (uA = HA*vA+wA, uB=HB*vB+wB)
// - les options de compilation conseillées par les auteurs d'Eigen

#include<vector>
#include<iostream>

#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/Sparse>

USING_PART_OF_NAMESPACE_EIGEN

class solverAC2
{
public:
  bool _LUFailed ;

  solverAC2(
    unsigned int d_in, // dimension (2 or 3)
    unsigned int NObj, // number of subsystems
    const unsigned int * ndof, // array of size NObj, the number of degree of freedom of each subsystem
    const double *const * MassMat, // array of pointers to the mass matrix of each subsystem
    const double * f_in, // the constant term in M*v+f=H^T*r
    unsigned int n_in, // number of contact points
    const double * mu_in, // array of size n giving the friction coeffs
    const double * E_in, // array of size n*d*d giving the n normals followed by the n tangent vectors (and by again n tangent vectors if d==3). Said otherwise, E is a (n*d) x d matrix, stored column-major, formed by n blocks of size d x d with each block being an orthogonal matrix (the transition matrix from the world space coordinates (x1, x2, x3) to the local coordinates (xN, xT1, xT2)
    const double * w_in, // array of size n*d, the constant term in u=H*v+w
    // Caution: ObjA and ObjB below are counted from ZERO !
    const int * ObjA, // array of size n, the first object involved in the i-th contact (must be an internal object)
    const int * const ObjB, // array of size n, the first object involved in the i-th contact (-1 for an external object)
    const double *const HA[], // array of size n, containing pointers to a dense, colum-major matrix of size d*ndof[ObjA[i]] corresponding to the H-matrix of objA (NULL for an external object)
    const double *const HB[], // array of size n, containing pointers to a dense, colum-major matrix of size d*ndof[ObjA[i]] corresponding to the H-matrix of objB (could be 'unsigned' but we keep an int for consistency with ObjA)
    double r0[], // length n*d : initialization for r (in world space coordinates) + used to return computed r
    double v0[] // length m: to return computed v
  );

  // constructor from the so-called "dual" formulation, i.e. the user inputs W,q directly (and not M,H,w,f)
  // here we assume that W is provided as an Eigen sparse matrix
  solverAC2(
    unsigned int d_in,
    unsigned int n_in,
    Eigen::SparseMatrix<double> * W_sparse_in,
    const double * q_in,
    const double * mu_in,
    const double * E_in,
    double *r0);

  // constructor from file
  solverAC2(
    unsigned int d_in, // dimension (2 or 3)
    unsigned int n_in, // number of contact points
    unsigned int m_in,  // number of degrees of freedom
    const char* filename, // file
    double r0[], // length n*d : initialization for r (in world space coordinates) + used to return computed r
    double v0[] // length m: to return computed v
  );

  ~solverAC2();

  // setters / getters
  unsigned int getNIterMax()
  {
    return NIterMax;
  };
  unsigned int getNIterMax_LineSearch()
  {
    return NIterMax_LineSearch;
  };
  double getTol()
  {
    return Tol;
  };
  void setNIterMax(unsigned int N)
  {
    NIterMax = N;
  };
  void setNIterMax_LineSearch(unsigned int N)
  {
    NIterMax_LineSearch = N;
  };
  void setTol(double t)
  {
    Tol = t;
  };

  // solve functions
  double solve();  // legacy: calls solveDampedNewton
  double solveDampedNewton();  // Newton step + Goldstein-Price line-search
  double solvePureNewton();  // Newton step
  double solveFixedPoint();  // fixed-point iteration
  double solveDampedHybrid();  // Newton step + GP line-search, unless Newton step computation fails, then fixed-point iter
  double solvePureHybrid();  // Newton step (without LS), unless Newton step computation fails, then fixed-point iter

  // bidouille affreuse
  double getRidOfPenetration();

  // print content
  //void print_M () {std::cout << "M =\n" << M.format(printFormat) << std::endl;} ;
  void print_f()
  {
    std::cout << "f =\n" << f.format(printFormat) << std::endl;
  } ;
  void print_H()
  {
    std::cout << "H =\n" << H.format(printFormat) << std::endl;
  } ;
  void print_w()
  {
    std::cout << "w =\n" << w.format(printFormat) << std::endl;
  };
  void print_E()
  {
    std::cout << "E =\n" << E.format(printFormat) << std::endl;
  };
  void print_mu()
  {
    std::cout << "mu =\n" << mu.format(printFormat) << std::endl;
  }; // autant stocker mu comme un VectorXd

  // file input/output functions
  int dumpToFile(double res, const char *filename);
  int readFromFile(const char *filename);
  Eigen::IOFormat printFormat; // stores the format for printing matrices (you can modify its public attributes)

protected:

  // data
  const unsigned int d;
  // m should be 'const' but is it possible ? (m results from sum(ndof))
  // trick : declare m as a reference to a const _m; _m will be mutable
  // and m will be const, if we use only m we should avoid to mutate _m

  unsigned int _m;
  const unsigned int & m;
  const unsigned int n;

  //const double *const mu;
  VectorXd mu;

  Eigen::SparseMatrix<double> * M_sparse; // mass matrix
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::UmfPack> * luOfM_sparse; // should be a block-Cholesky
  Eigen::SparseMatrix<double> * W_sparse;// "Delassus" operator

  VectorXd f;
  MatrixXd H; // keep this until lu/sparse accepts sparse mat as rhs (then use: Eigen::SparseMatrix<double> * H_sparse;)
  MatrixXd Htrans;
  VectorXd w;
  MatrixXd E;
  VectorXd q;
  Eigen::SparseMatrix<double> Id; // identity matrix

  bool primalAvailable; // do we now M,H,f,w (primal) or only W,q (dual) ?

  // current iterate
  double *v; // we keep v and r as double* rather than VectorXd to manipulate
  // directly the data v0 and r0 given by the user (is it clever ?)
  double *r;
  VectorXd u;
  VectorXd dr;
  VectorXd rhs_for_lusolve; // not very useful, could be removed

  // params rhoN and rhoT for AC function
  std::vector<double> rhoN;
  std::vector<double> rhoT;

  // AC function (storage)
  VectorXd fAC;
  MatrixXd dfAC_du_packed;
  MatrixXd dfAC_dr_packed;
  Eigen::SparseMatrix<double> * dfAC_du_sparse;
  Eigen::SparseMatrix<double> * dfAC_dr_sparse;
  VectorXd G;
  Eigen::SparseMatrix<double> * dG_dr_sparse;

  // AC function (computation)
  void compute_fAC(bool computeJacobianAsWell = false);
  bool compute_dr();  // compute a Newton step
  bool lineSearchGP(double &t_out, double t_in = 1);  // compute a step size which decreases norm(G)

  // params
  unsigned int NIterMax;
  unsigned int NIterMax_LineSearch;
  double Tol;

  // helper functions
  template<typename T>
  void printVector(const std::vector<T>& vec, const char * name);
  void printSparse(Eigen::SparseMatrix<double> & M, std::ostream& out);

public: // TODO repasser protected
  // projectors
  void projK(VectorXd& x);  // onto K_mu
  void projKstar(VectorXd& x);  // onto K_mu^star = -K_{1/mu}
  void projCylinder(VectorXd &x, const VectorXd &s);
  // soqp (2nd order cone quadratic programming): min 1/2 r'*W*r + q'*r   s.t.   r \in K
  void soqp(VectorXd &x, VectorXd &qs); // x: init, and return value
  // cylqp (quadratic programming over a cylinder): min 1/2 r'*W*r + q'*r   s.t.   r \in Cyl(s) (cf projCylinder)
  void cylqp(VectorXd &x,  VectorXd &s); // x: init, and return value
  double solveCylqp();
  double solveSoqp();
};
#endif //_SOLVERAC2_H

