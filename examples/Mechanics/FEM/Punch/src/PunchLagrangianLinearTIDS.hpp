/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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



#include "SiconosKernel.hpp"

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;



class LinearElacticMaterial
{
  double _E;
  double _nu;
  double _rho;
  double _lambda;
  double _mu;

public:
  LinearElacticMaterial(){};

  
  LinearElacticMaterial(double E, double nu, double rho)
  {
    _E = E;
    _nu = nu;
    _rho = rho;

    _lambda = _E * _nu /( (1.0 +_nu) * (1.0 - 2* _nu));
    _mu = _E /(2.0*(1+_nu));
  };

  double mu()
  {
    return _mu;
  };
  double lambda()
  {
    return _lambda;
  };

};


class PunchLagrangianLinearTIDS : public LagrangianLinearTIDS
{

  LinearElacticMaterial _mat;
  
  SP::SimpleMatrix create_matrix_from_mfem(SparseMatrix A)
  {
    int * Ai  = A.GetI();
    int * Aj  = A.GetJ();
    double * Ax  = A.GetData();
    int size = A.Size();

    int nnz = Ai[size];
    
    SP::SimpleMatrix M(new SimpleMatrix(size,size,Siconos::SPARSE,nnz));
   
    for (int row =0; row < size ; row++)
    {
      for (int k = Ai[row], end = Ai[row+1]; k < end; k++)
      {
        M->setValue(row, Aj[k], Ax[k]);
      }
    }

    //M->display();
    return M;
  };

public:
  PunchLagrangianLinearTIDS( const char *mesh_file, LinearElacticMaterial mat)
  {
    _mat =mat;

    int order = 1;
    bool static_cond = false;



   
    // 2. Read the mesh from the given mesh file. We can handle triangular,
    //    quadrilateral, tetrahedral or hexahedral elements with the same code.
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    if (mesh->bdr_attributes.Max() < 2)
    {
      cerr << "\nInput mesh should have at least two materials and "
           << "two boundary attributes! (See schematic in ex2.cpp)\n"
           << endl;
    }

       // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 5,000
   //    elements.
   // {
   //    int ref_levels =
   //       (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);
   //    for (int l = 0; l < ref_levels; l++)
   //    {
   //       mesh->UniformRefinement();
   //    }
   // }
       // 5. Define a finite element space on the mesh. Here we use vector finite
   //    elements, i.e. dim copies of a scalar finite element space. The vector
   //    dimension is specified by the last argument of the FiniteElementSpace
   //    constructor. For NURBS meshes, we use the (degree elevated) NURBS space
   //    associated with the mesh nodes.
   FiniteElementCollection *fec;
   FiniteElementSpace *fespace;
   if (mesh->NURBSext)
   {
      fec = NULL;
      fespace = mesh->GetNodes()->FESpace();
   }
   else
   {
      fec = new H1_FECollection(order, dim);
      fespace = new FiniteElementSpace(mesh, fec, dim);
   }


   unsigned int nDof = fespace->GetTrueVSize();
   _ndof = nDof;

   cout << "Number of finite element unknowns: " << fespace->GetTrueVSize()
        << endl << "Assembling: " << flush;
    // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking only
   //    boundary attribute 1 from the mesh as essential and converting it to a
   //    list of true dofs.
   Array<int> ess_tdof_list, ess_bdr(mesh->bdr_attributes.Max());
   ess_bdr = 0;
   ess_bdr[0] = 1;
   fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   std::cout  << "ess_tdof_list.Print()" << std::endl;
   ess_tdof_list.Print();
   std::cout  << "ess_bdr.Print()" << std::endl;
   ess_bdr.Print();
   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system. In this case, b_i equals the boundary integral
   //    of f*phi_i where f represents a "pull down" force on the Neumann part
   //    of the boundary and phi_i are the basis functions in the finite element
   //    fespace. The force is defined by the VectorArrayCoefficient object f,
   //    which is a vector of Coefficient objects. The fact that f is non-zero
   //    on boundary attribute 2 is indicated by the use of piece-wise constants
   //    coefficient for its last component.
   VectorArrayCoefficient f(dim);
   for (int i = 0; i < dim-1; i++)
   {
      f.Set(i, new ConstantCoefficient(0.0));
   }
   {
      Vector pull_force(mesh->bdr_attributes.Max());
      pull_force = 0.0;
      pull_force(1) = -1.0e-2;
      f.Set(dim-1, new PWConstCoefficient(pull_force));
   }

   LinearForm *b = new LinearForm(fespace);
   b->AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(f));
   cout << "r.h.s. ... " << flush;
   b->Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.0;

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the linear elasticity integrator with piece-wise
   //    constants coefficient lambda and mu.
   Vector lambda(mesh->attributes.Max());
   lambda = 1.0;
   lambda(0) = lambda(1)*50;
   PWConstCoefficient lambda_func(lambda);
   Vector mu(mesh->attributes.Max());
   mu = 1.0;
   mu(0) = mu(1)*50;
   PWConstCoefficient mu_func(mu);

   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func,mu_func));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   cout << "matrix ... " << flush;
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   cout << "done." << endl;

   //A.Print();
   cout << "Size of linear system: " << A.Height() << endl;



   BilinearForm *m = new BilinearForm(fespace);
   m->AddDomainIntegrator(new VectorMassIntegrator());
   m->Assemble();



   // shift the eigenvalue corresponding to eliminated dofs to a large value
   // m->EliminateEssentialBCDiag(ess_bdr, numeric_limits<double>::min());
   m->Finalize();

   SparseMatrix M;
   m->FormLinearSystem(ess_tdof_list, x, *b, M, X, B);
   //M.Print();

   double position_init=0.0;
   double velocity_init=-1.0;

   // -- Initial positions and velocities --
   SP::SiconosVector q0(new SiconosVector(nDof,position_init));
   SP::SiconosVector v0(new SiconosVector(nDof,velocity_init));

   _init(q0,v0);
   _K = create_matrix_from_mfem(A);
   _mass = create_matrix_from_mfem(M);


   
  };





};

TYPEDEF_SPTR(PunchLagrangianLinearTIDS)
