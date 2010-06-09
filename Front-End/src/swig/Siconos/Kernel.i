// -*- C++ -*-
// Siconos-Front-End version 3.2.0, Copyright INRIA 2005-2010.
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.	
// Siconos is a free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// Siconos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Siconos; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr 
//	

// Siconos.i - SWIG interface for Siconos
%module Kernel

 // where is doxygen feature ?
%feature("autodoc", "1");

%{
#define SWIG_FILE_WITH_INIT
#include <sstream>
#if defined(Py_COMPLEXOBJECT_H)
#undef c_sum
#undef c_diff
#undef c_neg
#undef c_prod
#undef c_quot
#undef c_pow
#undef c_abs
#endif
#include "SiconosKernel.hpp"
#include "Disk.hpp"
#include "Circle.hpp"
#include "DiskDiskR.hpp"
#include "DiskPlanR.hpp"
#include "DiskMovingPlanR.hpp"
#include "CircleCircleR.hpp"
#include "ExternalBody.hpp"
#include "SpaceFilter.hpp"
#include "SiconosBodies.hpp"
#include "addons.hpp"
%} 

// mandatory !
%rename (lambda_) lambda;

// shared ptr management
%include "boost_shared_ptr.i"

// numpy macros
%include numpy.i 	

%init %{
  import_array();
%}

// handle standard exceptions
%include "exception.i"
%exception
{
  try
  {
    $action
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch (const std::out_of_range& e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
}

// handle stl data types
%include "stl.i"

// 1. Vector and Matrix <=> numpy array (dense only)


// numpy array to SP::SimpleVector (here a SiconosVector is always a
// SimpleVector)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosVector>)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in,fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> (PyArrayObject* array=NULL, int is_new_object)
{
  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,1) ||
      !require_native(array) || !require_contiguous(array)) SWIG_fail;

  SP::SimpleVector tmp;
  tmp.reset(new SimpleVector(array_size(array,0)));
  // copy : with SimpleVector based on resizable std::vector there is
  // no other way
  memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
  $1 = tmp;
 }

// numpy array to SP::SimpleMatrix (here a SiconosMatrix is always a
// SimpleMatrix)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosMatrix>)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in) boost::shared_ptr<SiconosMatrix> (PyArrayObject* array=NULL, int is_new_object) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,2) ||
      !require_native(array) || !require_contiguous(array)) SWIG_fail;

  SP::SimpleMatrix tmp;
  tmp.reset(new SimpleMatrix(array_size(array,0), array_size(array,1)));
  // copy : with SimpleMatrix based on resizable std::vector
  memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
  $1 = tmp;
 }


// from C++ to python 
%template() boost::shared_ptr<SiconosVector>;
%template() boost::shared_ptr<SimpleVector>;
%template() boost::shared_ptr<SiconosMatrix>;
%template() boost::shared_ptr<SimpleMatrix>;

%typemap(out) boost::shared_ptr<SiconosVector>
{
  npy_intp this_vector_dim[1];
  this_vector_dim[0]=$1->size();
  $result = PyArray_SimpleNewFromData(1,this_vector_dim,NPY_DOUBLE,$1->getArray());
}

%typemap(out) boost::shared_ptr<SiconosMatrix>
{
  npy_intp this_matrix_dim[2];
  this_matrix_dim[0]=$1->size(0);
  this_matrix_dim[1]=$1->size(1);
  $result = PyArray_SimpleNewFromData(2,this_matrix_dim,NPY_DOUBLE,$1->getArray());
  PyArray_UpdateFlags((PyArrayObject *)$result, NPY_FORTRAN);
}



// SiconosMatrix in api mean SimpleMatrix here
%apply (boost::shared_ptr<SiconosVector>) { (SP::SimpleVector) };
%apply (boost::shared_ptr<SiconosVector>) { (boost::shared_ptr<SimpleVector>) };
%apply (boost::shared_ptr<SiconosVector>) { (SP::SiconosVector) };

%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SimpleMatrix) };
%apply (boost::shared_ptr<SiconosMatrix>) { (boost::shared_ptr<SimpleMatrix>) };
%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SiconosMatrix) }; 


// 2. try to hide SP::Type on python side

%define SP_TYPE(TYPE, BASE)
  
DEFINE_SPTR(TYPE);

typedef SP::##TYPE XSPtr##TYPE;

SWIG_SHARED_PTR_DERIVED(XSPtr##TYPE, BASE, TYPE);

%enddef

%include "SiconosPointers.hpp"

%define REGISTER(X,B)
SP_TYPE(X,B)
%include("X##.hpp")
%enddef

SP_TYPE(NonSmoothLaw,NonSmoothLaw);
SP_TYPE(NewtonImpactNSL,NonSmoothLaw);
SP_TYPE(NewtonImpactFrictionNSL,NonSmoothLaw);



SP_TYPE(NonSmoothDynamicalSystem,NonSmoothDynamicalSystem);
SP_TYPE(Topology,Topology);

SP_TYPE(DynamicalSystem,DynamicalSystem);
SP_TYPE(LagrangianDS,DynamicalSystem);
SP_TYPE(LagrangianLinearTIDS,LagrangianDS);

SP_TYPE(Relation,Relation);
SP_TYPE(UnitaryRelation,UnitaryRelation);
SP_TYPE(LagrangianR,Relation);
SP_TYPE(LagrangianLinearTIR,LagrangianR);
SP_TYPE(LagrangianRheonomousR,LagrangianR);
SP_TYPE(LagrangianScleronomousR,LagrangianR);

SP_TYPE(Interaction,Interaction)

SP_TYPE(TimeDiscretisation,TimeDiscretisation);

SP_TYPE(OneStepNSProblem,OneStepNSProblem);
SP_TYPE(LinearOSNS, OneStepNSProblem);
SP_TYPE(LCP,LinearOSNS);
SP_TYPE(FrictionContact,LinearOSNS);

SP_TYPE(OneStepIntegrator, OneStepIntegrator);
SP_TYPE(Moreau, OneStepIntegrator);

SP_TYPE(Simulation,Simulation);
SP_TYPE(TimeStepping, Simulation);

SP_TYPE(Model, Model);

SP_TYPE(CircularDS,LagrangianDS);
SP_TYPE(Disk,CircularDS);
SP_TYPE(Circle,CircularDS);
SP_TYPE(ExternalBody,LagrangianDS);
SP_TYPE(DiskDiskR,LagrangianScleronomousR);
SP_TYPE(DiskPlanR,LagrangianScleronomousR);
SP_TYPE(DiskMovingPlanR,LagrangianRheonomousR);
SP_TYPE(SiconosBodies,SiconosBodies);

SP_TYPE(SpaceFilter,SpaceFilter)
SP_TYPE(Hashed,Hashed)
SP_TYPE(HashedDisk,Hashed)
SP_TYPE(HashedCircle,Hashed)
SP_TYPE(HashedSphereLDS,Hashed)

SP_TYPE(SiconosMatrix,SiconosMatrix)
SP_TYPE(SimpleMatrix,SiconosMatrix)
SP_TYPE(SiconosVector,SiconosVector)
SP_TYPE(SimpleVector,SiconosVector)

%include "InteractionsSet.hpp"
%include "SiconosSet.hpp"

%template (InteractionsSet) SiconosSet<Interaction,double*>;

// boost >= 1.40
%import "boost/smart_ptr/enable_shared_from_this.hpp"
%template (sharedModel) boost::enable_shared_from_this<Model>;

SP_TYPE(InteractionsSet,InteractionsSet)


// dummy namespaces to make swig happy
namespace boost 
{
  namespace numeric
  {
    namespace ublas
    {
    }
    namespace bindings
    {
      namespace atlas
      {
      }
    }
  }
};


namespace SP
{
};

// note Visitor etc. => module Visitor, Type if ever needed

%ignore PURE_DEF;
%ignore VIRTUAL_ACCEPT_VISITORS;
%ignore ACCEPT_STD_VISITORS;
%ignore ACCEPT_VISITORS;

%ignore nullDeleter;

// defined in SimpleVector.cpp
%ignore setBlock;
%ignore add;
%ignore sub;
%ignore axpby;
%ignore axpy;
%ignore inner_prod;
%ignore outer_prod;
%ignore scal;
%ignore subscal;
%ignore cross_product;
%ignore operator +;
%ignore operator [];
%ignore operator -;
%ignore operator *;
%ignore operator /;
%ignore operator ==;
%ignore operator =;

// defined in SimpleMatrix.cpp 
%ignore private_addprod;
%ignore private_prod;
%ignore prod;
%ignore axpy_prod;
%ignore subprod;
%ignore axpy_prod;
%ignore gemv;
%ignore gemm;

%include "Tools.hpp"

DEFINE_SPTR(SiconosVector);
DEFINE_SPTR(SiconosMatrix);

%include "addons.hpp"


%include "SiconosVector.hpp"
%include "SiconosMatrix.hpp"
%include "SimpleVector.hpp"
%include "SimpleMatrix.hpp"
%include "DynamicalSystem.hpp"
%include "NonSmoothDynamicalSystem.hpp"
%include "Topology.hpp"
%include "LagrangianDS.hpp"
%include "LagrangianLinearTIDS.hpp"

%include "Relation.hpp"
%include "UnitaryRelation.hpp"
%include "LagrangianR.hpp"
%include "LagrangianLinearTIR.hpp"
%include "LagrangianScleronomousR.hpp"
%include "LagrangianRheonomousR.hpp"
%include "Interaction.hpp"
%include "Model.hpp"

%include "OneStepIntegrator.hpp"
%include "Moreau.hpp"

%include "OneStepNSProblem.hpp"
%include "LinearOSNS.hpp"
%include "LCP.hpp"
%include "FrictionContact.hpp"

%include "Simulation.hpp"
%include "TimeStepping.hpp"

%include "NonSmoothLaw.hpp"
%include "NewtonImpactNSL.hpp"
%include "NewtonImpactFrictionNSL.hpp"

%include "TimeDiscretisation.hpp"

%include "CircularDS.hpp"
%include "Disk.hpp"
%include "Circle.hpp"
%include "DiskDiskR.hpp"
%include "DiskPlanR.hpp"
%include "DiskMovingPlanR.hpp"
%include "CircleCircleR.hpp"
%include "ExternalBody.hpp"
%include "SpaceFilter.hpp"
%include "SiconosBodies.hpp"
%import "SiconosGraph.hpp"

%import "boost/config.hpp"
%import "boost/graph/graph_utility.hpp"

 // boost >= 1.40
 //%import "boost/graph/adjacency_list.hpp"

%fragment("StdSequenceTraits");

%fragment("StdMapTraits");

// this fail with swig::type_name
//%template (VMap) std::map<boost::shared_ptr<DynamicalSystem>,
//                          graph_traits<adjacency_list<
//                                         listS, listS, undirectedS, 
//                                         property <vertex_bundle_t, boost::shared_ptr<DynamicalSystem>, 
//                                                   property < vertex_color_t , 
//                                                              default_color_type , 
//                                                              property < vertex_index_t, size_t > > >, 
//                                         property <edge_bundle_t, boost::shared_ptr<UnitaryRelation>, 
//                                                   property < edge_color_t , 
//                                                              default_color_type , 
//                                                              property < edge_index_t, size_t > > > > >::vertex_descriptor >;
//

%include "SimulationTypeDef.hpp"

%template (dsv) std::vector<boost::shared_ptr<DynamicalSystem> >;

%template (dsi) std::pair<unsigned int, unsigned int >;

%template (dsp) std::pair<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<DynamicalSystem> >;

%template (dspv) std::vector<std::pair<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<DynamicalSystem> > >;

%template (dsiv) std::vector<std::pair<unsigned int, unsigned int > >;

%template (dsg) SiconosGraph<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<UnitaryRelation> >;


//%template (urv) std::vector<boost::shared_ptr<UnitaryRelation> >;
