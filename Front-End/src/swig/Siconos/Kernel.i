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

// SWIG interface for Siconos Kernel
%module(directors="1") Kernel

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


%include "KernelTypes.i"

// 2. try to hide SP::Type on python side

%define SP_TYPE(TYPE)
%shared_ptr(TYPE);
%enddef



SP_TYPE(NonSmoothLaw);
SP_TYPE(NewtonImpactNSL);
SP_TYPE(NewtonImpactFrictionNSL);

SP_TYPE(NonSmoothDynamicalSystem);
SP_TYPE(Topology);

SP_TYPE(DynamicalSystem);
SP_TYPE(LagrangianDS);
SP_TYPE(LagrangianLinearTIDS);

SP_TYPE(Relation);
SP_TYPE(UnitaryRelation);
SP_TYPE(LagrangianR);
SP_TYPE(LagrangianLinearTIR);
SP_TYPE(LagrangianRheonomousR);
SP_TYPE(LagrangianScleronomousR);

SP_TYPE(Interaction)

SP_TYPE(TimeDiscretisation);

SP_TYPE(OneStepNSProblem);
SP_TYPE(LinearOSNS);
SP_TYPE(LCP);
SP_TYPE(FrictionContact);

SP_TYPE(OneStepIntegrator);
SP_TYPE(Moreau);

SP_TYPE(Simulation);
SP_TYPE(TimeStepping);

SP_TYPE(Model);

SP_TYPE(CircularDS);
SP_TYPE(Disk);
SP_TYPE(Circle);
SP_TYPE(ExternalBody);
SP_TYPE(DiskDiskR);
SP_TYPE(DiskPlanR);
SP_TYPE(DiskMovingPlanR);
SP_TYPE(SiconosBodies);

SP_TYPE(SpaceFilter);

SP_TYPE(SiconosMatrix);
SP_TYPE(SimpleMatrix);
SP_TYPE(SiconosVector);
SP_TYPE(SimpleVector);
SP_TYPE(BlockVector);

SP_TYPE(InteractionsSet)


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
%ignore operator+;
%ignore operator[];
%ignore operator-;
%ignore operator*;
%ignore operator/;
%ignore operator==;
%ignore operator=;

// defined in SimpleMatrix.cpp 
%ignore private_addprod;
%ignore private_prod;
%ignore prod;
%ignore axpy_prod;
%ignore subprod;
%ignore axpy_prod;
%ignore gemv;
%ignore gemm;

%ignore getWMap;
%ignore getWBoundaryConditionsMap;
%ignore getDSBlocks;
%ignore ACCEPT_SP_VISITORS;


%include "SiconosPointers.hpp"
%import "SimulationTypeDef.hpp" 

%include "InteractionsSet.hpp"
%include "SiconosSet.hpp"

%import "SiconosGraph.hpp"

%import "boost/config.hpp"
%import "boost/graph/graph_utility.hpp"

// boost >= 1.40
%import "boost/version.hpp"
#if (BOOST_VERSION >= 104000)
%import "boost/smart_ptr/enable_shared_from_this.hpp"
#else
%import "boost/enable_shared_from_this.hpp"
#endif



%include "Tools.hpp"
%include "addons.hpp"

%include "SiconosAlgebra.hpp"
%include "SiconosVector.hpp"
%include "BlockVector.hpp"
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



%fragment("StdSequenceTraits");

%fragment("StdMapTraits");


%template (InteractionsSet) SiconosSet<Interaction,double*>;

%template (sharedModel) boost::enable_shared_from_this<Model>;

%template (dsv) std::vector<boost::shared_ptr<DynamicalSystem> >;

%template (dsi) std::pair<unsigned int, unsigned int >;

%template (dsp) std::pair<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<DynamicalSystem> >;

%template (dspv) std::vector<std::pair<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<DynamicalSystem> > >;

%template (dsiv) std::vector<std::pair<unsigned int, unsigned int > >;


