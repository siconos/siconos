/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "LagrangianDS.h"
#include "LagrangianDSXML.h"
#include "BlockVector.h"
#include "BlockMatrix.h"

using namespace std;

// Private function to set linked with members of Dynamical top class
void LagrangianDS::connectToDS()
{
  // dim
  n = 2 * ndof;

  // All links between DS and LagrangianDS class members are pointer links, which means
  // that no useless memory is allocated when connection is established.
  // One exception: zero and identity matrices, used to filled in M and jacobianXF.

  // Initial conditions

#ifndef WithSmartPtr
  x0 = new BlockVector(q0, velocity0);
  isAllocatedIn["x0"] = true;
#else
  x0.reset(new BlockVector(q0, velocity0));
#endif

  // Current State: \f$ x \f$ and rhs = \f$ \dot x \f$

#ifndef WithSmartPtr
  x[0] = new BlockVector(q[0], q[1]);
  isAllocatedIn["x"] = true;
  x[1] = new BlockVector(q[1], q[2]);
  isAllocatedIn["rhs"] = true;
#else
  x[0].reset(new BlockVector(q[0], q[1]));
  x[1].reset(new BlockVector(q[1], q[2]));
#endif

  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related functions.
}


LagrangianDS::LagrangianDS(int newNumber, const SiconosVector& newQ0, const SiconosVector& newVelocity0):
  DynamicalSystem(LNLDS, newNumber, 2 * newQ0.size()), ndof(newQ0.size()),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL)
{
  // -- Memory allocation for vector and matrix members --
  // Initial conditions
  q0.reset(new SimpleVector(newQ0));
  velocity0.reset(new SimpleVector(newVelocity0));

  // Current state
  q.resize(3);
  q[0].reset(new SimpleVector(*q0));
  q[1].reset(new SimpleVector(*velocity0));
  q[2].reset(new SimpleVector(ndof));

  // set allocation flags: true for required input, false for others
  initPluginFlags(false);
  jacobianFInt.resize(2);
  jacobianNNL.resize(2);
  p.resize(3);
  computeJacobianFIntPtr.resize(2);
  computeJacobianNNLPtr.resize(2);
}

void LagrangianDS::initPluginFlags(bool val)
{
  isPlugin["mass"] = val;
  isPlugin["fInt"] = val;
  isPlugin["fExt"] = val;
  isPlugin["NNL"] = val;
  isPlugin["jacobianFInt0"] = val;
  isPlugin["jacobianFInt1"] = val;
  isPlugin["jacobianNNL0"] = val;
  isPlugin["jacobianNNL1"] = val;
}
// -- Default constructor --
LagrangianDS::LagrangianDS():
  DynamicalSystem(LNLDS), ndof(0), computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL)
{
  // Protected constructor - Only call from derived class(es).

  initPluginFlags(false);

#ifndef WithSmartPtr
  initAllocationFlags(false);
  q.resize(3, NULL);
  p.resize(3, NULL);
#else
  q.resize(3);
  p.resize(3);
#endif


  // !!! No plug-in connection !!!
}

// --- Constructor from an xml file ---
#ifndef WithSmartPtr
LagrangianDS::LagrangianDS(DynamicalSystemXML * dsXML, SP::NonSmoothDynamicalSystem newNsds):
  DynamicalSystem(dsXML, newNsds), ndof(0), q0(NULL), velocity0(NULL), qMemory(NULL), velocityMemory(NULL), mass(NULL), fInt(NULL), fExt(NULL),
  NNL(NULL), fL(NULL), computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL)
#else
LagrangianDS::LagrangianDS(DynamicalSystemXMLSPtr dsXML, SP::NonSmoothDynamicalSystem newNsds):
  DynamicalSystem(LNLDS), ndof(0), computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL)
#endif
{
  // -- Lagrangian  xml object --
  LagrangianDSXMLSPtr lgptr = boost::static_pointer_cast <LagrangianDSXML>(dsxml);

  // === Initial conditions ===
  // Warning: ndof is given by q0.size() !!
  if (! lgptr->hasQ0())
    RuntimeException::selfThrow("LagrangianDS:: xml constructor, q0 is a required input");

#ifndef WithSmartPtr
  q0 = new SimpleVector(lgptr->getQ0()); // required node in xml file
#else
  q0.reset(new SimpleVector(lgptr->getQ0())); // required node in xml file
#endif

  ndof = q0->size();

  if (! lgptr->hasVelocity0())
    RuntimeException::selfThrow("LagrangianDS:: xml constructor, v0 is a required input");

#ifndef WithSmartPtr
  velocity0 = new SimpleVector(lgptr->getVelocity0()); // required node in xml file
#else
  velocity0.reset(new SimpleVector(lgptr->getVelocity0())); // required node in xml file
#endif

  if (velocity0->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS::xml constructor - size of input velocity0 differs from ndof");


  // --- Current State (optional input) ---

#ifndef WithSmartPtr
  q.resize(3, NULL);

  if (lgptr->hasQ())
    q[0] = new SimpleVector(lgptr->getQ()); // Means that first q is different from q0 ?? Strange case ...
  else
    q[0] = new SimpleVector(*q0);           // initialize q with q0
  if (lgptr->hasVelocity())
    q[1] = new SimpleVector(lgptr->getVelocity0()); // same remark as for q
  else
    q[1] = new SimpleVector(*velocity0);

  q[2] = new SimpleVector(ndof);

  p.resize(3, NULL);

  initAllocationFlags(); // Default configuration
#else
  q.resize(3);

  if (lgptr->hasQ())
    q[0].reset(new SimpleVector(lgptr->getQ())); // Means that first q is different from q0 ?? Strange case ...
  else
    q[0].reset(new SimpleVector(*q0));           // initialize q with q0
  if (lgptr->hasVelocity())
    q[1].reset(new SimpleVector(lgptr->getVelocity0())); // same remark as for q
  else
    q[1].reset(new SimpleVector(*velocity0));

  q[2].reset(new SimpleVector(ndof));

  p.resize(3);
#endif

  // Memories
  if (lgptr->hasQMemory())   // qMemory
  {

#ifndef WithSmartPtr
    qMemory = new SiconosMemory(lgptr->getQMemoryXML());
    isAllocatedIn["qMemory"] = true;
#else
    qMemory.reset(new SiconosMemory(lgptr->getQMemoryXML()));
#endif
  }


  if (lgptr->hasVelocityMemory())   // velocityMemory
  {

#ifndef WithSmartPtr
    velocityMemory = new SiconosMemory(lgptr->getVelocityMemoryXML());
    isAllocatedIn["velocityMemory"] = true;
#else
    velocityMemory.reset(new SiconosMemory(lgptr->getVelocityMemoryXML()));
#endif

  }

  initPluginFlags(false);

  string plugin;
  // mass
  if (lgptr->isMassPlugin()) // if mass is plugged
  {
    plugin = lgptr->getMassPlugin();
    setComputeMassFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
  }
#ifndef WithSmartPtr
  else mass = new SimpleMatrix(lgptr->getMassMatrix());
#else
  else mass.reset(new SimpleMatrix(lgptr->getMassMatrix()));
#endif

  // === Optional inputs ===

#ifndef WithSmartPtr
  jacobianFInt.resize(2, NULL);
  jacobianNNL.resize(2, NULL);
#else
  jacobianFInt.resize(2);
  jacobianNNL.resize(2);
#endif

  // fInt
  if (lgptr->hasFInt())  // if fInt is given
  {
    if (lgptr->isFIntPlugin()) // if fInt is plugged
    {
      plugin = lgptr->getFIntPlugin();
      setComputeFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {

#ifndef WithSmartPtr
      fInt = new SimpleVector(lgptr->getFIntVector());
      isAllocatedIn["fInt"] = true;
#else
      fInt.reset(new SimpleVector(lgptr->getFIntVector()));
#endif

    }
#ifndef WithSmartPtr
    computeJacobianFIntPtr.resize(2, NULL);
#else
    computeJacobianFIntPtr.resize(2);
#endif

  }

  // fExt
  if (lgptr->hasFExt())  // if fExt is given
  {
    if (lgptr->isFExtPlugin())// if fExt is plugged
    {
      plugin = lgptr->getFExtPlugin();
      setComputeFExtFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {

#ifndef WithSmartPtr
      fExt = new SimpleVector(lgptr->getFExtVector());
      isAllocatedIn["fExt"] = true;
#else
      fExt.reset(new SimpleVector(lgptr->getFExtVector()));
#endif

    }
  }

  // NNL
  if (lgptr ->hasNNL())// if NNL is given
  {
    if (lgptr->isNNLPlugin())// if NNL is plugged
    {
      plugin = lgptr->getNNLPlugin();
      setComputeNNLFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {

#ifndef WithSmartPtr
      NNL = new SimpleVector(lgptr->getNNLVector());
      isAllocatedIn["NNL"] = true;
#else
      NNL.reset(new SimpleVector(lgptr->getNNLVector()));
#endif

    }

#ifndef WithSmartPtr
    computeJacobianNNLPtr.resize(2, NULL);
#else
    computeJacobianNNLPtr.resize(2);
#endif

  }

  for (unsigned int i = 0; i < 2; ++i)
  {
    // jacobian(s) of fInt
    string name = "jacobianFInt";
    if (lgptr ->hasJacobianFInt(i))// if given
    {
      name += toString<unsigned int>(i);
      if (lgptr->isJacobianFIntPlugin(i))// if is plugged
      {
        plugin = lgptr->getJacobianFIntPlugin(i);
        setComputeJacobianFIntFunction(i, cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else
      {
#ifndef WithSmartPtr
        jacobianFInt[i] = new SimpleMatrix(lgptr->getJacobianFIntMatrix(i));
        isAllocatedIn[name] = true;
#else
        jacobianFInt[i].reset(new SimpleMatrix(lgptr->getJacobianFIntMatrix(i)));
#endif

      }
    }

    // jacobian of NNL
    name = "jacobianNNL";
    if (lgptr -> hasJacobianNNL(i)) // if given
    {
      name += toString<unsigned int>(i);
      if (lgptr->isJacobianNNLPlugin(i))// if is plugged
      {
        plugin = lgptr->getJacobianNNLPlugin(i);
        setComputeJacobianNNLFunction(i, cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else
      {

#ifndef WithSmartPtr
        jacobianNNL[i] = new SimpleMatrix(lgptr->getJacobianNNLMatrix(i));
        isAllocatedIn[name] = true;
#else
        jacobianNNL[i].reset(new SimpleMatrix(lgptr->getJacobianNNLMatrix(i)));
#endif

      }
    }
  }
}

// From a set of data; Mass filled-in directly from a siconosMatrix -
// This constructor leads to the minimum Lagrangian System form: \f$ M\ddot q = p \f$
LagrangianDS::LagrangianDS(int newNumber, const SiconosVector& newQ0, const SiconosVector& newVelocity0, const SiconosMatrix& newMass):
  DynamicalSystem(LNLDS, newNumber, 2 * newQ0.size()), ndof(newQ0.size()),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL)
{
  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

#ifndef WithSmartPtr
  // Mass matrix
  mass = new SimpleMatrix(newMass);

  // Initial conditions
  q0 = new SimpleVector(newQ0);
  velocity0 = new SimpleVector(newVelocity0);

  // Current state
  q.resize(3, NULL);
  q[0] = new SimpleVector(*q0);
  q[1] = new SimpleVector(*velocity0);
  q[2] = new SimpleVector(ndof);

  // set allocation flags: true for required input, false for others
  initAllocationFlags(); // Default
  initPluginFlags(false);

  jacobianFInt.resize(2, NULL);
  jacobianNNL.resize(2, NULL);
  p.resize(3, NULL);
  computeJacobianFIntPtr.resize(2, NULL);
  computeJacobianNNLPtr.resize(2, NULL);
#else
  // Mass matrix
  mass.reset(new SimpleMatrix(newMass));

  // Initial conditions
  q0.reset(new SimpleVector(newQ0));
  velocity0.reset(new SimpleVector(newVelocity0));

  // Current state
  q.resize(3);
  q[0].reset(new SimpleVector(*q0));
  q[1].reset(new SimpleVector(*velocity0));
  q[2].reset(new SimpleVector(ndof));

  // set allocation flags: true for required input, false for others
  initPluginFlags(false);

  jacobianFInt.resize(2);
  jacobianNNL.resize(2);
  p.resize(3);
  computeJacobianFIntPtr.resize(2);
  computeJacobianNNLPtr.resize(2);
#endif

}

// From a set of data - Mass loaded from a plugin
// This constructor leads to the minimum Lagrangian System form: \f$ M(q)\ddot q = p \f$
#ifndef WithSmartPtr
LagrangianDS::LagrangianDS(int newNumber, const SiconosVector& newQ0, const SiconosVector& newVelocity0, const string& massName):
  DynamicalSystem(LNLDS, newNumber, 2 * newQ0.size()), ndof(newQ0.size()), q0(NULL), velocity0(NULL), qMemory(NULL), velocityMemory(NULL), mass(NULL),
  fInt(NULL), fExt(NULL), NNL(NULL), fL(NULL), computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL)
#else
LagrangianDS::LagrangianDS(int newNumber, const SiconosVector& newQ0, const SiconosVector& newVelocity0, const string& massName):
  DynamicalSystem(LNLDS), ndof(0), computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeNNLPtr(NULL)
#endif
{
#ifndef WithSmartPtr
  // Initial conditions
  q0 = new SimpleVector(newQ0);
  velocity0 = new SimpleVector(newVelocity0);

  // Current state
  q.resize(3, NULL);
  q[0] = new SimpleVector(*q0);
  q[1] = new SimpleVector(*velocity0);
  q[2] = new SimpleVector(ndof);

  // Mass
  initPluginFlags(false);
  setComputeMassFunction(cShared.getPluginName(massName), cShared.getPluginFunctionName(massName));

  // set allocation flags: true for required input, false for others
  initAllocationFlags();// default

  jacobianFInt.resize(2, NULL);
  jacobianNNL.resize(2, NULL);
  p.resize(3, NULL);
  computeJacobianFIntPtr.resize(2, NULL);
  computeJacobianNNLPtr.resize(2, NULL);
#else
  // Initial conditions
  q0.reset(new SimpleVector(newQ0));
  velocity0.reset(new SimpleVector(newVelocity0));

  // Current state
  q.resize(3);
  q[0].reset(new SimpleVector(*q0));
  q[1].reset(new SimpleVector(*velocity0));
  q[2].reset(new SimpleVector(ndof));

  // Mass
  initPluginFlags(false);
  setComputeMassFunction(cShared.getPluginName(massName), cShared.getPluginFunctionName(massName));

  jacobianFInt.resize(2);
  jacobianNNL.resize(2);
  p.resize(3);
  computeJacobianFIntPtr.resize(2);
  computeJacobianNNLPtr.resize(2);
#endif

}

// Destructor
LagrangianDS::~LagrangianDS()
{
#ifndef WithSmartPtr
  if (isAllocatedIn["q"]) delete q[0];
  q[0] = NULL;
  if (isAllocatedIn["q0"]) delete q0;
  q0 = NULL;
  if (isAllocatedIn["qMemory"]) delete qMemory;
  qMemory = NULL;
  if (isAllocatedIn["velocity"]) delete q[1] ;
  q[1] = NULL;
  if (isAllocatedIn["velocity0"])delete velocity0 ;
  velocity0 = NULL;
  if (isAllocatedIn["velocityMemory"])delete velocityMemory;
  velocityMemory = NULL;
  if (isAllocatedIn["acceleration"]) delete q[2] ;
  q[2] = NULL;
  if (isAllocatedIn["p0"]) delete p[0] ;
  p[0] = NULL;
  if (isAllocatedIn["p1"]) delete p[1] ;
  p[1] = NULL;
  if (isAllocatedIn["p2"]) delete p[2] ;
  p[2] = NULL;
  if (isAllocatedIn["mass"]) delete mass;
  mass = NULL;
  if (isAllocatedIn["fInt"])delete fInt ;
  fInt = NULL;
  if (isAllocatedIn["fExt"])delete fExt ;
  fExt = NULL;
  if (isAllocatedIn["NNL"])delete NNL ;
  NNL = NULL;
  if (isAllocatedIn["fL"])delete fL ;
  fL = NULL;
  if (isAllocatedIn["jacobianFL0"]) delete jacobianFL[0]  ;
  if (isAllocatedIn["jacobianFL1"]) delete jacobianFL[1]  ;
  if (isAllocatedIn["jacobianFInt0"]) delete jacobianFInt[0] ;
  if (isAllocatedIn["jacobianFInt1"]) delete jacobianFInt[1] ;
  if (isAllocatedIn["jacobianNNL0"]) delete jacobianNNL[0] ;
  if (isAllocatedIn["jacobianNNL1"]) delete jacobianNNL[1] ;
#endif

  jacobianFL.clear();
  jacobianFInt.clear();
  jacobianNNL.clear();

#ifndef WithSmartPtr
  if (isAllocatedIn["jXb10"]) delete workMatrix["jacobianXBloc10"];
  workMatrix["jacobianXBloc10"] = NULL;
  if (isAllocatedIn["jXb11"]) delete workMatrix["jacobianXBloc11"];
  workMatrix["jacobianXBloc11"] = NULL;
#endif

}

bool LagrangianDS::checkDynamicalSystem()
{
  bool output = true;
  // ndof
  if (ndof == 0)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }

  // q0 and velocity0 != NULL
  if (! q0 || ! velocity0)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }

  // Mass
  if (! mass)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - Mass not set.");
    output = false;
  }

  // fInt
  if ((fInt && isPlugin["fInt"]) && (! jacobianFInt[0] || ! jacobianFInt[1]))
    // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined fInt but not its Jacobian (according to q and velocity).");
    output = false;
  }

  // NNL
  if ((NNL  && isPlugin["NNL"]) && (! jacobianNNL[0] || ! jacobianNNL[1]))
    // ie if NNL is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined NNL but not its Jacobian (according to q and velocity).");
    output = false;
  }

  if (isPlugin["f"] || isPlugin["jacobianXF"])
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - f and/or its Jacobian can not be plugged for a Lagrangian system.");
    output = false;
  }
  return output;
}

// TEMPORARY FUNCTION: Must be called before this->initialize
void LagrangianDS::initP(const string& simulationType)
{
  if (simulationType == "TimeStepping")
  {
#ifndef WithSmartPtr
    p[1] = new SimpleVector(ndof);
    isAllocatedIn["p1"] = true;
#else
    p[1].reset(new SimpleVector(ndof));
#endif

    p[2] = p[1];
    p[0] = p[1];

#ifndef WithSmartPtr
    isAllocatedIn["p2"] = false;
#endif

  }
  else if (simulationType == "EventDriven")
  {
#ifndef WithSmartPtr
    p[1] = new SimpleVector(ndof);
    isAllocatedIn["p1"] = true;
    p[2] = new SimpleVector(ndof);
    isAllocatedIn["p2"] = true;
#else
    p[1].reset(new SimpleVector(ndof));
    p[2].reset(new SimpleVector(ndof));
#endif



  }
}

void LagrangianDS::initFL()
{
  //   unsigned int count = 0;
  //   if(fInt!=NULL)
  //     ++count;
  //   if(NNL!=NULL)
  //     ++count;
  //   if(fExt!=NULL)
  //     ++count;
  // First case: fL is equal to a sum of members and needs memory allocation.
  //if (count > 1)
  //{ TEMP : find a way to avoid copy when count == 1. Pb: fL = - NNL or FInt

#ifndef WithSmartPtr
  fL = new SimpleVector(ndof);
  isAllocatedIn["fL"] = true;
#else
  fL.reset(new SimpleVector(ndof));
#endif

  //}
  //  else if (count == 1) // fL = fInt or fExt or NNL; no memory allocation, just a pointer link.
  //{
  //       if(fInt!=NULL)
  //  fL=fInt;
  //       else if(NNL!=NULL)
  //  fL=NNL;
  //       else if(fExt!=NULL)
  //  fL=fExt;
  //}
  // else nothing! (count = 0 and fL need not to be allocated).

  // === Now the jacobians ===
#ifndef WithSmartPtr
  jacobianFL.resize(2, NULL);
#else
  jacobianFL.resize(2);
#endif

  for (unsigned int i = 0; i < jacobianFInt.size(); ++i)
  {
    //       count = 0;
    //       if(jacobianFInt[i]!=NULL)
    //  ++count;
    //       if(jacobianNNL[i]!=NULL)
    //  ++count;
    //       // First case: jacobianFL is equal to a sum of members and needs memory allocation.
    //       if (count > 1)
    //  {
#ifndef WithSmartPtr
    string name = "jacobianFL" + toString<unsigned int>(i);
    jacobianFL[i] = new SimpleMatrix(ndof, ndof);
    isAllocatedIn[name] = true;
#else
    jacobianFL[i].reset(new SimpleMatrix(ndof, ndof));
#endif

    //     }
    //       else if (count == 1)
    //  {
    //    if(jacobianFInt[i]!=NULL)
    //      jacobianFL[i] = jacobianFInt[i];
    //    else if(jacobianNNL[i]!=NULL)
    //      jacobianFL[i] = jacobianNNL[i];
    //  }
    //       // else nothing!
  }
}

void LagrangianDS::initRhs(double time)
{

  // Solve Mq[2]=fL+p.
  *q[2] = *(p[2]); // Warning: r/p update is done in Interactions/Relations

  if (fL)
  {
    computeFL(time);
    *q[2] += *fL;
  }
  computeMass();
  // Copy of Mass into workMatrix for LU-factorization.

#ifndef WithSmartPtr
  workMatrix["invMass"] = new SimpleMatrix(*mass);
#else
  workMatrix["invMass"].reset(new SimpleMatrix(*mass));
#endif

  workMatrix["invMass"]->PLUForwardBackwardInPlace(*q[2]);

  bool flag1 = false, flag2 = false;
  if (jacobianFL[0])
  {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianFL(0, time);

#ifndef WithSmartPtr
    workMatrix["jacobianXBloc10"] = new SimpleMatrix(*jacobianFL[0]);
    isAllocatedIn["jXb10"] = true;
#else
    workMatrix["jacobianXBloc10"].reset(new SimpleMatrix(*jacobianFL[0]));
#endif

    workMatrix["invMass"]->PLUForwardBackwardInPlace(*workMatrix["jacobianXBloc10"]);
    flag1 = true;
  }

  if (jacobianFL[1])
  {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianFL(1, time);

#ifndef WithSmartPtr
    workMatrix["jacobianXBloc11"] = new SimpleMatrix(*jacobianFL[1]);
    isAllocatedIn["jXb11"] = true;
#else
    workMatrix["jacobianXBloc11"].reset(new SimpleMatrix(*jacobianFL[1]));
#endif

    workMatrix["invMass"]->PLUForwardBackwardInPlace(*workMatrix["jacobianXBloc11"]);
    flag2 = true;
  }

#ifndef WithSmartPtr
  workMatrix["zero-matrix"] = new SimpleMatrix(ndof, ndof, ZERO);
  workMatrix["Id-matrix"] = new SimpleMatrix(ndof, ndof, IDENTITY);

  if (flag1 && flag2)
    jacobianXRhs = new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["jacobianXBloc10"], workMatrix["jacobianXBloc11"]);
  else if (flag1) // flag2 = false
    jacobianXRhs = new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["jacobianXBloc10"], workMatrix["zero-matrix"]);
  else if (flag2) // flag1 = false
    jacobianXRhs = new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["zero-matrix"], workMatrix["jacobianXBloc11"]);
  else
    jacobianXRhs = new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["zero-matrix"], workMatrix["zero-matrix"]);
  isAllocatedIn["jacobianXRhs"] = true;
#else
  workMatrix["zero-matrix"].reset(new SimpleMatrix(ndof, ndof, ZERO));
  workMatrix["Id-matrix"].reset(new SimpleMatrix(ndof, ndof, IDENTITY));

  if (flag1 && flag2)
    jacobianXRhs.reset(new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["jacobianXBloc10"], workMatrix["jacobianXBloc11"]));
  else if (flag1) // flag2 = false
    jacobianXRhs.reset(new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["jacobianXBloc10"], workMatrix["zero-matrix"]));
  else if (flag2) // flag1 = false
    jacobianXRhs.reset(new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["zero-matrix"], workMatrix["jacobianXBloc11"]));
  else
    jacobianXRhs.reset(new BlockMatrix(workMatrix["zero-matrix"], workMatrix["Id-matrix"], workMatrix["zero-matrix"], workMatrix["zero-matrix"]));
#endif

}

void LagrangianDS::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  // Memory allocation for p[0], p[1], p[2].
  initP(simulationType);

  // set q and q[1] to q0 and velocity0, initialize acceleration.
  *q[0] = *q0;
  *q[1] = *velocity0;

  // If z is NULL (ie has not been set), we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! z)
  {

#ifndef WithSmartPtr
    z = new SimpleVector(1);
    isAllocatedIn["z"] = true;
#else
    z.reset(new SimpleVector(1));
#endif

  }

  // Memory allocation for fL and its jacobians.
  initFL();

  // Set links to variables of top-class DynamicalSystem.
  // Warning: this results only in pointers links. No more memory allocation for vectors or matrices.
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).

  bool res = checkDynamicalSystem();
  if (!res) cout << "Warning: your dynamical system seems to be uncomplete (check = false)" << endl;

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  initRhs(time);

}

// --- GETTERS/SETTERS ---

void LagrangianDS::setQ(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ: inconsistent input vector size ");

  if (! q[0])
  {

#ifndef WithSmartPtr
    q[0] = new SimpleVector(newValue);
    isAllocatedIn["q"] = true;
#else
    q[0].reset(new SimpleVector(newValue));
#endif

  }
  else
    *q[0] = newValue;
}

void LagrangianDS::setQPtr(SiconosVectorSPtr newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQPtr: inconsistent input vector size ");


#ifndef WithSmartPtr
  if (isAllocatedIn["q"]) delete q[0];
  isAllocatedIn["q"] = false;
#endif

  q[0] = newPtr;

}

void LagrangianDS::setQ0(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0: inconsistent input vector size ");

  if (! q0)
  {

#ifndef WithSmartPtr
    q0 = new SimpleVector(newValue);
    isAllocatedIn["q0"] = true;
#else
    q0.reset(new SimpleVector(newValue));
#endif

  }
  else
    *q0 = newValue;
}

void LagrangianDS::setQ0Ptr(SiconosVectorSPtr newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0Ptr: inconsistent input vector size ");

#ifndef WithSmartPtr
  if (isAllocatedIn["q0"]) delete q0;
  isAllocatedIn["q0"] = false;
#endif

  q0 = newPtr;

}

void LagrangianDS::setQMemory(const SiconosMemory& newValue)
{
  if (! qMemory)
  {

#ifndef WithSmartPtr
    qMemory = new SiconosMemory(newValue);
    isAllocatedIn["qMemory"] = true;
#else
    qMemory.reset(new SiconosMemory(newValue));
#endif

  }
  else
    *qMemory = newValue;
}

void LagrangianDS::setQMemoryPtr(SiconosMemorySPtr newPtr)
{
#ifndef WithSmartPtr
  if (isAllocatedIn["qMemory"]) delete qMemory;
  isAllocatedIn["qMemory"] = false;
#endif

  qMemory = newPtr;

}

void LagrangianDS::setVelocity(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity: inconsistent input vector size ");

  if (! q[1])
  {

#ifndef WithSmartPtr
    q[1] = new SimpleVector(newValue);
    isAllocatedIn["velocity"] = true;
#else
    q[1].reset(new SimpleVector(newValue));
#endif

  }
  else
    *q[1] = newValue;
}

void LagrangianDS::setVelocityPtr(SiconosVectorSPtr newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityPtr: inconsistent input vector size ");

#ifndef WithSmartPtr
  if (isAllocatedIn["velocity"]) delete q[1];
  isAllocatedIn["velocity"] = false;
#endif

  q[1] = newPtr;


}

void LagrangianDS::setVelocity0(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0: inconsistent input vector size ");

  if (! velocity0)
  {

#ifndef WithSmartPtr
    velocity0 = new SimpleVector(newValue);
    isAllocatedIn["velocity0"] = true;
#else
    velocity0.reset(new SimpleVector(newValue));
#endif

  }
  else
    *velocity0 = newValue;
}

void LagrangianDS::setVelocity0Ptr(SiconosVectorSPtr newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0Ptr: inconsistent input vector size ");

#ifndef WithSmartPtr
  if (isAllocatedIn["velocity0"]) delete velocity0;
  isAllocatedIn["velocity0"] = false;
#endif

  velocity0 = newPtr;

}

SiconosVectorSPtr LagrangianDS::getAccelerationPtr() const
{
  return q[2];
}

void LagrangianDS::setVelocityMemory(const SiconosMemory& newValue)
{
  if (! velocityMemory)
  {

#ifndef WithSmartPtr
    velocityMemory = new SiconosMemory(newValue);
    isAllocatedIn["velocityMemory"] = true;
#else
    velocityMemory.reset(new SiconosMemory(newValue));
#endif

  }
  else
    *velocityMemory = newValue;
}

void LagrangianDS::setVelocityMemoryPtr(SiconosMemorySPtr newPtr)
{

#ifndef WithSmartPtr
  if (isAllocatedIn["velocityMemory"]) delete velocityMemory;
  isAllocatedIn["velocityMemory"] = false;
#endif
  velocityMemory = newPtr;

}

void LagrangianDS::setP(const SiconosVector& newValue, unsigned int level)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setP: inconsistent input vector size ");

  if (! p[level])
  {

#ifndef WithSmartPtr
    p[level] = new SimpleVector(newValue);
    string stringValue;
    stringstream sstr;
    sstr << level;
    sstr >> stringValue;
    stringValue = "p" + stringValue;
    isAllocatedIn[stringValue] = true;
#else
    p[level].reset(new SimpleVector(newValue));
#endif

  }
  else
    *(p[level]) = newValue;
}

void LagrangianDS::setPPtr(SiconosVectorSPtr newPtr, unsigned int level)
{

  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setPPtr: inconsistent input vector size ");

#ifndef WithSmartPtr
  string stringValue;
  stringstream sstr;
  sstr << level;
  sstr >> stringValue;
  stringValue = "p" + stringValue;

  if (isAllocatedIn[stringValue]) delete p[level];
  p[level] = newPtr;
  isAllocatedIn[stringValue] = false;
#else
  p[level] = newPtr;
#endif

}

void LagrangianDS::setMass(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setMass: inconsistent input matrix size ");

  if (! mass)
  {

#ifndef WithSmartPtr
    mass = new SimpleMatrix(newValue);
    isAllocatedIn["mass"] = true;
#else
    mass.reset(new SimpleMatrix(newValue));
#endif

  }
  else
    *mass = newValue;
  isPlugin["mass"] = false;
}

void LagrangianDS::setMassPtr(SiconosMatrixSPtr newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setMassPtr: inconsistent input matrix size ");

#ifndef WithSmartPtr
  if (isAllocatedIn["mass"]) delete mass;
  isAllocatedIn["mass"] = false;
#endif

  mass = newPtr;

  isPlugin["mass"] = false;
}

void LagrangianDS::setFInt(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFInt: inconsistent dimensions with problem size for input vector FInt");

  if (! fInt)
  {

#ifndef WithSmartPtr
    fInt = new SimpleVector(newValue);
    isAllocatedIn["fInt"] = true;
#else
    fInt.reset(new SimpleVector(newValue));
#endif

  }
  else
    *fInt = newValue;
  isPlugin["fInt"] = false;
}

void LagrangianDS::setFIntPtr(SiconosVectorSPtr newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFIntPtr: inconsistent input matrix size ");


#ifndef WithSmartPtr
  if (isAllocatedIn["fInt"]) delete fInt;
  isAllocatedIn["fInt"] = false;
#endif

  fInt = newPtr;

  isPlugin["fInt"] = false;
}

void LagrangianDS::setFExt(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFExt: inconsistent dimensions with problem size for input vector FExt");

  if (! fExt)
  {

#ifndef WithSmartPtr
    fExt = new SimpleVector(newValue);
    isAllocatedIn["fExt"] = true;
#else
    fExt.reset(new SimpleVector(newValue));
#endif

  }
  else
    *fExt = newValue;
  isPlugin["fExt"] = false;
}

void LagrangianDS::setFExtPtr(SiconosVectorSPtr newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setFIntPtr: inconsistent input matrix size ");

#ifndef WithSmartPtr
  if (isAllocatedIn["fExt"]) delete fExt;
  isAllocatedIn["fExt"] = false;
#endif

  fExt = newPtr;
  isPlugin["fExt"] = false;
}

void LagrangianDS::setNNL(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setNNL: inconsistent dimensions with problem size for input vector NNL");

  if (! NNL)
  {

#ifndef WithSmartPtr
    NNL = new SimpleVector(ndof);
    isAllocatedIn["NNL"] = true;
#else
    NNL.reset(new SimpleVector(ndof));
#endif

  }
  else
    *NNL = newValue;
  isPlugin["NNL"] = false;
}

void LagrangianDS::setNNLPtr(SiconosVectorSPtr newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setNNLPtr: inconsistent input matrix size ");

#ifndef WithSmartPtr
  if (isAllocatedIn["NNL"]) delete NNL;
  isAllocatedIn["NNL"] = false;
#endif

  NNL = newPtr;
  isPlugin["NNL"] = false;
}

void LagrangianDS::setJacobianFInt(unsigned int i, const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianFInt: inconsistent dimensions with problem size for input matrix JacobianQFInt");

  string name = "jacobianFInt" + toString<unsigned int>(i);
  if (! jacobianFInt[i])
  {

#ifndef WithSmartPtr
    jacobianFInt[i] = new SimpleMatrix(newValue);
    isAllocatedIn[name] = true;
#else
    jacobianFInt[i].reset(new SimpleMatrix(newValue));
#endif

  }
  else
    *jacobianFInt[i] = newValue;
  isPlugin[name] = false;
}

void LagrangianDS::setJacobianFIntPtr(unsigned int i, SiconosMatrixSPtr newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianFIntPtr: inconsistent input matrix size ");

  string name = "jacobianFInt" + toString<unsigned int>(i);

#ifndef WithSmartPtr
  if (isAllocatedIn[name]) delete jacobianFInt[i];
  isAllocatedIn[name] = false;
#endif

  jacobianFInt[i] = newPtr;
  isPlugin[name] = false;
}

void LagrangianDS::setJacobianNNL(unsigned int i, const SiconosMatrix& newValue)
{
  if (newValue.size(0) != ndof || newValue.size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianNNL: inconsistent dimensions with problem size for input matrix JacobianQNNL");

  string name = "jacobianNNL" + toString<unsigned int>(i);
  if (! jacobianNNL[i])
  {

#ifndef WithSmartPtr
    jacobianNNL[i] = new SimpleMatrix(newValue);
    isAllocatedIn[name] = true;
#else
    jacobianNNL[i].reset(new SimpleMatrix(newValue));
#endif
  }
  else
    *jacobianNNL[i] = newValue;
  isPlugin[name] = false;
}

void LagrangianDS::setJacobianNNLPtr(unsigned int i, SiconosMatrixSPtr newPtr)
{
  if (newPtr->size(0) != ndof || newPtr->size(1) != ndof)
    RuntimeException::selfThrow("LagrangianDS - setJacobianNNLPtr: inconsistent input matrix size ");

  string name = "jacobianNNL" + toString<unsigned int>(i);

#ifndef WithSmartPtr
  if (isAllocatedIn[name]) delete jacobianNNL[i];
  isAllocatedIn[name] = false;
#endif

  jacobianNNL[i] = newPtr;
  isPlugin[name] = false;
}

// --- Plugins related functions ---
void LagrangianDS::setComputeMassFunction(const string& pluginPath, const string& functionName)
{
  if (! mass)
  {

#ifndef WithSmartPtr
    mass = new SimpleMatrix(ndof, ndof);
    isAllocatedIn["mass"] = true;
#else
    mass.reset(new SimpleMatrix(ndof, ndof));
#endif

  }

  cShared.setFunction(&computeMassPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["mass"] = plugin + ":" + functionName;
  isPlugin["mass"] = true;
}

void LagrangianDS::setComputeFIntFunction(const string& pluginPath, const string& functionName)
{
  if (! fInt)
  {

#ifndef WithSmartPtr
    fInt = new SimpleVector(ndof);
    isAllocatedIn["fInt"] = true;
#else
    fInt.reset(new SimpleVector(ndof));
#endif

  }

  cShared.setFunction(&computeFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["fInt"] = plugin + ":" + functionName;
  isPlugin["fInt"] = true;
}

void LagrangianDS::setComputeFExtFunction(const string& pluginPath, const string& functionName)
{
  if (! fExt)
  {

#ifndef WithSmartPtr
    fExt = new SimpleVector(ndof);
    isAllocatedIn["fExt"] = true;
#else
    fExt.reset(new SimpleVector(ndof));
#endif

  }

  cShared.setFunction(&computeFExtPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["fExt"] = plugin + ":" + functionName;
  isPlugin["fExt"] = true;
}

void LagrangianDS::setComputeNNLFunction(const string& pluginPath, const string& functionName)
{
  if (! NNL)
  {

#ifndef WithSmartPtr
    NNL = new SimpleVector(ndof);
    isAllocatedIn["NNL"] = true;
#else
    NNL.reset(new SimpleVector(ndof));
#endif

  }

  cShared.setFunction(&computeNNLPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames["NNL"] = plugin + ":" + functionName;
  isPlugin["NNL"] = true;
}

void LagrangianDS::setComputeJacobianFIntFunction(unsigned int i, const string& pluginPath, const string& functionName)
{
  string name = "jacobianFInt" + toString<unsigned int>(i);

  if (! jacobianFInt[i])
  {

#ifndef WithSmartPtr
    jacobianFInt[i] = new SimpleMatrix(ndof, ndof);
    isAllocatedIn[name] = true;
#else
    jacobianFInt[i].reset(new SimpleMatrix(ndof, ndof));
#endif

  }

  cShared.setFunction(&computeJacobianFIntPtr[i], pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames[name] = plugin + ":" + functionName;
  isPlugin[name] = true;
}

void LagrangianDS::setComputeJacobianNNLFunction(unsigned int i, const string& pluginPath, const string& functionName)
{
  string name = "jacobianNNL" + toString<unsigned int>(i);
  if (! jacobianNNL[i])
  {

#ifndef WithSmartPtr
    jacobianNNL[i] = new SimpleMatrix(ndof, ndof);
    isAllocatedIn[name] = true;
#else
    jacobianNNL[i].reset(new SimpleMatrix(ndof, ndof));
#endif

  }

  cShared.setFunction(& (computeJacobianNNLPtr[i]), pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  pluginNames[name] = plugin + ":" + functionName;
  isPlugin[name] = true;
}

void LagrangianDS::setJacobianNNLFunction(unsigned int i, FPtr5 myF)
{
  computeJacobianNNLPtr[i] = myF;
}

void LagrangianDS::computeMass()
{
  if (isPlugin["mass"])
  {
    if (computeMassPtr == NULL)
      RuntimeException::selfThrow("computeMass() is not linked to a plugin function");
    computeMassPtr(ndof, &(*q[0])(0), &(*mass)(0, 0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeMass(SiconosVectorSPtr q2)
{
  if (isPlugin["mass"])
  {
    if (computeMassPtr == NULL)
      RuntimeException::selfThrow("computeMass() is not linked to a plugin function");

    computeMassPtr(ndof, &(*q2)(0), &(*mass)(0, 0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeFInt(double time)
{
  if (isPlugin["fInt"])
  {
    if (computeFIntPtr == NULL)
      RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

    computeFIntPtr(time, ndof, &(*q[0])(0), &(*q[1])(0), &(*fInt)(0), z->size(), &(*z)(0));
  }
}
void LagrangianDS::computeFInt(double time, SiconosVectorSPtr q2, SiconosVectorSPtr velocity2)
{
  if (isPlugin["fInt"])
  {
    if (computeFIntPtr == NULL)
      RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

    computeFIntPtr(time, ndof, &(*q2)(0), &(*velocity2)(0), &(*fInt)(0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeFExt(double time)
{
  if (isPlugin["fExt"])
  {
    if (computeFExtPtr == NULL)
      RuntimeException::selfThrow("computeFExt() is not linked to a plugin function");

    computeFExtPtr(time, ndof, &(*fExt)(0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeNNL()
{
  if (isPlugin["NNL"])
  {
    if (computeNNLPtr == NULL)
      RuntimeException::selfThrow("computeQ() is not linked to a plugin function");
    computeNNLPtr(ndof, &(*q[0])(0), &(*q[1])(0), &(*NNL)(0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeNNL(SiconosVectorSPtr q2, SiconosVectorSPtr velocity2)
{
  if (isPlugin["NNL"])
  {
    if (computeNNLPtr == NULL)
      RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

    computeNNLPtr(ndof, &(*q2)(0), &(*velocity2)(0), &(*NNL)(0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeJacobianFInt(unsigned int i, double time)
{
  string name = "jacobianFInt" + toString<unsigned int>(i);

  if (isPlugin[name])
  {
    if (computeJacobianFIntPtr[i] == NULL)
      RuntimeException::selfThrow("computeJacobianFInt(i,time) is not linked to a plugin function. i=" + i);

    (computeJacobianFIntPtr[i])(time, ndof, &(*q[0])(0), &(*q[1])(0), &(*jacobianFInt[i])(0, 0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeJacobianFInt(unsigned int i, double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  string name = "jacobianFInt" + toString<unsigned int>(i);
  if (isPlugin[name])
  {
    if (computeJacobianFIntPtr[i] == NULL)
      RuntimeException::selfThrow("computeJacobianFInt(i, ...) is not linked to a plugin function. i=" + i);

    computeJacobianFIntPtr[i](time, ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianFInt[i])(0, 0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeJacobianNNL(unsigned int i)
{
  string name = "jacobianNNL" + toString<unsigned int>(i);
  if (isPlugin[name])
  {
    if (computeJacobianNNLPtr[i] == NULL)
      RuntimeException::selfThrow("computeJacobianNNL(i) is not linked to a plugin function. i=" + i);

    computeJacobianNNLPtr[i](ndof, &(*q[0])(0), &(*q[1])(0), &(*jacobianNNL[i])(0, 0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeJacobianNNL(unsigned int i, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  string name = "jacobianNNL" + toString<unsigned int>(i);
  if (isPlugin[name])
  {
    if (computeJacobianNNLPtr[i] == NULL)
      RuntimeException::selfThrow("computeJacobianNNL(i, ...) is not linked to a plugin function. i=" + i);

    computeJacobianNNLPtr[i](ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianNNL[i])(0, 0), z->size(), &(*z)(0));
  }
}

void LagrangianDS::computeRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  *q[2] = *(p[2]); // Warning: r/p update is done in Interactions/Relations

  if (fL)
  {
    computeFL(time);
    *q[2] += *fL;
  }

  // mass and inv(mass) computatiton
  if (!isDSup) // if it is necessary to re-compute mass, FInt ..., ie if they have not been compiled during the present time step
    computeMass();


  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // -- Case 1: if mass is constant, then a copy of imass is LU-factorized during initialization and saved into workMatrix["invMass"].
  // -- Case 2: mass is not constant, we copy it into workMatrix["invMass"]
  // Then we proceed with PLUForwardBackward.
  if (isPlugin["mass"])
    *workMatrix["invMass"] = *mass;

  workMatrix["invMass"]->PLUForwardBackwardInPlace(*q[2]);

}

void LagrangianDS::computeJacobianXRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  if (!isDSup)
    computeMass();

  if (isPlugin["mass"])
    *workMatrix["invMass"] = *mass;

  if (jacobianFL[0])
  {
    SiconosMatrixSPtr bloc10 = jacobianXRhs->getBlockPtr(1, 0);
    computeJacobianFL(0, time);
    *bloc10 = *jacobianFL[0];
    workMatrix["invMass"]->PLUForwardBackwardInPlace(*bloc10);
  }

  if (jacobianFL[1])
  {
    SiconosMatrixSPtr bloc11 = jacobianXRhs->getBlockPtr(1, 1);
    computeJacobianFL(1, time);
    *bloc11 = *jacobianFL[1];
    workMatrix["invMass"]->PLUForwardBackwardInPlace(*bloc11);
  }
}

void LagrangianDS::computeFL(double time)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (fL)
  {
    // 1 - Computes the required functions
    if (isPlugin["fInt"])
      computeFInt(time);

    if (isPlugin["fExt"])
      computeFExt(time);

    if (isPlugin["NNL"])
      computeNNL();

    // 2 - set fL = fExt - fInt - NNL

    if (
#ifndef WithSmartPtr
      isAllocatedIn["fL"]
#else
      fL.use_count() == 1
#endif
    )
    {
      //if not that means that fL is already (pointer-)connected with
      // either fInt, NNL OR fExt.
      fL->zero();

      if (fInt)
        *fL -= *fInt;

      if (fExt)
        *fL += *fExt;

      if (NNL)
        *fL -= *NNL;
    }
  }
  // else nothing.
}

void LagrangianDS::computeFL(double time, SiconosVectorSPtr q2, SiconosVectorSPtr v2)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (fL)
  {
    // 1 - Computes the required functions
    if (isPlugin["fInt"])
      computeFInt(time, q2, v2);

    if (isPlugin["fExt"])
      computeFExt(time);

    if (isPlugin["NNL"])
      computeNNL(q2, v2);

    // 2 - set fL = fExt - fInt - NNL
    if (
#ifndef WithSmartPtr
      isAllocatedIn["fL"]
#else
      fL.use_count() == 1
#endif
    )
    {
      //if not that means that fL is already (pointer-)connected with
      // either fInt, NNL OR fExt.
      fL->zero();

      if (fInt)
        *fL -= *fInt;

      if (fExt)
        *fL += *fExt;

      if (NNL)
        *fL -= *NNL;
    }
  }
  // else nothing.
}

void LagrangianDS::computeJacobianFL(unsigned int i, double time)
{
  if (jacobianFL[i])
  {
    string name = "jacobianFInt" + toString<unsigned int>(i);
    if (isPlugin[name])
      computeJacobianFInt(i, time);
    name = "jacobianNNL" + toString<unsigned int>(i);
    if (isPlugin[name])
      computeJacobianNNL(i);
    if (
#ifndef WithSmartPtr
      isAllocatedIn[name]
#else
      jacobianFL[i].use_count() == 1
#endif
    )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      jacobianFL[i]->zero();
      if (jacobianFInt[i])
        *jacobianFL[i] -= *jacobianFInt[i];
      if (jacobianNNL[i])
        *jacobianFL[i] -= *jacobianNNL[i];
    }
  }
  //else nothing.
}

void LagrangianDS::saveSpecificDataToXML()
{
  // --- other data ---
  if (dsxml == NULL)
    RuntimeException::selfThrow("LagrangianDS::saveDSToXML - object DynamicalSystemXML does not exist");

  LagrangianDSXMLSPtr lgptr = boost::static_pointer_cast <LagrangianDSXML>(dsxml);
  lgptr->setMassPlugin(pluginNames["mass"]);
  lgptr->setQ(*q[0]);
  lgptr->setQ0(*q0);
  lgptr->setQMemory(*qMemory);
  lgptr->setVelocity(*q[1]);
  lgptr->setVelocity0(*velocity0);
  lgptr->setVelocityMemory(*velocityMemory);

  // FExt
  if (lgptr->hasFExt())
  {
    if (!lgptr->isFExtPlugin())
      lgptr->setFExtVector(*fExt);
  }
  else
    lgptr->setFExtPlugin(pluginNames["fExt"]);

  // FInt
  if (lgptr->hasFInt())
  {
    if (!lgptr->isFIntPlugin())
      if (fInt->size() > 0)
        lgptr->setFIntVector(*fInt);
      else cout << "Warning : FInt can't be saved, the FInt vector is not defined." << endl;
  }
  else
    lgptr->setFIntPlugin(pluginNames["fInt"]);

  for (unsigned int i = 0; i < 2; ++i)
  {
    // Jacobian FInt
    if (lgptr->hasJacobianFInt(i))
    {
      if (!lgptr->isJacobianFIntPlugin(i))
        lgptr->setJacobianFIntMatrix(i, *jacobianFInt[i]);
    }
    else
    {
      string name = "jacobianFInt" + toString<unsigned int>(i);
      lgptr->setJacobianFIntPlugin(i, pluginNames[name]);
    }

    // JacobianQNNL
    if (lgptr->hasJacobianNNL(i))
    {
      if (!lgptr->isJacobianNNLPlugin(i))
        lgptr->setJacobianNNLMatrix(i, *jacobianNNL[i]);
    }
    else
    {
      string name = "jacobianNNL" + toString<unsigned int>(i);
      lgptr->setJacobianNNLPlugin(i, pluginNames[name]);
    }
  }
  // NNL
  if (lgptr->hasNNL())
  {
    if (!lgptr->isNNLPlugin())
      lgptr->setNNLVector(*NNL);
  }
  else
    lgptr->setNNLPlugin(pluginNames["NNL"]);
}

void LagrangianDS::display() const
{
  cout << "=====> Lagrangian System display (number: " << number << ")." << endl;
  cout << "- ndof : " << ndof << endl;
  cout << "- q " << endl;
  if (q[0]) q[0]->display();
  else cout << "-> NULL" << endl;
  cout << "- q0 " << endl;
  if (q0) q0->display();
  cout << "- v " << endl;
  if (q[1]) q[1]->display();
  else cout << "-> NULL" << endl;
  cout << "- v0 " << endl;
  if (velocity0) velocity0->display();
  else cout << "-> NULL" << endl;
  cout << "- p " << endl;
  if (p[2]) p[2]->display();
  else cout << "-> NULL" << endl;
  cout << "===================================== " << endl;
}

// --- Functions for memory handling ---
void LagrangianDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" << endl;
  else
  {
#ifndef WithSmartPtr
    if (isAllocatedIn["qMemory"]) delete qMemory;
    qMemory = new SiconosMemory(steps);
    isAllocatedIn["qMemory"] = true;
    if (isAllocatedIn["velocityMemory"]) delete velocityMemory;
    velocityMemory = new SiconosMemory(steps);
    isAllocatedIn["velocityMemory"] = true;
#else
    qMemory.reset(new SiconosMemory(steps));
    velocityMemory.reset(new SiconosMemory(steps));
#endif

    swapInMemory();
  }
}

void LagrangianDS::swapInMemory()
{

  xMemory->swap(x[0]);
  qMemory->swap(q[0]);
  velocityMemory->swap(q[1]);
  // initialization of the reaction force due to the non smooth law
  p[1]->zero();

#ifndef WithSmartPtr
  if (isAllocatedIn["p2"]) p[2]->zero();
#endif

}

LagrangianDS* LagrangianDS::convert(DynamicalSystem* ds)
{
  LagrangianDS* lnlds = dynamic_cast<LagrangianDS*>(ds);
  return lnlds;
}

double LagrangianDS::dsConvergenceIndicator()
{
  double dsCvgIndic;
  SimpleVector diff(q[0]->size());
  // Compute difference between present and previous Newton steps
  SiconosVectorSPtr valRef = workVector["NewtonCvg"];
  sub(*(q[0]), *valRef, diff);
  if (valRef->norm2() != 0)
    dsCvgIndic = diff.norm2() / (valRef->norm2());
  else
    dsCvgIndic = diff.norm2();
  return (dsCvgIndic);
}

void LagrangianDS::computeQFree(double time, unsigned int level, SP::SiconosVector qFreeOut)
{
  // to compute qFree, derivative number level. Mainly used in EventDriven to compute second derivative
  // of q for Output y computation.

  if (qFreeOut->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS::computeQFree - Wrong size for output (different from ndof)");


  if (level != 2)
    RuntimeException::selfThrow("LagrangianDS::computeQFree - Only implemented for second derivative.");

  // Warning: we suppose that all other operators are up to date (FInt, FExt ...)

  qFreeOut->zero();
  if (fInt)
    *qFreeOut -= *fInt;
  if (fExt)
    *qFreeOut += *fExt;
  if (NNL)
    *qFreeOut -= *NNL;

  workMatrix["invMass"]->PLUForwardBackwardInPlace(*qFreeOut);
}

void LagrangianDS::resetNonSmoothPart()
{
  p[1]->zero();

#ifndef WithSmartPtr
  if (isAllocatedIn["p2"]) p[2]->zero();
#endif

}

void LagrangianDS::computePostImpactVelocity()
{
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  workMatrix["invMass"]->PLUForwardBackwardInPlace(*p[1]);
  *q[1] += *p[1];  // v+ = v- + p
}
