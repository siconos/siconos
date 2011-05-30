/* Siconos-Kernel, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY ory FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "GenericMechanical.hpp"
#include "FrictionContactXML.hpp"
#include "Topology.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonImpactFrictionNSL.hpp"

using namespace std;
using namespace RELATION;

//#define GMP_DEBUG

GenericMechanical::GenericMechanical(int FC3D_Solver_Id):
  LinearOSNS()
{
  _MStorageType = SICONOS_SPARSE;
  _pnumerics_GMP = buildEmptyGenericMechanicalProblem();
  genericMechnicalProblem_setDefaultSolverOptions(&*_numerics_solver_options, FC3D_Solver_Id);
}


void GenericMechanical::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  LinearOSNS::initialize(sim);




}

void GenericMechanical::computeDiagonalUnitaryBlock(const UnitaryRelationsGraph::VDescriptor& vd)
{
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());
  //bool isTimeInvariant = simulation()->model()->nonSmoothDynamicalSystem()->topology()->isTimeInvariant();

  /*Build the corresponding numerics problems*/

  SP::DynamicalSystem DS1 = indexSet->properties(vd).source;
  SP::DynamicalSystem DS2 = indexSet->properties(vd).target;
  SP::UnitaryRelation UR = indexSet->bundle(vd);

#ifdef GMP_DEBUG
  printf("GenericMechanical::computeUnitaryBlock: add problem of type ");
#endif
  if (!_hasBeUpdated)
  {
    int size = UR->getNonSmoothLawSize();
    if (Type::value(*(UR->interaction()->nonSmoothLaw()))
        == Type::EqualityConditionNSL)
    {
      addProblem(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_EQUALITY, size);
#ifdef GMP_DEBUG
      printf(" Type::EqualityConditionNSL\n");
#endif
      //pAux->size= UR->getNonSmoothLawSize();
    }
    else if (Type::value(*(UR->interaction()->nonSmoothLaw()))
             == Type::NewtonImpactNSL)
    {
#ifdef GMP_DEBUG
      printf(" Type::NewtonImpactNSL\n");
#endif
    }
    else if (Type::value(*(UR->interaction()->nonSmoothLaw()))
             == Type::NewtonImpactFrictionNSL)
    {
      FrictionContactProblem * pAux =
        (FrictionContactProblem *)addProblem(_pnumerics_GMP, SICONOS_NUMERICS_PROBLEM_FC3D, size);
      SP::NewtonImpactFrictionNSL nsLaw =
        boost::static_pointer_cast<NewtonImpactFrictionNSL> (UR->interaction()->nonSmoothLaw());
      pAux->dimension = 3;
      pAux->numberOfContacts = 1;
      *(pAux->mu) = nsLaw->mu();
#ifdef GMP_DEBUG
      printf(" Type::NewtonImpactFrictionNSL\n");
#endif
    }
    else
    {
      RuntimeException::selfThrow("GenericMechanical::computeDiagonalUnitaryBlock- not yet implemented for that NSLAW type");
    }
  }
  LinearOSNS::computeDiagonalUnitaryBlock(vd);
}

void GenericMechanical::computeUnitaryBlock(const UnitaryRelationsGraph::EDescriptor& ed)
{
  LinearOSNS::computeUnitaryBlock(ed);
}

int GenericMechanical::compute(double time)
{
  int info = 0;
  // --- Prepare data for GenericMechanical computing ---
  preCompute(time);
  _hasBeUpdated = true;
  /*
    La matrice _M est construite.  Ici, il faut construire les
    sous-problemes, c'est a dire completer les champs des
    NumericsProblem (_mu, _e, _en, les dimentions...).  Il faut aussi
    remplir la sous matrice M du sous-probleme.  Pour cela, on peut
    boucler sur les interactions et completer le membres
    _numerics_problem.problems[i] and
    _numerics_problem.problemsType[i].
   */

  //......


  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  if (_sizeOutput != 0)
  {
    // The GenericMechanical Problem in Numerics format

    _pnumerics_GMP->M = &*_M->getNumericsMatrix();
    _pnumerics_GMP->q = &*_q->getArray();

    // Call Numerics Driver for GenericMechanical
    //    display();
    info = genericMechanical_driver(_pnumerics_GMP,
                                    &*_z->getArray() ,
                                    &*_w->getArray() ,
                                    &*_numerics_solver_options);
    //printf("GenericMechanical::compute : R:\n");
    //_z->display();
    postCompute();

  }
  else
  {
#ifdef GMP_DEBUG
    printf("GenericMechanical::compute : sizeoutput is null\n");
#endif
  }

  return info;
}

void GenericMechanical::display() const
{
  cout << "===== " << "Generic mechanical Problem " << endl;
  LinearOSNS::display();
}


void  GenericMechanical::updateUnitaryBlocks()
{
  if (!_hasBeUpdated)
  {
    //    printf("GenericMechanical::updateUnitaryBlocks : must be updated\n");
    freeGenericMechanicalProblem(_pnumerics_GMP, NUMERICS_GMP_FREE_GMP);
    _pnumerics_GMP = buildEmptyGenericMechanicalProblem();
  }
  LinearOSNS::updateUnitaryBlocks();
}
void  GenericMechanical::computeAllUnitaryBlocks()
{
  assert(0);
  //printf("GenericMechanical::updateUnitaryBlocks : free and build a new GMP\n");
  freeGenericMechanicalProblem(_pnumerics_GMP, NUMERICS_GMP_FREE_GMP);
  _pnumerics_GMP = buildEmptyGenericMechanicalProblem();

  LinearOSNS::computeAllUnitaryBlocks();
}
GenericMechanical::~GenericMechanical()
{
  freeGenericMechanicalProblem(_pnumerics_GMP, NUMERICS_GMP_FREE_GMP);
  _pnumerics_GMP = 0;
  deleteSolverOptions(&*_numerics_solver_options);
}


