/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "MLCPProjectOnConstraints.hpp"
#include "MixedComplementarityConditionNSL.hpp"
#include "EqualityConditionNSL.hpp"
#include "Simulation.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
using namespace std;
using namespace RELATION;




// Constructor from a set of data
MLCPProjectOnConstraints::MLCPProjectOnConstraints(const int newNumericsSolverId):
  MLCP(newNumericsSolverId)
{

}


void MLCPProjectOnConstraints::display() const
{
  cout << "======= MLCPProjectOnConstraints of size " << _sizeOutput << " with: " << endl;
  cout << "======= m " << _m << " _n " << _n << endl;
  LinearOSNS::display();
}

MLCPProjectOnConstraints* MLCPProjectOnConstraints::convert(OneStepNSProblem* osnsp)
{
  MLCPProjectOnConstraints* mlcp = dynamic_cast<MLCPProjectOnConstraints*>(osnsp);
  return mlcp;
}

void MLCPProjectOnConstraints::computeUnitaryBlock(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2)
{
  computeOptions(UR1, UR2);
  // Computes matrix _unitaryBlocks[UR1][UR2] (and allocates memory if
  // necessary) if UR1 and UR2 have commond DynamicalSystem.  How
  // _unitaryBlocks are computed depends explicitely on the type of
  // Relation of each UR.

  // Get DS common between UR1 and UR2
  DynamicalSystemsSet commonDS;
  intersection(*UR1->dynamicalSystems(), *UR2->dynamicalSystems(), commonDS);
  assert(!commonDS.isEmpty()) ;

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)
  unsigned int nslawSize1 = UR1->getNonSmoothLawSizeProjectOnConstraints();
  unsigned int nslawSize2 = UR2->getNonSmoothLawSizeProjectOnConstraints();
  // Check allocation
  if (! _unitaryBlocks[UR1][UR2])
  {
    _unitaryBlocks[UR1][UR2].reset(new SimpleMatrix(nslawSize1, nslawSize2));
    _unitaryBlocks[UR1][UR2]->zero();
  }
  // Get the W and Theta maps of one of the Unitary Relation -
  // Warning: in the current version, if OSI!=Moreau, this fails.
  // If OSI = MOREAU, centralUnitaryBlocks = W if OSI = LSODAR,
  // centralUnitaryBlocks = M (mass matrices)
  MapOfDSMatrices centralUnitaryBlocks;
  getOSIMaps(UR1, centralUnitaryBlocks);

  SP::SiconosMatrix currentUnitaryBlock = _unitaryBlocks[UR1][UR2];

  SP::SiconosMatrix leftUnitaryBlock, rightUnitaryBlock;

  unsigned int sizeDS;
  RELATION::TYPES relationType1, relationType2;
  double h = simulation()->timeDiscretisation()->currentTimeStep();

  // General form of the unitaryBlock is : unitaryBlock =
  // a*extraUnitaryBlock + b * leftUnitaryBlock * centralUnitaryBlocks
  // * rightUnitaryBlock a and b are scalars, centralUnitaryBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.
  relationType1 = UR1->getRelationType();
  relationType2 = UR2->getRelationType();

  currentUnitaryBlock->zero();


  // loop over the common DS
  for (DSIterator itDS = commonDS.begin(); itDS != commonDS.end(); itDS++)
  {
    SP::DynamicalSystem ds = *itDS;
    Type::Siconos dsType = Type::value(*ds);
    if (dsType != Type::NewtonEulerDS)
      RuntimeException::selfThrow("MLCPProjectOnConstraints::computeUnitaryBlock - ds is not from NewtonEulerDS.");
    sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();

    // get _unitaryBlocks corresponding to the current DS
    // These _unitaryBlocks depends on the relation type.
    leftUnitaryBlock.reset(new SimpleMatrix(nslawSize1, sizeDS));
    UR1->getLeftUnitaryBlockForDSProjectOnConstraints(*itDS, leftUnitaryBlock);

    if (relationType1 == Lagrangian || relationType2 == Lagrangian || relationType1 == NewtonEuler || relationType2 == NewtonEuler)
    {
      if (UR1 == UR2)
      {
        SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
        //
        //        cout<<"LinearOSNS : leftUBlock\n";
        //        leftUnitaryBlock->display();
        work->trans();
        //        cout<<"LinearOSNS : leftUBlock'\n";
        //        work->display();
        //        cout<<"LinearOSNS::computeUnitaryBlock leftUnitaryBlock"<<endl;
        //        leftUnitaryBlock->display();

        //*currentUnitaryBlock +=  *leftUnitaryBlock ** work;
        prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
        //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftUnitaryBlock,*work,1.0,*currentUnitaryBlock);
        //*currentUnitaryBlock *=h;
        //        cout<<"LinearOSNS::computeUnitaryBlock unitaryBlock"<<endl;
        //        currentUnitaryBlock->display();

      }
      else
      {
        rightUnitaryBlock.reset(new SimpleMatrix(nslawSize2, sizeDS));
        UR2->getLeftUnitaryBlockForDSProjectOnConstraints(*itDS, rightUnitaryBlock);
        // Warning: we use getLeft for Right unitaryBlock
        // because right = transpose(left) and because of
        // size checking inside the getBlock function, a
        // getRight call will fail.
        rightUnitaryBlock->trans();

        //*currentUnitaryBlock +=  *leftUnitaryBlock ** work;
        prod(*leftUnitaryBlock, *rightUnitaryBlock, *currentUnitaryBlock, false);
      }
    }

    else RuntimeException::selfThrow("LinearOSNS::computeUnitaryBlock not yet implemented for relation of type " + relationType1);
  }
}
void MLCPProjectOnConstraints::computeqBlock(SP::UnitaryRelation UR, unsigned int pos)
{
  unsigned int sizeY = UR->getNonSmoothLawSizeProjectOnConstraints();
  for (int i = 0; i < sizeY; i++)
    _q->setValue(pos + i, UR->interaction()->y(0)->getValue(UR->getRelativePosition() + i));

}
void MLCPProjectOnConstraints::postCompute()
{

  // This function is used to set y/lambda values using output from
  // lcp_driver (w,z).  Only UnitaryRelations (ie Interactions) of
  // indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());

  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;

  // === Loop through "active" Unitary Relations (ie present in
  // indexSets[1]) ===

  unsigned int pos = 0;

  UnitaryRelationsGraph::VIterator ui, uiend;
  for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    SP::UnitaryRelation ur = indexSet->bundle(*ui);
    // Get the relative position of UR-unitaryBlock in the vector w
    // or z
    pos = _M->getPositionOfUnitaryBlock(ur);

    // Get Y and Lambda for the current Unitary Relation
    y = ur->y(levelMin());
    lambda = ur->lambda(levelMin());
    // Copy _w/_z values, starting from index pos into y/lambda.

    //      setBlock(*_w, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved in y !!
    setBlock(*_z, lambda, lambda->size(), pos, 0);

    SP::Relation R = ur->interaction()->relation();
    Type::Siconos RType = Type::value(*R);
    if (RType != Type::NewtonEulerR)
      RuntimeException::selfThrow("MLCPProjectOnConstraints::postCompute - R is not from NewtonEulerR.");
    SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(R));
    /*Now, update the ds's dof throw the relation*/
    SP::SiconosMatrix J = ner->jachq();
    SP::SimpleMatrix aux(new SimpleMatrix(*J));
    aux->trans();
    prod(*aux, *lambda, *(ner->getq()), false);
  }

}
