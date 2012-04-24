/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "LagrangianR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonEulerFrom1DLocalFrameR.hpp"
#include "OSNSMatrixProjectOnConstraints.hpp"
#include "LagrangianLinearTIDS.hpp"

using namespace std;
using namespace RELATION;
using namespace Siconos;
//#define MLCPPROJ_DEBUG

void MLCPProjectOnConstraints::initOSNSMatrix()
{
  _M.reset(new OSNSMatrixProjectOnConstraints(0, 0, _MStorageType));
  _n = 0;
  _m = 0;
  _curBlock = 0;
}


// Constructor from a set of data
MLCPProjectOnConstraints::MLCPProjectOnConstraints(const int newNumericsSolverId, double alphaval):
  MLCP(newNumericsSolverId), _alpha(alphaval)
{
  _levelMin = 1;
  _levelMax = 1;
}


void MLCPProjectOnConstraints::display() const
{
  cout << "======= MLCPProjectOnConstraints of size " << _sizeOutput << " with: " << endl;
  cout << "======= m " << _m << " _n " << _n << endl;
  LinearOSNS::display();
}
void MLCPProjectOnConstraints::updateUnitaryBlocks()
{
  // The present functions checks various conditions and possibly
  // compute unitaryBlocks matrices.
  //
  // Let URi and URj be two Unitary Relations.
  //
  // Things to be checked are:
  //  1 - is the topology time invariant?
  //  2 - does unitaryBlocks[URi][URj] already exists (ie has been
  //  computed in a previous time step)?
  //  3 - do we need to compute this unitaryBlock? A unitaryBlock is
  //  to be computed if URi and URj are in IndexSet1 AND if URi and
  //  URj have common DynamicalSystems.
  //
  // The possible cases are:
  //
  //  - If 1 and 2 are true then it does nothing. 3 is not checked.
  //  - If 1 == true, 2 == false, 3 == false, it does nothing.
  //  - If 1 == true, 2 == false, 3 == true, it computes the
  //    unitaryBlock.
  //  - If 1==false, 2 is not checked, and the unitaryBlock is
  //    computed if 3==true.
  //

  // Get index set from Simulation
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(_levelMin);


  // It seems that index() in not update in Index(0)
  // see comment in void Simulation::updateIndexSets()
  if (_levelMin == 0)
  {
    indexSet->update_vertices_indices();
    indexSet->update_edges_indices();
  }
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();

  // we put diagonal informations on vertices
  // self loops with bgl are a *nightmare* at the moment
  // (patch 65198 on standard boost install)

  if (indexSet->properties().symmetric)
  {
    RuntimeException::selfThrow
    (" MLCPProjectOnConstraints::updateUnitaryBlocks() - not yet implemented for symmetric case");
  }
  else // not symmetric => follow out_edges for each vertices
  {
    if (!_hasBeUpdated)
    {
      //      printf("MLCPProjectOnConstraints::updateUnitaryBlocks must be updated.\n");
      _n = 0;
      _m = 0;
      _curBlock = 0;
    }
    UnitaryRelationsGraph::VIterator vi, viend;
    for (boost::tie(vi, viend) = indexSet->vertices();
         vi != viend; ++vi)
    {

      SP::UnitaryRelation UR = indexSet->bundle(*vi);
      unsigned int nslawSize = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
                               (_M)->computeSizeForProjection(UR->interaction());

      if (! properties(*indexSet).blockProj[*vi])
      {
        properties(*indexSet).blockProj[*vi].reset(new SimpleMatrix(nslawSize, nslawSize));
      }

      if (!isLinear || !_hasBeUpdated)
      {
        computeDiagonalUnitaryBlock(*vi);
      }

      /* unitaryBlock must be zeroed at init */
      std::vector<bool> initialized;
      initialized.resize(indexSet->edges_number());
      std::fill(initialized.begin(), initialized.end(), false);

      /* on a undirected graph, out_edges gives all incident edges */
      UnitaryRelationsGraph::OEIterator oei, oeiend;
      for (boost::tie(oei, oeiend) = indexSet->out_edges(*vi);
           oei != oeiend; ++oei)
      {

        /* on adjoint graph there is at most 2 edges between source and target */
        UnitaryRelationsGraph::EDescriptor ed1, ed2;
        boost::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));

        assert(*oei == ed1 || *oei == ed2);

        /* the first edge as the lower index */
        assert(indexSet->index(ed1) <= indexSet->index(ed2));

        SP::UnitaryRelation UR1 = indexSet->bundle(indexSet->source(*oei));
        SP::UnitaryRelation UR2 = indexSet->bundle(indexSet->target(*oei));

        // Memory allocation if needed
        unsigned int nslawSize1 = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
                                  (_M)->computeSizeForProjection(UR1->interaction());
        unsigned int nslawSize2 = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
                                  (_M)->computeSizeForProjection(UR2->interaction());
        unsigned int isrc = indexSet->index(indexSet->source(*oei));
        unsigned int itar = indexSet->index(indexSet->target(*oei));

        SP::SiconosMatrix currentUnitaryBlock;

        if (itar > isrc) // upper block
        {
          if (! properties(*indexSet).upper_blockProj[ed1])
          {
            properties(*indexSet).upper_blockProj[ed1].reset(new SimpleMatrix(nslawSize1, nslawSize2));
            if (ed2 != ed1)
              properties(*indexSet).upper_blockProj[ed2] = properties(*indexSet).upper_blockProj[ed1];
          }
          currentUnitaryBlock = properties(*indexSet).upper_blockProj[ed1];

        }
        else  // lower block
        {
          if (! properties(*indexSet).lower_blockProj[ed1])
          {
            properties(*indexSet).lower_blockProj[ed1].reset(new SimpleMatrix(nslawSize1, nslawSize2));
            if (ed2 != ed1)
              properties(*indexSet).lower_blockProj[ed2] = properties(*indexSet).lower_blockProj[ed1];
          }
          currentUnitaryBlock = properties(*indexSet).lower_blockProj[ed1];
        }


        if (!initialized[indexSet->index(ed1)])
        {
          initialized[indexSet->index(ed1)] = true;
          currentUnitaryBlock->zero();
        }


        if (!isLinear || !_hasBeUpdated)
        {
          if (isrc != itar)
            computeUnitaryBlock(*oei);
        }

      }
    }
  }
#ifdef MLCPPROJ_DEBUG
  displayBlocks(indexSet);
#endif

}
void MLCPProjectOnConstraints::displayBlocks(SP::UnitaryRelationsGraph indexSet)
{

  std::cout <<  "MLCPProjectOnConstraints::displayBlocks(SP::UnitaryRelationsGraph indexSet) " << std::endl;
  std::cout << "                          indexSet :" << indexSet << std::endl;


  UnitaryRelationsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = indexSet->vertices();
       vi != viend; ++vi)
  {
    SP::UnitaryRelation UR = indexSet->bundle(*vi);
    std::cout << "                          vertex :" << *vi << std::endl;
    std::cout << "                          bundle :" << indexSet->bundle(*vi) << std::endl;

    if (properties(*indexSet).blockProj[*vi])
    {
      std::cout << "                          blockProj ";
      properties(*indexSet).blockProj[*vi]->display();
    }

    UnitaryRelationsGraph::OEIterator oei, oeiend;



    for (boost::tie(oei, oeiend) = indexSet->out_edges(*vi);
         oei != oeiend; ++oei)
    {
      unsigned int isrc = indexSet->index(indexSet->source(*oei));
      unsigned int itar = indexSet->index(indexSet->target(*oei));
      std::cout << "                          isrc :" << isrc << std::endl;
      std::cout << "                          itar :" << itar << std::endl;


      UnitaryRelationsGraph::EDescriptor ed1, ed2;
      std::cout << "                          outedges :" << *oei << std::endl;
      boost::tie(ed1, ed2) = indexSet->edges(indexSet->source(*oei), indexSet->target(*oei));
      std::cout << "                          edges(ed1,ed2) :" << ed1 << " " << ed2  << std::endl;
      std::cout << "                          (ed1).upper_blockProj : ";
      if (properties(*indexSet).upper_blockProj[ed1])
      {
        std::cout << properties(*indexSet).upper_blockProj[ed1] << "   :" ;
        properties(*indexSet).upper_blockProj[ed1]->display();
      }
      else
        std::cout << "NULL " << std::endl;

      std::cout << "                          (ed1).lower_blockProj : ";
      if (properties(*indexSet).lower_blockProj[ed1])
      {
        std::cout << properties(*indexSet).lower_blockProj[ed1] << "   :" ;
        properties(*indexSet).lower_blockProj[ed1]->display();
      }
      else
        std::cout << "NULL " << std::endl;

      std::cout << "                          (ed2).upper_blockProj : ";
      if (properties(*indexSet).upper_blockProj[ed2])
      {
        std::cout << properties(*indexSet).upper_blockProj[ed2] << "   :" ;
        properties(*indexSet).upper_blockProj[ed2]->display();
      }
      else
        std::cout << "NULL" << std::endl;

      std::cout << "                          (ed2).lower_blockProj : ";
      if (properties(*indexSet).lower_blockProj[ed2])
      {
        std::cout << properties(*indexSet).lower_blockProj[ed2] << "   :" ;
        properties(*indexSet).lower_blockProj[ed2]->display();
      }
      else
        std::cout << "NULL" << std::endl;
    }

  }
}
void MLCPProjectOnConstraints::updateUnitaryBlocksOLD()
{

  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());

  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();
  //  cout<<"isLinear: "<<isLinear<<" hasTopologyChanged: "<<hasTopologyChanged<<"hasBeUpdated: "<<_hasBeUpdated<<endl;


  if (indexSet->properties().symmetric)
  {
    RuntimeException::selfThrow
    ("MLCPProjectOnConstraints::updateUnitaryBlocks() - symmetric case for the indexSet is not yet implemented");
  }
  else // not symmetric => follow out_edges for each vertices
  {
    if (!_hasBeUpdated || !isLinear)
    {
      if (!_hasBeUpdated)
      {
        //      printf("MLCPProjectOnConstraints::updateUnitaryBlocks must be updated.\n");
        _n = 0;
        _m = 0;
        _curBlock = 0;
      }
      UnitaryRelationsGraph::VIterator vi, viend;
      for (boost::tie(vi, viend) = indexSet->vertices();
           vi != viend; ++vi)
      {
        SP::UnitaryRelation UR = indexSet->bundle(*vi);
        unsigned int sizeY = 0;
        sizeY = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
                (_M)->computeSizeForProjection(UR->interaction());

        // #ifdef MLCPPROJ_DEBUG
        //       cout<<"\nMLCPProjectOnConstraints::updateUnitaryBlocks()"<<endl;
        //       std::cout << "indexSet :"<< indexSet << std::endl;
        //       indexSet->display();
        //       std::cout << "vi :"<< *vi << std::endl;
        //       std::cout << "properties(*indexSet).blockProj[*vi]: before"<< properties(*indexSet).blockProj[*vi] << std::endl;
        // #endif

        if (! properties(*indexSet).blockProj[*vi])
        {
          properties(*indexSet).blockProj[*vi].reset(new SimpleMatrix(sizeY, sizeY));
        }
        // #ifdef MLCPPROJ_DEBUG
        //       std::cout << "properties(*indexSet).blockProj[*vi]: after"<< properties(*indexSet).blockProj[*vi] << std::endl;
        // #endif

        computeDiagonalUnitaryBlock(*vi);
      }





      UnitaryRelationsGraph::EIterator ei, eiend;
      for (boost::tie(ei, eiend) = indexSet->edges();
           ei != eiend; ++ei)
      {
        SP::UnitaryRelation UR1 = indexSet->bundle(indexSet->source(*ei));
        SP::UnitaryRelation UR2 = indexSet->bundle(indexSet->target(*ei));
        unsigned int sizeY1 = 0;
        unsigned int sizeY2 = 0;
        sizeY1 = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
                 (_M)->computeSizeForProjection(UR1->interaction());
        sizeY2 = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
                 (_M)->computeSizeForProjection(UR2->interaction());

        // Memory allocation if needed
        unsigned int isrc = indexSet->index(indexSet->source(*ei));
        unsigned int itar = indexSet->index(indexSet->target(*ei));
        if (itar > isrc) // upper block
        {
          if (! properties(*indexSet).upper_blockProj[*ei])
          {
            properties(*indexSet).upper_blockProj[*ei].reset(new SimpleMatrix(sizeY1, sizeY2));
          }
        }
        else  // lower block
        {
          if (! properties(*indexSet).lower_blockProj[*ei])
          {
            properties(*indexSet).lower_blockProj[*ei].reset(new SimpleMatrix(sizeY1, sizeY2));
          }
        }

        // Computation of the diagonal block
        computeUnitaryBlock(*ei);

        // allocation for transposed block
        // should be avoided
        if (itar > isrc) // upper block has been computed
        {
          // if (!properties(*indexSet).lower_blockProj[*ei])
          //   {
          //     properties(*indexSet).lower_blockProj[*ei].
          //  reset(new SimpleMatrix(properties(*indexSet).upper_blockProj[*ei]->size(1),
          //             properties(*indexSet).upper_blockProj[*ei]->size(0)));
          //   }
          properties(*indexSet).lower_blockProj[*ei].reset(new SimpleMatrix(*(properties(*indexSet).upper_blockProj[*ei])));
          properties(*indexSet).lower_blockProj[*ei]->trans();
          //          properties(*indexSet).lower_blockProj[*ei]->trans(*properties(*indexSet).upper_blockProj[*ei]);
        }
        else
        {
          assert(itar < isrc);    // lower block has been computed
          // if (!properties(*indexSet).upper_blockProj[*ei])
          //   {
          //     properties(*indexSet).upper_blockProj[*ei].
          //  reset(new SimpleMatrix(properties(*indexSet).lower_blockProj[*ei]->size(1),
          //             properties(*indexSet).lower_blockProj[*ei]->size(0)));
          //   }
          properties(*indexSet).upper_blockProj[*ei].
          reset(new SimpleMatrix(*(properties(*indexSet).lower_blockProj[*ei])));
          properties(*indexSet).upper_blockProj[*ei]->trans();
        }
        // #ifdef MLCPPROJ_DEBUG
        //             printf("MLCPP upper: %i %i\n",properties(*indexSet).upper_blockProj[*ei]->size(0),properties(*indexSet).upper_blockProj[*ei]->size(1));
        //             printf("MLCPP lower: %i %i\n",properties(*indexSet).lower_blockProj[*ei]->size(0),properties(*indexSet).lower_blockProj[*ei]->size(1));
        // #endif

      }

    }

  }

}

void MLCPProjectOnConstraints::computeDiagonalUnitaryBlock(const UnitaryRelationsGraph::VDescriptor& vd)
{
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());

  SP::DynamicalSystem DS1 = indexSet->properties(vd).source;
  SP::DynamicalSystem DS2 = indexSet->properties(vd).target;
  SP::UnitaryRelation UR = indexSet->bundle(vd);

  unsigned int sizeY = 0;
  sizeY = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
          (_M)->computeSizeForProjection(UR->interaction());


#ifdef MLCPPROJ_DEBUG
  cout << "\nMLCPProjectOnConstraints::computeDiagonalUnitaryBlock" << endl;
  std::cout << "levelMin()" << levelMin() << std::endl;
  std::cout << "indexSet :" << indexSet << std::endl;
  std::cout << "vd :" << vd << std::endl;
  indexSet->display();
  // std::cout << "DS1 :" << std::endl;
  // DS1->display();
  // std::cout << "DS2 :" << std::endl;
  // DS2->display();
#endif

  assert(properties(*indexSet).blockProj[vd]->size(0) == sizeY);
  assert(properties(*indexSet).blockProj[vd]->size(1) == sizeY);

  SP::SiconosMatrix currentUnitaryBlock = properties(*indexSet).blockProj[vd];
  if (!_hasBeUpdated)
    computeOptions(UR, UR);
  // Computes matrix _unitaryBlocks[UR1][UR2] (and allocates memory if
  // necessary) if UR1 and UR2 have commond DynamicalSystem.  How
  // _unitaryBlocks are computed depends explicitely on the type of
  // Relation of each UR.

  // Warning: we suppose that at this point, all non linear
  // operators (G for lagrangian relation for example) have been
  // computed through plug-in mechanism.

  // Get the W and Theta maps of one of the Unitary Relation -
  // Warning: in the current version, if OSI!=Moreau, this fails.
  // If OSI = MOREAU, centralUnitaryBlocks = W if OSI = LSODAR,
  // centralUnitaryBlocks = M (mass matrices)
  MapOfDSMatrices centralUnitaryBlocks;
  getOSIMaps(UR, centralUnitaryBlocks);

  SP::SiconosMatrix leftUnitaryBlock, rightUnitaryBlock, leftUnitaryBlock1;


  // General form of the unitaryBlock is : unitaryBlock =
  // a*extraUnitaryBlock + b * leftUnitaryBlock * centralUnitaryBlocks
  // * rightUnitaryBlock a and b are scalars, centralUnitaryBlocks a
  // matrix depending on the integrator (and on the DS), the
  // simulation type ...  left, right and extra depend on the relation
  // type and the non smooth law.




  currentUnitaryBlock->zero();

  // loop over the common DS
  bool endl = false;
  for (SP::DynamicalSystem ds = DS1; !endl; ds = DS2)
  {
    assert(ds == DS1 || ds == DS2);
    endl = (ds == DS2);

    if (Type::value(*ds) == Type::LagrangianLinearTIDS ||
        Type::value(*ds) == Type::LagrangianDS)
    {
      if (UR->getRelationType() != Lagrangian)
      {
        RuntimeException::selfThrow(
          "MLCPProjectOnConstraints::computeDiagonalUnitaryBlock - relation is not of type Lagrangian with a LagrangianDS.");
      }


      SP::LagrangianDS lds = (boost::static_pointer_cast<LagrangianDS>(ds));
      unsigned int sizeDS = lds->getDim();
      leftUnitaryBlock.reset(new SimpleMatrix(sizeY, sizeDS));
      UR->getLeftUnitaryBlockForDS(ds, leftUnitaryBlock);

      if (lds->boundaryConditions()) // V.A. Should we do that ?
      {
        for (vector<unsigned int>::iterator itindex =
               lds->boundaryConditions()->velocityIndices()->begin() ;
             itindex != lds->boundaryConditions()->velocityIndices()->end();
             ++itindex)
        {
          // (sizeY,sizeDS));
          SP::SiconosVector coltmp(new SimpleVector(sizeY));
          coltmp->zero();
          leftUnitaryBlock->setCol(*itindex, *coltmp);
        }
      }
      // (UR1 == UR2)
      SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
      //
      //        cout<<"LinearOSNS : leftUBlock\n";
      //        work->display();
      work->trans();
      //        cout<<"LinearOSNS::computeUnitaryBlock leftUnitaryBlock"<<endl;
      //        leftUnitaryBlock->display();
      centralUnitaryBlocks[ds]->PLUForwardBackwardInPlace(*work);
      //*currentUnitaryBlock +=  *leftUnitaryBlock ** work;


      prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
      //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftUnitaryBlock,*work,1.0,*currentUnitaryBlock);
      //*currentUnitaryBlock *=h;
    }
    else if (Type::value(*ds) == Type::NewtonEulerDS)
    {

      if (UR->getRelationType() != NewtonEuler)
      {
        RuntimeException::selfThrow("MLCPProjectOnConstraints::computeDiagonalUnitaryBlock - relation is not from NewtonEulerR.");
      }
      SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(ds));
#ifdef MLCPPROJ_WITH_CT

      unsigned int sizeDS = neds->getDim();
      SP::SimpleMatrix T = neds->T();
      SP::SimpleMatrix workT(new SimpleMatrix(*T));
      workT->trans();
      SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
      prod(*workT, *T, *workT2, true);
      leftUnitaryBlock.reset(new SimpleMatrix(sizeY, sizeDS));
      UR->getLeftUnitaryBlockForDS(ds, leftUnitaryBlock);
      SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
      cout << "LinearOSNS : leftUBlock\n";
      work->display();
      work->trans();
      cout << "LinearOSNS::computeUnitaryBlock workT2" << endl;
      workT2->display();
      workT2->PLUForwardBackwardInPlace(*work);
      prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
#else
      if (0) //(boost::static_pointer_cast<NewtonEulerR> UR->interaction()->relation())->_isConstact){
      {
        unsigned int sizeDS = neds->getDim();
        SP::SimpleMatrix T = neds->T();
        SP::SimpleMatrix workT(new SimpleMatrix(*T));
        workT->trans();
        SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
        prod(*workT, *T, *workT2, true);
        leftUnitaryBlock1.reset(new SimpleMatrix(sizeY, sizeDS));
        UR->getLeftUnitaryBlockForDS(ds, leftUnitaryBlock);
        leftUnitaryBlock.reset(new SimpleMatrix(1, sizeDS));
        for (unsigned int ii = 0; ii < sizeDS; ii++)
          leftUnitaryBlock->setValue(1, ii, leftUnitaryBlock1->getValue(1, ii));

        SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
        //cout<<"LinearOSNS : leftUBlock\n";
        //work->display();
        work->trans();
        //cout<<"LinearOSNS::computeUnitaryBlock workT2"<<endl;
        //workT2->display();
        workT2->PLUForwardBackwardInPlace(*work);
        prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
      }
      else
      {
        unsigned int sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
        leftUnitaryBlock.reset(new SimpleMatrix(sizeY, sizeDS));
        UR->getLeftUnitaryBlockForDSProjectOnConstraints(ds, leftUnitaryBlock);
        SP::SiconosMatrix work(new SimpleMatrix(*leftUnitaryBlock));
        //cout<<"LinearOSNS sizeY="<<sizeY<<": leftUBlock\n";
        //work->display();
        work->trans();
        prod(*leftUnitaryBlock, *work, *currentUnitaryBlock, false);
      }

    }
    else
      RuntimeException::selfThrow("MLCPProjectOnConstraints::computeDiagonalUnitaryBlock - ds is not from NewtonEulerDS neither a LagrangianDS.");



#endif
#ifdef MLCPPROJ_DEBUG
      std::cout << "MLCPProjectOnConstraints::computeDiagonalUnitaryBlock DiagunitaryBlock " << std::endl;
      currentUnitaryBlock->display();
#endif
    }
  }

  void MLCPProjectOnConstraints::computeUnitaryBlock(const UnitaryRelationsGraph::EDescriptor& ed)
  {

    // Computes matrix _unitaryBlocks[UR1][UR2] (and allocates memory if
    // necessary) if UR1 and UR2 have commond DynamicalSystem.  How
    // _unitaryBlocks are computed depends explicitely on the type of
    // Relation of each UR.

    // Warning: we suppose that at this point, all non linear
    // operators (G for lagrangian relation for example) have been
    // computed through plug-in mechanism.

    // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)
    SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(_levelMin);

    SP::DynamicalSystem ds = indexSet->bundle(ed);
    SP::UnitaryRelation UR1 = indexSet->bundle(indexSet->source(ed));
    SP::UnitaryRelation UR2 = indexSet->bundle(indexSet->target(ed));

    unsigned int index1 = indexSet->index(indexSet->source(ed));
    unsigned int index2 = indexSet->index(indexSet->target(ed));

    unsigned int sizeY1 = 0;
    sizeY1 = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
             (_M)->computeSizeForProjection(UR1->interaction());
    unsigned int sizeY2 = 0;
    sizeY2 = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
             (_M)->computeSizeForProjection(UR2->interaction());



    /*
      DynamicalSystemsSet commonDS;
      intersection(*UR1->dynamicalSystems(),*UR2->dynamicalSystems(), commonDS);
      assert (!commonDS.isEmpty()) ;
      for (DSIterator itDS = commonDS.begin(); itDS!=commonDS.end(); itDS++)
      {
      assert (*itDS == ds);
      }
    */
    MapOfDSMatrices centralUnitaryBlocks;
    getOSIMaps(UR1, centralUnitaryBlocks);

    SP::SiconosMatrix currentUnitaryBlock;

    assert(index1 != index2);

    if (index2 > index1) // upper block
    {
      //     if (! indexSet->properties(ed).upper_block)
      //     {
      //       indexSet->properties(ed).upper_block.reset(new SimpleMatrix(sizeY1, sizeY2));
      //     }

      assert(properties(*indexSet).upper_blockProj[ed]->size(0) == sizeY1);
      assert(properties(*indexSet).upper_blockProj[ed]->size(1) == sizeY2);

      currentUnitaryBlock = properties(*indexSet).upper_blockProj[ed];
    }
    else  // lower block
    {
      //     if (! indexSet->properties(ed).lower_block)
      //     {
      //       indexSet->properties(ed).lower_block.reset(new SimpleMatrix(sizeY1, sizeY2));
      //     }

      assert(properties(*indexSet).lower_blockProj[ed]->size(0) == sizeY1);
      assert(properties(*indexSet).lower_blockProj[ed]->size(1) == sizeY2);

      currentUnitaryBlock = properties(*indexSet).lower_blockProj[ed];
    }


    SP::SiconosMatrix leftUnitaryBlock, rightUnitaryBlock;

    RELATION::TYPES relationType1, relationType2;

    // General form of the unitaryBlock is : unitaryBlock =
    // a*extraUnitaryBlock + b * leftUnitaryBlock * centralUnitaryBlocks
    // * rightUnitaryBlock a and b are scalars, centralUnitaryBlocks a
    // matrix depending on the integrator (and on the DS), the
    // simulation type ...  left, right and extra depend on the relation
    // type and the non smooth law.
    relationType1 = UR1->getRelationType();
    relationType2 = UR2->getRelationType();
    if (relationType1 == NewtonEuler &&
        relationType2 == NewtonEuler)
    {
      assert(UR1 != UR2);
      currentUnitaryBlock->zero();
#ifdef MLCPPROJ_WITH_CT
      unsigned int sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getDim();
      leftUnitaryBlock.reset(new SimpleMatrix(sizeY1, sizeDS));
      UR1->getLeftUnitaryBlockForDS(ds, leftUnitaryBlock);
      SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(ds));
      SP::SimpleMatrix T = neds->T();
      SP::SimpleMatrix workT(new SimpleMatrix(*T));
      workT->trans();
      SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
      prod(*workT, *T, *workT2, true);
      rightUnitaryBlock.reset(new SimpleMatrix(sizeY2, sizeDS));
      UR2->getLeftUnitaryBlockForDS(ds, rightUnitaryBlock);
      rightUnitaryBlock->trans();
      workT2->PLUForwardBackwardInPlace(*rightUnitaryBlock);
      prod(*leftUnitaryBlock, *rightUnitaryBlock, *currentUnitaryBlock, false);

#else

      unsigned int sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
      leftUnitaryBlock.reset(new SimpleMatrix(sizeY1, sizeDS));
      UR1->getLeftUnitaryBlockForDSProjectOnConstraints(ds, leftUnitaryBlock);
      SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(ds));
      rightUnitaryBlock.reset(new SimpleMatrix(sizeY2, sizeDS));
      UR2->getLeftUnitaryBlockForDSProjectOnConstraints(ds, rightUnitaryBlock);
      rightUnitaryBlock->trans();
      prod(*leftUnitaryBlock, *rightUnitaryBlock, *currentUnitaryBlock, false);
    }
#endif
      else if (relationType1 == Lagrangian &&
               relationType2 == Lagrangian)
      {
        unsigned int sizeDS =  ds->getDim();
        leftUnitaryBlock.reset(new SimpleMatrix(sizeY1, sizeDS));
        UR1->getLeftUnitaryBlockForDS(ds, leftUnitaryBlock);

        Type::Siconos dsType = Type::value(*ds);
        if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
        {
          SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

          if (d->boundaryConditions()) // V.A. Should we do that ?
          {
            for (vector<unsigned int>::iterator itindex =
                   d->boundaryConditions()->velocityIndices()->begin() ;
                 itindex != d->boundaryConditions()->velocityIndices()->end();
                 ++itindex)
            {
              // (sizeY1,sizeDS));
              SP::SiconosVector coltmp(new SimpleVector(sizeY1));
              coltmp->zero();
              leftUnitaryBlock->setCol(*itindex, *coltmp);
            }
          }
        }

        // UR1 != UR2
        rightUnitaryBlock.reset(new SimpleMatrix(sizeY2, sizeDS));
        UR2->getLeftUnitaryBlockForDS(ds, rightUnitaryBlock);

        // Warning: we use getLeft for Right unitaryBlock
        // because right = transpose(left) and because of
        // size checking inside the getBlock function, a
        // getRight call will fail.
        rightUnitaryBlock->trans();
        centralUnitaryBlocks[ds]->PLUForwardBackwardInPlace(*rightUnitaryBlock);
        //*currentUnitaryBlock +=  *leftUnitaryBlock ** work;
        prod(*leftUnitaryBlock, *rightUnitaryBlock, *currentUnitaryBlock, false);
#ifdef MLCPPROJ_DEBUG
        std::cout << "MLCPProjectOnConstraints::computeUnitaryBlock : currentUnitaryBlock" << std::endl;
        currentUnitaryBlock->display();
#endif
      }

      else
        RuntimeException::selfThrow("MLCPProjectOnConstraints::computeUnitaryBlock not yet implemented for relation of type " + relationType1);

    }

    void MLCPProjectOnConstraints::computeqBlock(SP::UnitaryRelation UR, unsigned int pos)
    {

      assert(UR);
      unsigned int sizeY = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
                           (_M)->computeSizeForProjection(UR->interaction());
      for (unsigned int i = 0; i < sizeY; i++)
        _q->setValue(pos + i, UR->interaction()->y(0)->getValue(UR->getRelativePosition() + i));
#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::computeqBlock, _q from y(0)\n");
      _q->display();
#endif
    }

    void MLCPProjectOnConstraints::postCompute()
    {
      _hasBeUpdated = true;
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
      /** We chose to do a small step _alpha in view of stabilized the algorithm.*/
#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postCompute damping value = %f\n", _alpha);
#endif
      (*_z) *= _alpha;
      unsigned int pos = 0;
#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postCompute _z\n");
      _z->display();
      display();
#endif



      UnitaryRelationsGraph::VIterator ui, uiend;
#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postCompute BEFORE UPDATE:\n");
      for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
      {
        SP::UnitaryRelation ur = indexSet->bundle(*ui);
        RELATION::TYPES relationType = ur->getRelationType();
        if (relationType == NewtonEuler)
        {
          SP::Relation R = ur->interaction()->relation();
          SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(R));
          ner->computeh(0);
          ner->getq()->display();
        }
        else if (relationType == Lagrangian)
        {
          ur->interaction()->relation()->computeOutput(0.0, 0);
          for (DSIterator it = ur->interaction()->dynamicalSystemsBegin();
               it != ur->interaction()->dynamicalSystemsEnd();
               ++it)
          {
            SP::LagrangianDS lds =  boost::static_pointer_cast<LagrangianDS>(*it);
            lds->q()->display();
          }
        }
        else
        {
          RuntimeException::selfThrow("MLCPProjectOnConstraints::postCompute - relation type is not from Lagrangian type neither NewtonEuler.");
        }
      }
#endif


      for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
      {

        SP::UnitaryRelation ur = indexSet->bundle(*ui);
        // Get the relative position of UR-unitaryBlock in the vector w
        // or z
        pos = _M->getPositionOfUnitaryBlock(ur);
        RELATION::TYPES relationType = ur->getRelationType();
        if (relationType == NewtonEuler)
        {
          postComputeNewtonEulerR(ur, pos);
        }
        else if (relationType == Lagrangian)
        {
          postComputeLagrangianR(ur, pos);

        }
        else
        {
          RuntimeException::selfThrow("MLCPProjectOnConstraints::computeUnitaryBlock - relation type is not from Lagrangian type neither NewtonEuler.");
        }

      }
#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postCompute AFTER UPDATE:\n");
      for (boost::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
      {
        SP::UnitaryRelation ur = indexSet->bundle(*ui);
        RELATION::TYPES relationType = ur->getRelationType();
        if (relationType == NewtonEuler)
        {
          SP::Relation R = ur->interaction()->relation();
          SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(R));
          ner->computeh(0);
          ner->getq()->display();
        }
        else if (relationType == Lagrangian)
        {
          ur->interaction()->relation()->computeOutput(0.0, 0);
          for (DSIterator it = ur->interaction()->dynamicalSystemsBegin();
               it != ur->interaction()->dynamicalSystemsEnd();
               ++it)
          {
            SP::LagrangianDS lds =  boost::static_pointer_cast<LagrangianDS>(*it);
            lds->q()->display();
          }
        }
        else
        {
          RuntimeException::selfThrow("MLCPProjectOnConstraints::postCompute - relation type is not from Lagrangian type neither NewtonEuler.");
        }
      }
#endif


    }

    void MLCPProjectOnConstraints::postComputeLagrangianR(SP::UnitaryRelation ur, unsigned int pos)
    {

      SP::LagrangianR  lr = boost::static_pointer_cast<LagrangianR>(ur->interaction()->relation());
#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postComputeLagrangian lr->jachq \n");
      lr->jachq()->display();
#endif

#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postComputeLagrangianR q before update\n");
      for (DSIterator it = ur->interaction()->dynamicalSystemsBegin();
           it != ur->interaction()->dynamicalSystemsEnd();
           ++it)
      {
        SP::LagrangianDS lds =  boost::static_pointer_cast<LagrangianDS>(*it);
        lds->q()->display();
      }
#endif
      unsigned int sizeY = ur->interaction()->nonSmoothLaw()->size();
      SP::SimpleVector aBuff(new SimpleVector(sizeY));
      setBlock(*_z, aBuff, sizeY, pos, 0);
#ifdef MLCPPROJ_DEBUG
      printf("MLCPP lambda of ur is pos =%i :", pos);
      aBuff->display();
#endif



      // aBuff should nornally be in lambda[corectlevel]
      // The update of the position in DS should be made in Moreau::upateState or ProjectedMoreau::updateState
      SP::SiconosMatrix J = lr->jachq();
      SP::SimpleMatrix aux(new SimpleMatrix(*J));
      aux->trans();
      prod(*aux, *aBuff, *(lr->q()), false);
#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postComputeLagrangianR _z\n");
      _z->display();
      printf("MLCPProjectOnConstraints::postComputeLagrangianR updated\n");
      (lr->q())->display();
#endif



      //RuntimeException::selfThrow("MLCPProjectOnConstraints::postComputeLagrangianR() - not yet implemented");
    }
    void MLCPProjectOnConstraints::postComputeNewtonEulerR(SP::UnitaryRelation ur, unsigned int pos)
    {
      // Get Y and Lambda for the current Unitary Relation
      //y = ur->y(levelMin());
      //lambda = ur->lambda(levelMin());
      // Copy _w/_z values, starting from index pos into y/lambda.

      //      setBlock(*_w, y, y->size(), pos, 0);// Warning: yEquivalent is
      // saved in y !!
      //setBlock(*_z, lambda, lambda->size(), pos, 0);

      unsigned int sizeY = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
                           (_M)->computeSizeForProjection(ur->interaction());

      SP::NewtonEulerR ner = (boost::static_pointer_cast<NewtonEulerR>(ur->interaction()->relation()));
#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postCompute ner->jachq \n");
      ner->jachq()->display();
      ner->jachqT()->display();
#endif

#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postCompute q before update\n");
      (ner->getq())->display();
#endif


      SP::SimpleVector aBuff(new SimpleVector(sizeY));
      setBlock(*_z, aBuff, sizeY, pos, 0);
      /*Now, update the ds's dof throw the relation*/
      //proj_with_q SP::SiconosMatrix J=ner->jachqProj();
#ifdef MLCPPROJ_DEBUG
      printf("MLCPP lambda of ur is pos =%i :", pos);
      aBuff->display();
#endif

      DSIterator itDS;
      Index coord(8);
#ifdef MLCPPROJ_WITH_CT
      unsigned int k = 0;
      unsigned int NumDS = 0;
      SP::SiconosMatrix J = ner->jachqT();
      //printf("J\n");
      //J->display();

      SP::SimpleMatrix aux(new SimpleMatrix(*J));
      aux->trans();

      itDS = ur->interaction()->dynamicalSystemsBegin();
      while (itDS != ur->interaction()->dynamicalSystemsEnd())
      {
        SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(*itDS));
        k += neds->getqDim();
        itDS++;
      }
      SP::SimpleVector deltaqG(new SimpleVector(k));
      // prod(*aux,*aBuff,*vel,true);
      k = 0;
      deltaqG->zero();

      itDS = ur->interaction()->dynamicalSystemsBegin();
      while (itDS != ur->interaction()->dynamicalSystemsEnd())
      {
        Type::Siconos dsType = Type::value(**itDS);
        if (dsType != Type::NewtonEulerDS)
          RuntimeException::selfThrow("MLCPProjectOnConstraint::postCompute- ds is not from NewtonEulerDS.");
        SP::NewtonEulerDS neds = (boost::static_pointer_cast<NewtonEulerDS>(*itDS));

        SP::SimpleVector vel(new SimpleVector(neds->getDim()));
        SP::SimpleVector deltaQ(new SimpleVector(neds->getqDim()));
        coord[0] = k;
        coord[1] = k + neds->getDim();
        coord[2] = 0;
        coord[3] = aux->size(1);
        coord[4] = 0;
        coord[5] = aux->size(1);
        coord[6] = 0;
        coord[7] = neds->getDim();
        subprod(*aux   , *aBuff, *vel, coord, true);
        SP::SimpleMatrix T = neds->T();
        SP::SimpleMatrix workT(new SimpleMatrix(*T));
        workT->trans();
        SP::SimpleMatrix workT2(new SimpleMatrix(6, 6));
        prod(*workT, *T, *workT2, true);

        workT2->PLUForwardBackwardInPlace(*vel);
        prod(*T, *vel, *deltaQ);

        for (unsigned int ii = 0; ii < neds->getqDim(); ii++)
        {
          ner->getq()->setValue(k + ii,
                                deltaQ->getValue(ii) +
                                ner->getq()->getValue(k + ii));
          deltaqG->setValue(k + ii, deltaQ->getValue(ii) + deltaqG->getValue(k + ii));
        }

        k += neds->getDim();
        itDS++;
#ifdef MLCPPROJ_DEBUG
        printf("MLCPProjectOnConstraints::postCompute ds %d, q updated\n", NumDS);
        (ner->getq())->display();
#endif
        NumDS++;

      }
      printf("MLCPProjectOnConstraints deltaqG:\n");
      deltaqG->display();
#else

      SP::SiconosMatrix J = ner->jachq();
      SP::SimpleMatrix aux(new SimpleMatrix(*J));
      aux->trans();
      prod(*aux, *aBuff, *(ner->getq()), false);
#endif


#ifdef MLCPPROJ_DEBUG
      printf("MLCPProjectOnConstraints::postComputeNewtonR _z\n");
      _z->display();
      printf("MLCPProjectOnConstraints::postComputeNewtonR q updated\n");
      (ner->getq())->display();
#endif
    }

    void MLCPProjectOnConstraints::computeOptions(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2)
    {
      //  printf("MLCPProjectOnConstraints::computeOptions\n");
      // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)
      RELATION::TYPES relationType1;
      relationType1 = UR1->getRelationType();
      // Retrieve size of Y (projected variable)
      unsigned int sizeY1;
      sizeY1 = boost::static_pointer_cast<OSNSMatrixProjectOnConstraints>
               (_M)->computeSizeForProjection(UR1->interaction());

      // Compute the number of equalities
      unsigned int equalitySize1 =  sizeY1; //default behavior

      if (Type::value(*(UR1->interaction()->nonSmoothLaw())) == Type::NewtonImpactFrictionNSL ||
          Type::value(*(UR1->interaction()->nonSmoothLaw())) == Type::NewtonImpactNSL)
      {
        equalitySize1 = 0;
      }
      else if (Type::value(*(UR1->interaction()->nonSmoothLaw()))
               == Type::MixedComplementarityConditionNSL)
      {
        equalitySize1 =
          MixedComplementarityConditionNSL::convert(UR1->interaction()->nonSmoothLaw())->getEqualitySize();
      }

      // Compute the number of inequalities
      unsigned int inequalitySize1 =  sizeY1 - equalitySize1;



      if (UR1 == UR2)
      {
        //UR1->getExtraUnitaryBlock(currentUnitaryBlock);
        _m += inequalitySize1;
        _n += equalitySize1;
        //    _m=0;
        //_n=6;
        if (_curBlock > MLCP_NB_BLOCKS - 2)
          printf("MLCP.cpp : number of block to small, memory crach below!!!\n");
        /*add an equality block.*/

        // #ifdef MLCPPROJ_DEBUG
        //   printf("MLCPProjectOnConstraints::computeOptions()\n");
        // #endif

        if (equalitySize1 > 0)
        {
          _numerics_problem.blocksLine[_curBlock + 1] = _numerics_problem.blocksLine[_curBlock] + equalitySize1;
          _numerics_problem.blocksIsComp[_curBlock] = 0;
          // #ifdef MLCPPROJ_DEBUG
          //       std::cout << "_curBlock : " << _curBlock <<std::endl;
          //       std::cout << "_numerics_problem.blocksLine["<<_curBlock+1 <<" ] : " << _numerics_problem.blocksLine[_curBlock+1] <<std::endl;
          //       std::cout << "_numerics_problem.blocksIsComp["<<_curBlock <<" ] : " << _numerics_problem.blocksIsComp[_curBlock] <<std::endl;
          // #endif

          _curBlock++;
        }
        /*add a complementarity block.*/
        if (inequalitySize1 > 0)
        {
          _numerics_problem.blocksLine[_curBlock + 1] = _numerics_problem.blocksLine[_curBlock] + inequalitySize1;
          _numerics_problem.blocksIsComp[_curBlock] = 1;
          // #ifdef MLCPPROJ_DEBUG
          //       std::cout << "_curBlock : " << _curBlock <<std::endl;
          //       std::cout << "_numerics_problem.blocksLine["<<_curBlock+1<< "] : " << _numerics_problem.blocksLine[_curBlock+1] <<std::endl;
          //       std::cout << "_numerics_problem.blocksIsComp["<<_curBlock<< "] : " << _numerics_problem.blocksIsComp[_curBlock] <<std::endl;
          // #endif

          _curBlock++;

        }
      }
      // #ifdef MLCPPROJ_DEBUG
      //   std::cout << "_m : " << _m <<std::endl;
      //   std::cout << "_n : " << _n <<std::endl;
      // #endif
    }
