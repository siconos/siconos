#include "LagrangianLinearR.h"
//
#include "LagrangianDS.h"



using namespace std;

LagrangianLinearR::LagrangianLinearR():
  Relation(), H(NULL), b(NULL), isHAllocatedIn(false), isBAllocatedIn(false)
{
  relationType = LAGRANGIANLINEARRELATION;
}

LagrangianLinearR::LagrangianLinearR(RelationXML* relxml):
  Relation(relxml), H(NULL), b(NULL),
  isHAllocatedIn(true), isBAllocatedIn(true)
{
  relationType = LAGRANGIANLINEARRELATION;
  if (relxml != NULL)
  {
    // nothing about sizeH or sizeB in xml ... to review
    int row = ((static_cast<LagrangianLinearRXML*>(relationxml))->getH()).size(0);
    int col = ((static_cast<LagrangianLinearRXML*>(relationxml))->getH()).size(1);
    H = new SiconosMatrix(row, col);
    *H = (static_cast<LagrangianLinearRXML*>(relationxml))->getH();
    int size = ((static_cast<LagrangianLinearRXML*>(relationxml))->getB()).size();
    b = new SimpleVector(size);
    *b = (static_cast<LagrangianLinearRXML*>(relationxml))->getB();
  }
  else RuntimeException::selfThrow("LagrangianLinearR::xml constructor xml file=NULL");
}

LagrangianLinearR::LagrangianLinearR(SiconosMatrix* newH, SimpleVector* newB):
  Relation(), H(NULL), b(NULL), isHAllocatedIn(true), isBAllocatedIn(true)
{
  relationType = LAGRANGIANLINEARRELATION;
  H = new SiconosMatrix(newH->size(0), newH->size(1));
  *H = *newH;
  b = new SimpleVector(newB->size());
  *b = *newB;
}

LagrangianLinearR::~LagrangianLinearR()
{
  if (isHAllocatedIn) delete H;
  H = NULL;
  if (isBAllocatedIn) delete b;
  b = NULL;
}

SiconosMatrix LagrangianLinearR::getHRelatingToDS(const int& position)
{
  if (interaction->getDynamicalSystems()[ position ]->getType() != LNLDS
      && interaction->getDynamicalSystems()[ position ]->getType() != LTIDS)
    RuntimeException::selfThrow("LagrangianLinearR::getHRelatingToDS : Error! LagrangianLinear Relation linked to a Dynamical System which is not lagrangian");
  else
  {
    int row, col, gap;
    row = H->size(0);
    col = static_cast<LagrangianDS*>(interaction->getDynamicalSystems()[ position ])->getNdof();

    SiconosMatrix newH(row, col);

    /*
     * the gap is used to select the good part of the H matrix, according to the right DynamicalSystem
     */
    gap = col * position;
    for (int i = 0; i < row; i++)
      for (int j = 0; j < col; j++)
        newH(i, j) = (*H)(i, j + gap);
    return newH;
  }
}

void LagrangianLinearR::computeOutput(const double& time)
{
  IN("LagrangianLinearR::computeOutput\n");

  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();

  DynamicalSystem *ds1 , *ds2;
  SimpleVector *y = interaction->getYPtr();
  SimpleVector *yDot = interaction->getYDotPtr();

  if (vDS.size() == 2)
  {
    ds1 = vDS[0];
    ds2 = vDS[1];
    if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
    {
      LagrangianDS *d1 = static_cast<LagrangianDS*>(ds1);
      LagrangianDS *d2 = static_cast<LagrangianDS*>(ds2);

      CompositeVector q;
      q.add(*(d1->getQPtr()));
      q.add(*(d2->getQPtr()));
      //    cout<<"LagrangianLinearR::computeOutput ###### H then q ################"<<endl;
      //    this->h.display();
      //    q.display();
      //    cout<<"/LagrangianLinearR::computeOutput ##################################"<<endl;
      *y = (*H * q) + *b;

      CompositeVector vel;
      vel.add(*(d1->getVelocityPtr()));
      vel.add(*(d2->getVelocityPtr()));
      *yDot = (*H * vel);
    }
    else
    {
      //    SiconosVector x(*(vDS[0]->getXPtr()), false);
      //    x.add(*( vDS[1]->getXPtr()));
      // To be Finished
      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for this type of dynamical system " + vDS[0]->getType());
    }
  }
  else if (vDS.size() == 1)
  {
    ds1 = vDS[0];
    if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
    {
      LagrangianDS *d1 = static_cast<LagrangianDS*>(ds1);
      /*SiconosVector*/
      SimpleVector q(*(d1->getQPtr())/*, false*/);

      *y = (*H * q) + *b;

      /*SiconosVector*/
      SimpleVector vel(*(d1->getVelocityPtr())/*, false*/);

      *yDot = (*H * vel);
    }
    else
    {
      //SiconosVector x(*(vDS[0]->getXPtr()), false);
      // To be Finished
      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for this type of dynamical system " + vDS[0]->getType());
    }
  }
  else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");

  OUT("LagrangianLinearR::computeOutput\n");
}

void LagrangianLinearR::computeFreeOutput(const double& time)
{
  IN("LagrangianLinearR::computeFreeOutput\n");

  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();

  DynamicalSystem *ds1 , *ds2;
  SiconosVector * y = interaction->getYPtr();
  SiconosVector * yDot = interaction->getYDotPtr();

  if (vDS.size() == 2)
  {
    ds1 = vDS[0];
    ds2 = vDS[1];
    if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
    {
      LagrangianDS *d1 = static_cast<LagrangianDS*>(ds1);
      LagrangianDS *d2 = static_cast<LagrangianDS*>(ds2);
      //        SiconosVector qfree(*(d1->getQFreePtr()), false);
      CompositeVector qfree;
      qfree.add(*(d1->getQFreePtr()));
      qfree.add(*(d2->getQFreePtr()));
      *y = (*H * qfree) + *b;
      //        SiconosVector velfree(*(d1->getVelocityFreePtr()), false);
      CompositeVector velfree;
      velfree.add(*(d1->getVelocityFreePtr()));
      velfree.add(*(d2->getVelocityFreePtr()));
      *yDot = (*H * velfree);
    }
    else
    {
      //      SiconosVector xfree(*(vDS[0]->getXFreePtr()), false);
      //      xfree.add(*( vDS[1]->getXFreePtr()));
      // To be Finished
      RuntimeException::selfThrow("LagrangianLinearR::computeFreeOutput not yet implemented for this type of dynamical system " + vDS[0]->getType());

    }
  }
  else if (vDS.size() == 1)
  {
    ds1 = vDS[0];
    if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
    {
      LagrangianDS *d1 = static_cast<LagrangianDS*>(ds1);
      //        SiconosVector qfree(*(d1->getQFreePtr()), false);
      SimpleVector *qfree = d1->getVelocityFreePtr();
      *y = (*H * *qfree) + *b;
      //        SiconosVector velfree(*(d1->getVelocityFreePtr()), false);
      SimpleVector *velfree = d1->getVelocityFreePtr();
      *yDot = (*H * *velfree);
    }
    else
    {
      //      SiconosVector xfree(*(vDS[0]->getXFreePtr()), false);
      // To be Finished
      RuntimeException::selfThrow("LagrangianLinearR::computeFreeOutput not yet implemented for this type of dynamical system " + vDS[0]->getType());
    }
  }
  else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");

  OUT("LagrangianLinearR::computeFreeOutput\n");
}

void LagrangianLinearR::computeInput(const double& time)
{
  IN("LagrangianLinearR::computeInput\n");

  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();

  DynamicalSystem *ds1 , *ds2;
  SiconosVector *lambda = interaction->getLambdaPtr();

  if (vDS.size() == 2)
  {
    ds1 = vDS[0];
    ds2 = vDS[1];
    if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
    {
      LagrangianDS *d1 = static_cast<LagrangianDS*>(ds1);
      LagrangianDS *d2 = static_cast<LagrangianDS*>(ds2);

      CompositeVector p;
      p.add(*(d1->getPPtr()));
      p.add(*(d2->getPPtr()));
      p += matTransVecMult(*H, *lambda);
    }
    else
    {
      // To be Finished
      RuntimeException::selfThrow("LagrangianLinearR::computeInput not yet implemented for this type of dynamical system " + vDS[0]->getType());
    }
  }
  else if (vDS.size() == 1)
  {
    ds1 = vDS[0];
    if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
    {
      LagrangianDS *d1 = static_cast<LagrangianDS*>(ds1);
      SimpleVector p(*(d1->getPPtr()));
      p = matTransVecMult(*H, *lambda);
    }
    else
    {
      // To be Finished
      RuntimeException::selfThrow("LagrangianLinearR::computeFreeOutput not yet implemented for this type of dynamical system " + vDS[0]->getType());
    }
  }
  else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");



  OUT("LagrangianLinearR::computeInput\n");
}

void LagrangianLinearR::display() const
{
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the LagrangianLinearR " << endl;
  cout << "| h " << endl;
  if (H != NULL) H->display();
  else cout << "->NULL" << endl;
  cout << "| b " << endl;
  if (b != NULL) b->display();
  else cout << "->NULL" << endl;
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LagrangianLinearR::saveRelationToXML()
{
  IN("LagrangianLinearR::saveRelationToXML\n");
  if (relationxml != NULL)
  {
    (static_cast<LagrangianLinearRXML*>(relationxml))->setH(*H) ;
    (static_cast<LagrangianLinearRXML*>(relationxml))->setB(*b) ;
  }
  else RuntimeException::selfThrow("LagrangianLinearR::saveRelationToXML - object RelationXML does not exist");
  OUT("LagrangianLinearR::saveRelationToXML\n");
}

LagrangianLinearR* LagrangianLinearR::convert(Relation *r)
{
  cout << "LagrangianLinearR::convert (Relation *r)" << endl;
  LagrangianLinearR* llr = dynamic_cast<LagrangianLinearR*>(r);
  return llr;
}
