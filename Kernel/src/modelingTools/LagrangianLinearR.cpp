#include "LagrangianLinearR.h"

#include "DynamicalSystem.h"
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"
#include "LinearSystemDS.h"

#include "check.h"


LagrangianLinearR::LagrangianLinearR()
{
  IN("LagrangianLinearR::LagrangianLinearR()\n");
  //  this->h = NULL;
  //  this->b = NULL;
  this->relationType = LAGRANGIANLINEARRELATION; //"LagrangianLinearR";
  OUT("LagrangianLinearR::LagrangianLinearR()\n");
}

LagrangianLinearR::LagrangianLinearR(RelationXML* relxml): Relation(relxml)
{
  //  this->h = NULL;
  //  this->b = NULL;
  this->relationType = LAGRANGIANLINEARRELATION;
}

LagrangianLinearR::~LagrangianLinearR()
{
  IN("LagrangianLinearR::~LagrangianLinearR()\n");
  OUT("LagrangianLinearR::~LagrangianLinearR()\n");
}

LagrangianLinearR::LagrangianLinearR(SiconosMatrix h, SimpleVector b)
{
  this->h = h;
  this->b = b;
  this->relationType = LAGRANGIANLINEARRELATION;
}


SiconosMatrix* LagrangianLinearR::getHPtr(void)
{
  return &this->h;
}

SiconosMatrix LagrangianLinearR::getHRelatingToDS(int position)
{
  if (this->interaction->getDynamicalSystems()[ position ]->getType() != LNLDS
      && this->interaction->getDynamicalSystems()[ position ]->getType() != LTIDS)
    RuntimeException::selfThrow("LagrangianLinearR::getHRelatingToDS : Error! LagrangianLinear Relation linked to a Dynamical System which is not lagrangian");
  else
  {
    int row, col, gap;
    row = this->h.size(0);
    col = static_cast<LagrangianDS*>(this->interaction->getDynamicalSystems()[ position ])->getNdof();

    SiconosMatrix H(row, col);

    /*
     * the gap is used to select the good part of the H matrix, according to the right DynamicalSystem
     */
    gap = col * position;
    for (int i = 0; i < row; i++)
      for (int j = 0; j < col; j++)
        H(i, j) = this->h(i, j + gap);
    return H;
  }
}

SiconosVector* LagrangianLinearR::getBPtr(void)
{
  return &this->b;
}


void LagrangianLinearR::computeOutput(double time)
{
  IN("LagrangianLinearR::computeOutput\n");

  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();

  DynamicalSystem *ds1 , *ds2;
  /*SiconosVector*/
  SimpleVector *y = this->interaction->getYPtr();

  /*SiconosVector*/
  SimpleVector *yDot = this->interaction->getYDotPtr();

  if (vDS.size() == 2)
  {
    ds1 = vDS[0];
    ds2 = vDS[1];
    // \todo : pretty strange to use LNLDS with H and b !
    if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
    {
      LagrangianDS *d1 = static_cast<LagrangianDS*>(ds1);
      LagrangianDS *d2 = static_cast<LagrangianDS*>(ds2);

      CompositeVector q;
      q.add(*(d1->getQPtr()));
      q.add(*(d2->getQPtr()));
      *y = (this->h * q) + this->b;

      CompositeVector vel;
      vel.add(*(d1->getVelocityPtr()));
      vel.add(*(d2->getVelocityPtr()));
      *yDot = (this->h * vel);
    }
    else
    {
      //      SiconosVector x(*(vDS[0]->getXPtr()), false);
      //      x.add(*( vDS[1]->getXPtr()));
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

      *y = (this->h * q) + this->b;

      /*SiconosVector*/
      SimpleVector vel(*(d1->getVelocityPtr())/*, false*/);

      *yDot = (this->h * vel);
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



void LagrangianLinearR::computeFreeOutput(double time)
{
  IN("LagrangianLinearR::computeFreeOutput\n");



  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();

  DynamicalSystem *ds1 , *ds2;
  SiconosVector *y = this->interaction->getYPtr();
  SiconosVector *yDot = this->interaction->getYDotPtr();


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
      *y = (this->h * qfree) + this->b;
      //        SiconosVector velfree(*(d1->getVelocityFreePtr()), false);
      CompositeVector velfree;
      velfree.add(*(d1->getVelocityFreePtr()));
      velfree.add(*(d2->getVelocityFreePtr()));
      *yDot = (this->h * velfree);
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
      *y = (this->h * *qfree) + this->b;
      //        SiconosVector velfree(*(d1->getVelocityFreePtr()), false);
      SimpleVector *velfree = d1->getVelocityFreePtr();
      *yDot = (this->h * *velfree);
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

void LagrangianLinearR::computeInput(double time)
{
  IN("LagrangianLinearR::computeInput\n");

  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();

  DynamicalSystem *ds1 , *ds2;
  //SiconosVector r;
  SiconosVector *lambda = this->interaction->getLambdaPtr();



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
      p += matTransVecMult(this->h, *lambda);

      //        cout<<" ### LagrangianLinearR::computeInput => 'p'"<<endl;
      //        p.display();
      //p = matTransVecMult(this->h, *lambda);
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
      //        SiconosVector p(*(d1->getPPtr()), false);
      SimpleVector p(*(d1->getPPtr()));
      //r = (this->h).multTranspose(*lambda);
      //p = lambda->matTransVecMult(this->h);
      p = matTransVecMult(this->h, *lambda);
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

void LagrangianLinearR::fillRelationWithRelationXML()
{
  OUT("LagrangianLinearR::fillRelationWithRelationXML\n");
  Relation::fillRelationWithRelationXML();
  if (this->relationxml != NULL)
  {
    this->h = (static_cast<LagrangianLinearRXML*>(this->relationxml))->getH();
    this->b = (static_cast<LagrangianLinearRXML*>(this->relationxml))->getB();

    //    this->display();
  }
  else RuntimeException::selfThrow("LagrangianLinearR::fillRelationWithRelationXML - object RelationXML does not exist");
}

void LagrangianLinearR::display() const
{
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the LagrangianLinearR " << endl;
  cout << "| h " << endl;
  (this->h).display();
  cout << "| b " << endl;
  this->b.display();
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LagrangianLinearR::saveRelationToXML()
{
  IN("LagrangianLinearR::saveRelationToXML\n");
  Relation::saveRelationToXML();
  if (this->relationxml != NULL)
  {
    (static_cast<LagrangianLinearRXML*>(this->relationxml))->setH(&(this->h));
    (static_cast<LagrangianLinearRXML*>(this->relationxml))->setB(&(this->b));

    //    this->display();
  }
  else RuntimeException::selfThrow("LagrangianLinearR::saveRelationToXML - object RelationXML does not exist");
  OUT("LagrangianLinearR::saveRelationToXML\n");
}

void LagrangianLinearR::createRelation(LagrangianLinearRXML * relationXML,
                                       SiconosMatrix* H, SiconosVector* b)//, Interaction * interaction)
{
  if (relationXML != NULL)
  {
    this->init();
    this->relationxml = relationXML;
    this->relationType = LAGRANGIANLINEARRELATION;
    this->fillRelationWithRelationXML();
  }
  else
  {
    this->h = *H;
    this->b = *b;
  }
}


LagrangianLinearR* LagrangianLinearR::convert(Relation *r)
{
  cout << "LagrangianLinearR::convert (Relation *r)" << endl;
  LagrangianLinearR* llr = dynamic_cast<LagrangianLinearR*>(r);
  return llr;
}
