#include "LagrangianLinearEC.h"
using namespace std;

LagrangianLinearEC::LagrangianLinearEC()
{
  IN("LagrangianLinearEC::LagrangianLinearEC()\n");
  this->type = LAGRANGIANLINEAREC;
  OUT("LagrangianLinearEC::LagrangianLinearEC()\n");
}

LagrangianLinearEC::LagrangianLinearEC(EqualityConstraintXML* ecxml): LagrangianEC(ecxml)
{
  //  this->h = NULL;
  //  this->b = NULL;
  this->type = LAGRANGIANLINEAREC;
}

LagrangianLinearEC::~LagrangianLinearEC()
{
  IN("LagrangianLinearEC::~LagrangianLinearEC()\n");
  OUT("LagrangianLinearEC::~LagrangianLinearEC()\n");
}

LagrangianLinearEC::LagrangianLinearEC(SiconosMatrix h, SimpleVector b)
{
  this->h = h;
  this->b = b;
  this->type = LAGRANGIANLINEAREC;
}


SiconosMatrix* LagrangianLinearEC::getHPtr(void)
{
  return &this->h;
}

SiconosVector* LagrangianLinearEC::getBPtr(void)
{
  return &this->b;
}

void LagrangianLinearEC::computeOutput(double time)
{
  IN("LagrangianLinearEC::computeOutput\n");
  /*
    vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();

    DynamicalSystem *ds1 ,*ds2;
    SimpleVector *y = this->interaction->getYPtr();

    SimpleVector *yDot = this->interaction->getYDotPtr();

    if (vDS.size() == 2)
    {
        ds1=vDS[0];
        ds2=vDS[1];
        if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
        {
          LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
          LagrangianDS *d2 = static_cast<LagrangianDS*> (ds2);

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
        RuntimeException::selfThrow("LagrangianLinearEC::computeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
      }
    }
    else if (vDS.size() == 1)
    {
        ds1=vDS[0];
      if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
        {
          LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
          SimpleVector q(*(d1->getQPtr()));

        *y = (this->h * q) + this->b;

          SimpleVector vel(*(d1->getVelocityPtr()));

        *yDot = (this->h * vel);
        }
      else
      {
        RuntimeException::selfThrow("LagrangianLinearEC::computeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
      }
    }
    else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");
  */
  OUT("LagrangianLinearEC::computeOutput\n");
}



void LagrangianLinearEC::computeFreeOutput(double time)
{
  IN("LagrangianLinearEC::computeFreeOutput\n");
  /*
    vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();

    DynamicalSystem *ds1 ,*ds2;
    SiconosVector *y = this->interaction->getYPtr();
    SiconosVector *yDot = this->interaction->getYDotPtr();


    if (vDS.size() == 2)
    {
        ds1=vDS[0];
        ds2=vDS[1];
        if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
        {
          LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
          LagrangianDS *d2 = static_cast<LagrangianDS*> (ds2);
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
        RuntimeException::selfThrow("LagrangianLinearEC::computeFreeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());

      }
    }
    else if (vDS.size() == 1)
    {
        ds1=vDS[0];
      if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
        {
          LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
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
        RuntimeException::selfThrow("LagrangianLinearEC::computeFreeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
      }
    }
    else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");
  */
  OUT("LagrangianLinearEC::computeFreeOutput\n");
}

void LagrangianLinearEC::computeInput(double time)
{
  IN("LagrangianLinearEC::computeInput\n");
  /*
    vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();

    DynamicalSystem *ds1 ,*ds2;
    //SiconosVector r;
    SiconosVector *lambda = this->interaction->getLambdaPtr();



    if (vDS.size() == 2)
    {
        ds1=vDS[0];
        ds2=vDS[1];
        if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
        {
          LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
          LagrangianDS *d2 = static_cast<LagrangianDS*> (ds2);

        CompositeVector p;
        p.add(*(d1->getPPtr()));
        p.add(*(d2->getPPtr()));
          p += matTransVecMult(this->h, *lambda);

  //        cout<<" ### LagrangianLinearEC::computeInput => 'p'"<<endl;
  //        p.display();
          //p = matTransVecMult(this->h, *lambda);
        }
      else
      {
        // To be Finished
        RuntimeException::selfThrow("LagrangianLinearEC::computeInput not yet implemented for this type of dynamical system "+vDS[0]->getType());

      }
    }
    else if (vDS.size() == 1)
    {
        ds1=vDS[0];
      if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
        {
          LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
  //        SiconosVector p(*(d1->getPPtr()), false);
        SimpleVector p(*(d1->getPPtr()));
        //r = (this->h).multTranspose(*lambda);
          //p = lambda->matTransVecMult(this->h);
          p = matTransVecMult(this->h, *lambda);
        }
      else
      {
        // To be Finished
        RuntimeException::selfThrow("LagrangianLinearEC::computeFreeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
      }
    }
    else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");
  */
  OUT("LagrangianLinearEC::computeInput\n");
}

void LagrangianLinearEC::fillEqualityConstraintWithEqualityConstraintXML()
{
  OUT("LagrangianLinearEC::fillEqualityConstraintWithEqualityConstraintXML\n");
  EqualityConstraint::fillEqualityConstraintWithEqualityConstraintXML();
  if (this->ecXML != NULL)
  {
    //this->h = (static_cast<LagrangianLinearECXML*>(this->ecXML))->getH();
    //this->b = (static_cast<LagrangianLinearECXML*>(this->ecXML))->getB();
  }
  else RuntimeException::selfThrow("LagrangianLinearEC::fillEqualityConstraintWithEqualityConstraintXML - object EEqualityConstraintXML does not exist");
}

void LagrangianLinearEC::display() const
{
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the LagrangianLinearEC " << endl;
  cout << "| h " << endl;
  (this->h).display();
  cout << "| b " << endl;
  this->b.display();
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LagrangianLinearEC::saveEqualityConstraintToXML()
{
  IN("LagrangianLinearEC::saveEqualityConstraintToXML\n");
  EqualityConstraint::saveEqualityConstraintToXML();
  if (this->ecXML != NULL)
  {
    //(static_cast<LagrangianLinearECXML*>(this->ecXML))->setH( &(this->h) );
    //(static_cast<LagrangianLinearECXML*>(this->ecXML))->setB( &(this->b) );

    //    this->display();
  }
  else RuntimeException::selfThrow("LagrangianLinearEC::saveEqualityConstraintToXML - object EqualityConstraintXML does not exist");
  OUT("LagrangianLinearEC::saveEqualityConstraintToXML\n");
}

void LagrangianLinearEC::createEqualityConstraint(EqualityConstraintXML* ecXML, SiconosMatrix* H, SiconosVector* b)
{
  if (ecXML != NULL)
  {
    this->init();
    this->ecXML = ecXML;
    this->type = LAGRANGIANLINEAREC;
    this->fillEqualityConstraintWithEqualityConstraintXML();
  }
  else
  {
    this->h = *H;
    this->b = *b;
  }
}


LagrangianLinearEC* LagrangianLinearEC::convert(EqualityConstraint *r)
{
  cout << "LagrangianLinearEC::convert (EqualityConstraint *r)" << endl;
  LagrangianLinearEC* llr = dynamic_cast<LagrangianLinearEC*>(r);
  return llr;
}
