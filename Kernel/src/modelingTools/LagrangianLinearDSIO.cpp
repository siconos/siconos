#include "LagrangianLinearDSIO.h"
#include "DynamicalSystem.h"
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"
#include "LinearSystemDS.h"

#include "check.h"


LagrangianLinearDSIO::LagrangianLinearDSIO()
{
  IN("LagrangianLinearDSIO::LagrangianLinearDSIO()\n");
  //  this->h = NULL;
  //  this->b = NULL;
  this->dsioType = LAGRANGIANLINEARDSIO; //"LagrangianLinearDSIO";
  OUT("LagrangianLinearDSIO::LagrangianLinearDSIO()\n");
}

LagrangianLinearDSIO::LagrangianLinearDSIO(DSInputOutputXML* dsioXML): DSInputOutput(dsioXML)
{
  //  this->h = NULL;
  //  this->b = NULL;
  this->dsioType = LAGRANGIANLINEARDSIO;
}

LagrangianLinearDSIO::~LagrangianLinearDSIO()
{}

SiconosMatrix* LagrangianLinearDSIO::getHPtr(void)
{
  return &this->H;
}

SiconosVector* LagrangianLinearDSIO::getBPtr(void)
{
  return &this->b;
}

void LagrangianLinearDSIO::computeOutput(double time)
{
  IN("LagrangianLinearDSIO::computeOutput\n");

  //  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();
  //
  //  DynamicalSystem *ds1 ,*ds2;
  //  /*SiconosVector*/SimpleVector *y = this->interaction->getYPtr();
  //
  //  /*SiconosVector*/ SimpleVector *yDot = this->interaction->getYDotPtr();
  //
  //  if (vDS.size() == 2)
  //  {
  //      ds1=vDS[0];
  //      ds2=vDS[1];
  //      if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
  //      {
  //        LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
  //        LagrangianDS *d2 = static_cast<LagrangianDS*> (ds2);
  //
  //        CompositeVector q;
  //        q.add(*(d1->getQPtr()));
  //      q.add(*(d2->getQPtr()));
  //        *y = (this->h * q) + this->b;
  //
  //      CompositeVector vel;
  //      vel.add(*(d1->getVelocityPtr()));
  //      vel.add(*(d2->getVelocityPtr()));
  //      *yDot = (this->h * vel);
  //    }
  //    else
  //    {
  ////      SiconosVector x(*(vDS[0]->getXPtr()), false);
  ////      x.add(*( vDS[1]->getXPtr()));
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearDSIO::computeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //    }
  //  }
  //  else if (vDS.size() == 1)
  //  {
  //      ds1=vDS[0];
  //    if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
  //      {
  //        LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
  //        /*SiconosVector*/SimpleVector q(*(d1->getQPtr())/*, false*/);
  //
  //      *y = (this->h * q) + this->b;
  //
  //        /*SiconosVector*/SimpleVector vel(*(d1->getVelocityPtr())/*, false*/);
  //
  //      *yDot = (this->h * vel);
  //      }
  //    else
  //    {
  //      //SiconosVector x(*(vDS[0]->getXPtr()), false);
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearDSIO::computeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //    }
  //  }
  //  else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");

  OUT("LagrangianLinearDSIO::computeOutput\n");
}



void LagrangianLinearDSIO::computeFreeOutput(double time)
{
  IN("LagrangianLinearDSIO::computeFreeOutput\n");
  //
  //
  //
  //  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();
  //
  //  DynamicalSystem *ds1 ,*ds2;
  //  SiconosVector *y = this->interaction->getYPtr();
  //  SiconosVector *yDot = this->interaction->getYDotPtr();
  //
  //
  //  if (vDS.size() == 2)
  //  {
  //      ds1=vDS[0];
  //      ds2=vDS[1];
  //      if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
  //      {
  //        LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
  //        LagrangianDS *d2 = static_cast<LagrangianDS*> (ds2);
  ////        SiconosVector qfree(*(d1->getQFreePtr()), false);
  //      CompositeVector qfree;
  //      qfree.add(*(d1->getQFreePtr()));
  //      qfree.add(*(d2->getQFreePtr()));
  //      *y = (this->h * qfree) + this->b;
  ////        SiconosVector velfree(*(d1->getVelocityFreePtr()), false);
  //      CompositeVector velfree;
  //      velfree.add(*(d1->getVelocityFreePtr()));
  //      velfree.add(*(d2->getVelocityFreePtr()));
  //      *yDot = (this->h * velfree);
  //      }
  //    else
  //    {
  ////      SiconosVector xfree(*(vDS[0]->getXFreePtr()), false);
  ////      xfree.add(*( vDS[1]->getXFreePtr()));
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearDSIO::computeFreeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //
  //    }
  //  }
  //  else if (vDS.size() == 1)
  //  {
  //      ds1=vDS[0];
  //    if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
  //      {
  //        LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
  ////        SiconosVector qfree(*(d1->getQFreePtr()), false);
  //      SimpleVector *qfree = d1->getVelocityFreePtr();
  //      *y = (this->h * *qfree) + this->b;
  ////        SiconosVector velfree(*(d1->getVelocityFreePtr()), false);
  //      SimpleVector *velfree = d1->getVelocityFreePtr();
  //      *yDot = (this->h * *velfree);
  //      }
  //    else
  //    {
  ////      SiconosVector xfree(*(vDS[0]->getXFreePtr()), false);
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearDSIO::computeFreeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //    }
  //  }
  //  else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");
  OUT("LagrangianLinearDSIO::computeFreeOutput\n");
}

void LagrangianLinearDSIO::computeInput(double time)
{
  IN("LagrangianLinearDSIO::computeInput\n");

  //  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();
  //
  //  DynamicalSystem *ds1 ,*ds2;
  //  //SiconosVector r;
  //  SiconosVector *lambda = this->interaction->getLambdaPtr();
  //
  //
  //
  //  if (vDS.size() == 2)
  //  {
  //      ds1=vDS[0];
  //      ds2=vDS[1];
  //      if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
  //      {
  //        LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
  //        LagrangianDS *d2 = static_cast<LagrangianDS*> (ds2);
  //
  //      CompositeVector p;
  //      p.add(*(d1->getPPtr()));
  //      p.add(*(d2->getPPtr()));
  //        p += matTransVecMult(this->h, *lambda);
  //
  ////        cout<<" ### LagrangianLinearDSIO::computeInput => 'p'"<<endl;
  ////        p.display();
  //        //p = matTransVecMult(this->h, *lambda);
  //      }
  //    else
  //    {
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearDSIO::computeInput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //
  //    }
  //  }
  //  else if (vDS.size() == 1)
  //  {
  //      ds1=vDS[0];
  //    if ((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS))
  //      {
  //        LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
  ////        SiconosVector p(*(d1->getPPtr()), false);
  //      SimpleVector p(*(d1->getPPtr()));
  //      //r = (this->h).multTranspose(*lambda);
  //        //p = lambda->matTransVecMult(this->h);
  //        p = matTransVecMult(this->h, *lambda);
  //      }
  //    else
  //    {
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearDSIO::computeFreeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //    }
  //  }
  //  else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");
  OUT("LagrangianLinearDSIO::computeInput\n");
}

void LagrangianLinearDSIO::fillDSInputOutputWithDSInputOutputXML()
{
  OUT("LagrangianLinearDSIO::fillDSInputOutputWithDSInputOutputXML\n");
  DSInputOutput::fillDSInputOutputWithDSInputOutputXML();
  if (this->dsioxml != NULL)
  {
    this->H = (static_cast<LagrangianLinearDSIOXML*>(this->dsioxml))->getH();
    this->b = (static_cast<LagrangianLinearDSIOXML*>(this->dsioxml))->getB();

    //    this->display();
  }
  else RuntimeException::selfThrow("LagrangianLinearDSIO::fillDSInputOutputWithDSInputOutputXML - object SInputOutputXML does not exist");
}

void LagrangianLinearDSIO::display() const
{
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the LagrangianLinearDSIO " << endl;
  cout << "| h " << endl;
  this->H.display();
  cout << "| b " << endl;
  this->b.display();
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LagrangianLinearDSIO::saveDSInputOutputToXML()
{
  IN("LagrangianLinearDSIO::saveDSInputOutputToXML\n");
  DSInputOutput::saveDSInputOutputToXML();
  if (this->dsioxml != NULL)
  {
    (static_cast<LagrangianLinearDSIOXML*>(this->dsioxml))->setH(&(this->H));
    (static_cast<LagrangianLinearDSIOXML*>(this->dsioxml))->setB(&(this->b));
  }
  else RuntimeException::selfThrow("LagrangianLinearDSIO::saveDSInputOutputToXML - object DSInputOutputXML does not exist");
  OUT("LagrangianLinearDSIO::saveDSInputOutputToXML\n");
}

void LagrangianLinearDSIO::createDSInputOutput(DSInputOutputXML * dsioXML,
    SiconosMatrix* H, SiconosVector* b)
{
  if (dsioXML != NULL)
  {
    this->init();
    this->dsioxml = dsioXML;
    this->dsioType = LAGRANGIANLINEARDSIO;
    this->fillDSInputOutputWithDSInputOutputXML();
  }
  else
  {
    this->H = *H;
    this->b = *b;
  }
}


LagrangianLinearDSIO* LagrangianLinearDSIO::convert(DSInputOutput *dsio)
{
  cout << "LagrangianLinearDSIO::convert (DSInputOutput *dsio)" << endl;
  LagrangianLinearDSIO* lldsio = dynamic_cast<LagrangianLinearDSIO*>(dsio);
  return lldsio;
}
