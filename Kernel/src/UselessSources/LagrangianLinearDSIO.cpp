/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
#include "LagrangianLinearDSIO.h"
using namespace std;

LagrangianLinearDSIO::LagrangianLinearDSIO()
{
  //  this->h = NULL;
  //  this->b = NULL;
  this->dsioType = LAGRANGIANLINEARDSIO; //"LagrangianLinearDSIO";
}

LagrangianLinearDSIO::LagrangianLinearDSIO(DSInputOutputXML* dsioXML): LagrangianDSIO()
{
  //  this->h = NULL;
  //  this->b = NULL;
  this->dsioType = LAGRANGIANLINEARDSIO;
}

LagrangianLinearDSIO::~LagrangianLinearDSIO()
{
  if (H != NULL) delete H;
  if (b != NULL) delete b;
}

SiconosMatrix* LagrangianLinearDSIO::getHPtr(void)
{
  return H;
}

SiconosVector* LagrangianLinearDSIO::getBPtr(void)
{
  return b;
}

void LagrangianLinearDSIO::computeOutput(double time)
{
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
  //        BlockVector q;
  //        q.add(*(d1->getQPtr()));
  //      q.add(*(d2->getQPtr()));
  //        *y = (this->h * q) + this->b;
  //
  //      BlockVector vel;
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
}



void LagrangianLinearDSIO::computeFreeOutput(double time)
{
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
  //      BlockVector qfree;
  //      qfree.add(*(d1->getQFreePtr()));
  //      qfree.add(*(d2->getQFreePtr()));
  //      *y = (this->h * qfree) + this->b;
  ////        SiconosVector velfree(*(d1->getVelocityFreePtr()), false);
  //      BlockVector velfree;
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
}

void LagrangianLinearDSIO::computeInput(double time)
{

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
  //      BlockVector p;
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
}

void LagrangianLinearDSIO::fillDSInputOutputWithDSInputOutputXML()
{
  DSInputOutput::fillDSInputOutputWithDSInputOutputXML();
  if (this->dsioxml != NULL)
  {
    H = new SimpleMatrix((static_cast<LagrangianLinearDSIOXML*>(this->dsioxml))->getH());
    b = new SimpleVector((static_cast<LagrangianLinearDSIOXML*>(this->dsioxml))->getB());

    //    this->display();
  }
  else RuntimeException::selfThrow("LagrangianLinearDSIO::fillDSInputOutputWithDSInputOutputXML - object SInputOutputXML does not exist");
}

void LagrangianLinearDSIO::display() const
{
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the LagrangianLinearDSIO " << endl;
  cout << "| h " << endl;
  this->H->display();
  cout << "| b " << endl;
  this->b->display();
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LagrangianLinearDSIO::saveDSInputOutputToXML()
{
  DSInputOutput::saveDSInputOutputToXML();
  if (this->dsioxml != NULL)
  {
    (static_cast<LagrangianLinearDSIOXML*>(this->dsioxml))->setH(H);
    (static_cast<LagrangianLinearDSIOXML*>(this->dsioxml))->setB(b);
  }
  else RuntimeException::selfThrow("LagrangianLinearDSIO::saveDSInputOutputToXML - object DSInputOutputXML does not exist");
}

void LagrangianLinearDSIO::createDSInputOutput(DSInputOutputXML * dsioXML,
    SiconosMatrix* newH, SiconosVector* newb)
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
    H = new SimpleMatrix(*newH);
    b = new SimpleVector(*newb);
  }
}


LagrangianLinearDSIO* LagrangianLinearDSIO::convert(DSInputOutput *dsio)
{
  LagrangianLinearDSIO* lldsio = dynamic_cast<LagrangianLinearDSIO*>(dsio);
  return lldsio;
}
