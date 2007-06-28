/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
#include "LagrangianLinearEC.h"
using namespace std;

LagrangianLinearEC::LagrangianLinearEC()
{
  this->type = LAGRANGIANLINEAREC;
}

LagrangianLinearEC::LagrangianLinearEC(EqualityConstraintXML* ecxml): LagrangianEC(ecxml)
{
  //  this->h = NULL;
  //  this->b = NULL;
  this->type = LAGRANGIANLINEAREC;
}

LagrangianLinearEC::~LagrangianLinearEC()
{
  if (h != NULL) delete h;
  if (b != NULL) delete b;
}

LagrangianLinearEC::LagrangianLinearEC(const SiconosMatrix& newh, const SimpleVector& newb)
{
  h = new SimpleMatrix(newh);
  b = new SimpleVector(newb);
  type = LAGRANGIANLINEAREC;
}


SiconosMatrix* LagrangianLinearEC::getHPtr(void)
{
  return h;
}

SiconosVector* LagrangianLinearEC::getBPtr(void)
{
  return b;
}

void LagrangianLinearEC::computeOutput(double time)
{
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

          BlockVector q;
          q.add(*(d1->getQPtr()));
        q.add(*(d2->getQPtr()));
          *y = (this->h * q) + this->b;

        BlockVector vel;
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
}



void LagrangianLinearEC::computeFreeOutput(double time)
{
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
        BlockVector qfree;
        qfree.add(*(d1->getQFreePtr()));
        qfree.add(*(d2->getQFreePtr()));
        *y = (this->h * qfree) + this->b;
  //        SiconosVector velfree(*(d1->getVelocityFreePtr()), false);
        BlockVector velfree;
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
}

void LagrangianLinearEC::computeInput(double time)
{
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

        BlockVector p;
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
}

void LagrangianLinearEC::fillEqualityConstraintWithEqualityConstraintXML()
{
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
  h->display();
  cout << "| b " << endl;
  b->display();
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LagrangianLinearEC::saveEqualityConstraintToXML()
{
  EqualityConstraint::saveEqualityConstraintToXML();
  if (this->ecXML != NULL)
  {
    //(static_cast<LagrangianLinearECXML*>(this->ecXML))->setH( &(this->h) );
    //(static_cast<LagrangianLinearECXML*>(this->ecXML))->setB( &(this->b) );

    //    this->display();
  }
  else RuntimeException::selfThrow("LagrangianLinearEC::saveEqualityConstraintToXML - object EqualityConstraintXML does not exist");
}

void LagrangianLinearEC::createEqualityConstraint(EqualityConstraintXML* ecXML, SiconosMatrix* H, SiconosVector* newb)
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
    h = new SimpleMatrix(*H);
    b = new SimpleVector(*newb);
  }
}


LagrangianLinearEC* LagrangianLinearEC::convert(EqualityConstraint *r)
{
  cout << "LagrangianLinearEC::convert (EqualityConstraint *r)" << endl;
  LagrangianLinearEC* llr = dynamic_cast<LagrangianLinearEC*>(r);
  return llr;
}
