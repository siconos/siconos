
#include "LagrangianLinearR.h"

#include "DynamicalSystem.h"
#include "LagrangianNLDS.h"
#include "LagrangianTIDS.h"
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
    col = static_cast<LagrangianNLDS*>(this->interaction->getDynamicalSystems()[ position ])->getNdof();

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
    if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
    {
      LagrangianNLDS *d1 = static_cast<LagrangianNLDS*>(ds1);
      LagrangianNLDS *d2 = static_cast<LagrangianNLDS*>(ds2);

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
      LagrangianNLDS *d1 = static_cast<LagrangianNLDS*>(ds1);
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
      LagrangianNLDS *d1 = static_cast<LagrangianNLDS*>(ds1);
      LagrangianNLDS *d2 = static_cast<LagrangianNLDS*>(ds2);
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
      LagrangianNLDS *d1 = static_cast<LagrangianNLDS*>(ds1);
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
      LagrangianNLDS *d1 = static_cast<LagrangianNLDS*>(ds1);
      LagrangianNLDS *d2 = static_cast<LagrangianNLDS*>(ds2);

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
      LagrangianNLDS *d1 = static_cast<LagrangianNLDS*>(ds1);
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

//$Log: LagrangianLinearR.cpp,v $
//Revision 1.23  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.22  2005/03/07 13:17:19  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.21  2005/03/01 10:38:19  jbarbier
//- RollingBalls sample is OK
//
//Revision 1.20  2005/02/28 16:22:33  jbarbier
//- rolling balls sample almost finished
//
//- in LCP, compute function now use all the interactions to make computations
//
//Revision 1.19  2005/02/15 15:15:32  charlety
//
//_ modified some very slow functions to increase performance
//
//Revision 1.18  2005/02/11 17:35:55  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.17  2005/02/10 10:35:19  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
//Revision 1.16  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.15  2005/02/04 09:44:43  jbarbier
//- rolling balls plugin renamed RBallPlugin
//
//- in Relation, some work to do to fill computOutput and Input functions !
//
//Revision 1.14  2005/02/04 07:46:20  jbarbier
//- last modification for RollingBalls
//
//Revision 1.13  2005/02/02 15:54:50  jbarbier
//- sample RollingBalls added
//
//- function getArray() added to SimpleVector to return the pointer on the array of double values
//
//Revision 1.12  2005/02/01 11:08:41  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.11  2005/01/31 16:26:19  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.10  2004/12/08 12:49:37  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.9  2004/09/16 11:35:24  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.8  2004/09/14 13:24:53  charlety
//
//_ changes in the interface of SiconosVector
//
//Revision 1.7  2004/09/10 11:26:12  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.6  2004/09/03 14:41:41  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.5  2004/08/17 15:12:37  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.4  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.3  2004/07/29 14:25:36  jbarbier
