//$Id: LagrangianTIDS.cpp,v 1.38 2005/03/23 15:03:54 jbarbier Exp $
#include "LagrangianTIDS.h"

#include "check.h"


LagrangianTIDS::LagrangianTIDS()
{
  IN("LagrangianTIDS::LagrangianTIDS()\n");
  this->DSType = LTIDS;
  LagrangianNLDS::init();
  this->init();
  this->K = SiconosMatrix::SiconosMatrix();
  this->C = SiconosMatrix::SiconosMatrix();
  this->DSType = LTIDS;
  OUT("LagrangianTIDS::LagrangianTIDS()\n");
}

LagrangianTIDS::LagrangianTIDS(DSXML* dsxml): LagrangianNLDS(dsxml)
{
  IN("LagrangianTIDS::LagrangianTIDS(DSXML* dsxml):LagrangianNLDS(dsxml)\n");

  this->init();
  this->K = SiconosMatrix::SiconosMatrix();
  this->C = SiconosMatrix::SiconosMatrix();
  this->DSType = LTIDS;

  OUT("LagrangianTIDS::LagrangianTIDS(DSXML* dsxml):LagrangianNLDS(dsxml)\n");
}


LagrangianTIDS::~LagrangianTIDS()
{
  IN("LagrangianTIDS::~LagrangianTIDS()\n");

  OUT("LagrangianTIDS::~LagrangianTIDS()\n");
}

SiconosMatrix* LagrangianTIDS::getKPtr(void)
{
  return &this->K;
}

SiconosMatrix* LagrangianTIDS::getCPtr(void)
{
  return &this->C;
}


void LagrangianTIDS::init()
{
  IN("LagrangianTIDS::init\n");
  //LagrangianNLDS::init();
  this->computeFExtPtr = LTIDSComputeFExt;
  OUT("LagrangianTIDS::init\n");
}


void LagrangianTIDS::fillDSWithDSXML()
{
  IN("LagrangianTIDS::fillDSWithDSXML\n");

  if (this->dsxml != NULL)
  {
    this->K = (static_cast <LagrangianTIDSXML*>(this->dsxml))->getK();
    this->C = (static_cast <LagrangianTIDSXML*>(this->dsxml))->getC();
    //this->init();
    //    this->display();
  }
  else RuntimeException::selfThrow("LagrangianTIDS::fillDSWithDSXML - object DSXML does not exist");
  OUT("LagrangianTIDS::fillDSWithDSXML\n");
}

void LagrangianTIDS::display() const
{
  IN("LagrangianTIDS::display\n");

  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the LangrangianTIDS " << endl;
  LagrangianNLDS::display();
  cout << "| Stiffness Matrix K : " << endl;
  K.display();
  cout << "| Viscosity Matrix C : " << endl;
  C.display();
  cout << "-----------------------------------------------------" << endl << endl;
  OUT("LagrangianTIDS::display\n");
}

void LagrangianTIDS::saveDSToXML()
{
  IN("LagrangianTIDS::saveDSToXML\n");
  LagrangianNLDS::saveDSToXML();

  if (this->dsxml != NULL)
  {
    (static_cast <LagrangianTIDSXML*>(this->dsxml))->setK(&(this->K));
    (static_cast <LagrangianTIDSXML*>(this->dsxml))->setC(&(this->C));
  }
  else RuntimeException::selfThrow("LagrangianTIDS::saveDSToXML - object DSXML does not exist");
  OUT("LagrangianTIDS::saveDSToXML\n");
}


extern "C" double FextFunction(double time)
{
  double t = time;

  double res = -0.0;
  return res;
}

void LTIDSComputeFExt(int* sizeOfq, double* time, double* qPtr, double* fExt)
{
  IN("LagrangianTIDS : friend LTIDSComputeFInt\n");

  const double m = 1; // ball mass
  const double g = 9.8; // gravity

  cout << " /!\\!!! POUET !!!! [LagrangianTIDS : friend LTIDSComputeFInt]" << endl;
  getchar();

  int i;
  int n = *sizeOfq;
  double t = *time;

  for (i = 0; i < n; i++)
  {
    fExt[i] = 0.0;
  }

  fExt[0] = -m * g + FextFunction(t);
  OUT("LagrangianTIDS : friend LTIDSComputeFInt\n");
}

void LagrangianTIDS::createDynamicalSystem(DSXML * dsXML, int number, int ndof,
    SiconosVector* q0, SiconosVector* velocity0, SiconosMatrix* mass,
    string fExt, SiconosMatrix* K, SiconosMatrix* C) //,NSDS * nsds)
{
  if (dsXML != NULL)
  {
    LagrangianNLDS::createDynamicalSystem(dsXML);
    //this->init();
    LagrangianNLDS::fillDSWithDSXML();
    this->fillDSWithDSXML();
    this->linkDSXML();
    this->DSType = LTIDS;
  }
  else
  {
    this->number = number;
    this->ndof = ndof;

    LagrangianNLDS::SiconosVectorSizeInit();

    this->q0 = *q0;
    this->q = *q0;
    this->velocity0 = *velocity0;
    this->velocity = *velocity0;

    LagrangianNLDS::CompositeVectorInit();

    //this->setComputeMassFunction(this->cShared.getPluginName( mass ), this->cShared.getPluginFunctionName( mass ));
    this->mass = *mass;

    this->setComputeFExtFunction(this->cShared.getPluginName(fExt), this->cShared.getPluginFunctionName(fExt));

    //    this->setComputeJacobianQFIntFunction(this->cShared.getPluginName( jacobianQFInt ), this->cShared.getPluginFunctionName( jacobianQFInt ));
    //    this->setComputeJacobianVelocityFIntFunction(this->cShared.getPluginName( jacobianVelocityFInt ), this->cShared.getPluginFunctionName( jacobianQFInt ));
    //    this->setComputeJacobianQQNLInertiaFunction(this->cShared.getPluginName( jacobianQQNLInertia ), this->cShared.getPluginFunctionName( jacobianQQNLInertia ));
    //    this->setComputeJacobianVelocityQNLInertiaFunction(this->cShared.getPluginName( jacobianVelocityQNLInertia ), this->cShared.getPluginFunctionName( jacobianVelocityQNLInertia ));
    //
    //    this->setComputeQNLInertiaFunction(this->cShared.getPluginName( QNLInertia ), this->cShared.getPluginFunctionName( QNLInertia ));

    this->K = *K;
    this->C = *C;
    this->DSType = LTIDS;
  }
}


LagrangianTIDS* LagrangianTIDS::convert(DynamicalSystem* ds)
{
  cout << "LagrangianTIDS::convert (DynamicalSystem* ds)" << endl;
  LagrangianTIDS* ltids = dynamic_cast<LagrangianTIDS*>(ds);
  return ltids;
}

//$Log: LagrangianTIDS.cpp,v $
//Revision 1.38  2005/03/23 15:03:54  jbarbier
//- adaptation to the LMGC90 tags in non smooth dynamical system and strategy
//
//Revision 1.37  2005/02/11 17:36:00  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.36  2005/01/31 16:26:20  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.35  2004/09/28 08:21:28  jbarbier
//
//- manual creation of the BouncingBall example successful
//
//Revision 1.34  2004/09/27 13:27:13  jbarbier
//
//- Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//
//- new required tags of the model : title, author, description, date, xmlSchema.
//They replace previous attributes author, description and date of the Model.
//
//Revision 1.33  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.32  2004/09/22 10:54:43  jbarbier
//- light modification according to the attribute mass of the lagrangian dynamical
//systems. The lagrangianNLDS take always an function from a plugin to compute the
//mass, whereas the lagrangianTIDS needs only a matrix.
//
//- xml input files have been modified in consequence
//
//Revision 1.31  2004/09/10 11:26:13  charlety
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
//Revision 1.30  2004/08/23 14:30:01  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.29  2004/08/17 15:12:38  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.28  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.27  2004/07/30 14:37:14  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.26  2004/07/23 14:39:25  jbarbier
//- createModel, createNSDS, createDynamicalSystem, createBoundaryCondition OK
//it's the first step, it do the same thing that before, but the process is
//unified and it must simply add the traitment for the creation of the nodes in
//the DOM tree
//
//Revision 1.25  2004/07/12 11:22:53  charlety
//
//_ the state vector x of the dynamical system is now plugged to q and the velocity when this system is a Lagrangian one.
//
//_ A function to force a SiconosVector to be composite has been written. This is a temporary solution. We should use the operator = and the constructor by copy instead. This problem will be fixed later in the summer.
//
//Revision 1.24  2004/07/09 11:14:53  charlety
//
//_ Added a constructor by copy and an operator = in class SiconosMemory
//_ the getters on memory in DynamicalSystems return now some pointers
//
//Revision 1.23  2004/07/05 12:38:08  charlety
//
//try of false plugin developed in LagrangianTIDS. The Moreau integrator considers it like a LagrangianNLDS, but this is not the plugin which is used to compute the external strength, but a function of the LagrangianTIDS.
//
//Revision 1.22  2004/07/02 14:40:20  jbarbier
//some IN and OUT added
//BoundaryConditon saveToXML ok
//problem after, during the call of the destructors
//
//Revision 1.21  2004/06/29 15:05:54  acary
//Change  in the Display method of the Dynamical system
//
//Revision 1.20  2004/06/29 08:49:57  acary
//Ajout des commentaires Doxygen et des tages CVS
//