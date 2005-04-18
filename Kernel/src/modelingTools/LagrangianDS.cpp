#include "LagrangianDS.h"
//#include "LagrangianDSXML.h"

#include "check.h"


LagrangianDS::LagrangianDS()/*:DynamicalSystem()*/
{
  IN("LagrangianDS::LagrangianDS()\n");
  this->init();
  this->DSType = LNLDS;

  OUT("LagrangianDS::LagrangianDS()\n");
}

LagrangianDS::LagrangianDS(DSXML * dsXML)
{
  if (dsXML != NULL)
  {
    this->init();
    this->DSType = LNLDS;
    this->dsxml = dsXML;
    this->fillDSWithDSXML();
    this->linkDSXML();
  }
  else
  {
    cout << "LagrangianDS::LagrangianDS - DSXML paramater musn't be NULL" << endl;
  }
}

LagrangianDS::LagrangianDS(int number, int ndof,
                           SiconosVector* q0, SiconosVector* velocity0, string mass,
                           string fInt, string fExt,
                           string jacobianQFInt, string jacobianVelocityFInt,
                           string jacobianQQNLInertia, string jacobianVelocityQNLInertia,
                           string QNLInertia)
{
  this->dsxml = NULL;

  this->DSType = LNLDS;
  this->number = number;
  this->ndof = ndof;
  this->n = 2 * this->ndof;

  this->init();
  SiconosVectorSizeInit();

  this->q0 = *q0;
  this->q = *q0;
  this->velocity0 = *velocity0;
  this->velocity = *velocity0;


  CompositeVectorInit();

  this->setComputeMassFunction(this->cShared.getPluginName(mass), this->cShared.getPluginFunctionName(mass));

  this->setComputeFIntFunction(this->cShared.getPluginName(fInt), this->cShared.getPluginFunctionName(fInt));
  this->setComputeFExtFunction(this->cShared.getPluginName(fExt), this->cShared.getPluginFunctionName(fExt));

  this->setComputeJacobianQFIntFunction(this->cShared.getPluginName(jacobianQFInt), this->cShared.getPluginFunctionName(jacobianQFInt));
  this->setComputeJacobianVelocityFIntFunction(this->cShared.getPluginName(jacobianVelocityFInt), this->cShared.getPluginFunctionName(jacobianQFInt));
  this->setComputeJacobianQQNLInertiaFunction(this->cShared.getPluginName(jacobianQQNLInertia), this->cShared.getPluginFunctionName(jacobianQQNLInertia));
  this->setComputeJacobianVelocityQNLInertiaFunction(this->cShared.getPluginName(jacobianVelocityQNLInertia), this->cShared.getPluginFunctionName(jacobianVelocityQNLInertia));

  this->setComputeQNLInertiaFunction(this->cShared.getPluginName(QNLInertia), this->cShared.getPluginFunctionName(QNLInertia));
}


LagrangianDS::~LagrangianDS()
{
  IN("LagrangianDS::~LagrangianDS()\n");

  //  if (this->x != NULL) delete x;
  //  if (this->x0 != NULL) delete x0;
  //  if (this->xDot != NULL) delete xDot;
  //  if (this->xFree != NULL) delete xFree;
  //
  OUT("LagrangianDS::~LagrangianDS()\n");
}


void LagrangianDS::initMemory(int steps)
{
  IN("LagrangianDS::initMemory\n");
  DynamicalSystem::initMemory(steps);

  qMemory = SiconosMemory::SiconosMemory(steps, this->qMemory.getSiconosMemoryXML());
  velocityMemory = SiconosMemory::SiconosMemory(steps, this->velocityMemory.getSiconosMemoryXML());

  OUT("LagrangianDS::initMemory\n");
}


void LagrangianDS::swapInMemory(void)
{
  IN("LagrangianDS::swapInMemory(void)\n");

  // This operation should be made only if necessary. See todo note.
  DynamicalSystem::swapInMemory();

  this->qMemory.swap(&q);
  this->velocityMemory.swap(&velocity);

  // initialization of the reaction force due to the non smooth law at the beginning of each time step
  this->p.zero();

  OUT("LagrangianDS::swapInMemory(void)\n");
}



SiconosMatrix* LagrangianDS::getMassPtr(void)
{
  return &(this->mass);
}

SimpleVector* LagrangianDS::getQPtr(void)
{
  return &this->q;
}

SimpleVector* LagrangianDS::getQ0Ptr(void)
{
  return &this->q0;
}


SiconosMemory* LagrangianDS::getQMemories(void)
{
  return &this->qMemory;
}



SimpleVector* LagrangianDS::getVelocityPtr(void)
{
  return &this->velocity;
}

SimpleVector* LagrangianDS::getVelocity0Ptr(void)
{
  return &this->velocity0;
}


SiconosMemory* LagrangianDS::getVelocityMemories(void)
{
  return &this->velocityMemory;
}


SimpleVector* LagrangianDS::getFIntPtr(void)
{
  return &this->fInt;
}

SimpleVector* LagrangianDS::getFExtPtr(void)
{
  return &this->fExt;
}


SimpleVector* LagrangianDS::getQNLInertiaPtr(void)
{
  return &this->QNLInertia;
}


SiconosMatrix* LagrangianDS::getJacobianQFIntPtr(void)
{
  return &this->jacobianQFInt;
}

SiconosMatrix* LagrangianDS::getJacobianVelocityFIntPtr(void)
{
  return &this->jacobianVelocityFInt;
}

SiconosMatrix* LagrangianDS::getJacobianQQNLInertiaPtr(void)
{
  return &this->jacobianQQNLInertia;
}

SiconosMatrix* LagrangianDS::getJacobianVelocityQNLInertiaPtr(void)
{
  return &this->jacobianVelocityQNLInertia;
}

//////////////////////////

void LagrangianDS::computeMass(double time)
{
  if (computeMassPtr == NULL)
    RuntimeException::selfThrow("computeMass() is not linked to a plugin function");

  int size = q.size();
  this->computeMassPtr(&size, &time, &q(0), &mass(0, 0));
}

void LagrangianDS::computeFInt(double time)
{
  if (computeFIntPtr == NULL)
    RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

  int size = this->q.size();
  this->computeFIntPtr(&size, &time, &this->q(0), &this->velocity(0), &this->fInt(0));
}

void LagrangianDS::computeFExt(double time)
{
  IN("LagrangianDS::computeFExt(double time)\n");
  if (computeFExtPtr == NULL)
    RuntimeException::selfThrow("computeFExt() is not linked to a plugin function");

  int size = q.size();

  this->computeFExtPtr(&size, &time, &q(0), &fExt(0));

  OUT("LagrangianDS::computeFExt(double time)\n");

}

void LagrangianDS::computeQNLInertia()
{
  if (computeQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

  int size = q.size();
  this->computeQNLInertiaPtr(&size, &q(0), &velocity(0), &QNLInertia(0));
}

void LagrangianDS::computeJacobianQFInt(double time)
{
  if (computeJacobianQFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

  int size = q.size();
  this->computeJacobianQFIntPtr(&size, &time, &q(0), &velocity(0), &jacobianQFInt(0, 0));
}

void LagrangianDS::computeJacobianVelocityFInt(double time)
{
  if (computeJacobianVelocityFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

  // to do
  //this->computeJacobianVelocityFIntPtr();
}

void LagrangianDS::computeJacobianQQNLInertia(double time)
{
  if (computeJacobianQQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQQNLInertia() is not linked to a plugin function");

  //to do
  //this->computeJacobianQQPtr(time);
}

void LagrangianDS::computeJacobianVelocityQNLInertia(double time)
{
  if (computeJacobianVelocityQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityQNLInertia() is not linked to a plugin function");

  //to do
  //this->computeJacobianVelocityQNLInertiaPtr(time);
}


//////


void LagrangianDS::setComputeMassFunction(string pluginPath, string functionName)
{
  IN("LagrangianDS::setComputeMassFunction\n");
  this->computeMassPtr = NULL;
  cShared.setFunction(&computeMassPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->massFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeMassFunction\n");

}

void LagrangianDS::setComputeFIntFunction(string pluginPath, string functionName)
{
  IN("LagrangianDS::setComputeFIntFunction\n");
  this->computeFIntPtr = NULL;
  cShared.setFunction(&computeFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->fIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeFIntFunction\n");
}

void LagrangianDS::setComputeFExtFunction(string pluginPath, string functionName)
{
  IN("LagrangianDS::setComputeFExtFunction\n");
  this->computeFExtPtr = NULL;
  cShared.setFunction(&computeFExtPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->fExtFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeFExtFunction\n");
}

void LagrangianDS::setComputeQNLInertiaFunction(string pluginPath, string functionName)
{
  IN("LagrangianDS::setComputeQNLInertiaFunction\n");
  this->computeQNLInertiaPtr = NULL;
  cShared.setFunction(&computeQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->QNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeQNLInertiaFunction\n");
}

void LagrangianDS::setComputeJacobianQFIntFunction(string pluginPath, string functionName)
{
  IN("LagrangianDS::setComputeJacobianQFIntFunction\n");
  this->computeJacobianQFIntPtr = NULL;
  cShared.setFunction(&computeJacobianQFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->jacobianQFIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianQFIntFunction\n");
}

void LagrangianDS::setComputeJacobianVelocityFIntFunction(string pluginPath, string functionName)
{
  IN("LagrangianDS::setComputeMassFunction\n");
  this->computeJacobianVelocityFIntPtr = NULL;
  cShared.setFunction(&computeJacobianVelocityFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->jacobianVelocityFIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeMassFunction\n");
}

void LagrangianDS::setComputeJacobianQQNLInertiaFunction(string pluginPath, string functionName)
{
  IN("LagrangianDS::setComputeJacobianQQNLInertiaFunction\n");
  this->computeJacobianQQNLInertiaPtr = NULL;
  cShared.setFunction(&computeJacobianQQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->jacobianQQNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianQQNLInertiaFunction\n");
}

void LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction(string pluginPath, string functionName)
{
  IN("LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction\n");
  this->computeJacobianVelocityQNLInertiaPtr = NULL;
  cShared.setFunction(&computeJacobianVelocityQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->jacobianVelocityQNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction\n");
}


//////////////////////////


void LagrangianDS::fillDSWithDSXML()
{
  string plugin;

  IN("LagrangianDS::fillDSWithDSXML\n");
  if (this->dsxml != NULL)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(this->dsxml);

    this->ndof = lgptr->getNdof();
    this->n = 2 * this->ndof;

    this->r = SimpleVector(this->n);

    this->q = SimpleVector(this->ndof);
    this->q0 =  SimpleVector(this->ndof);
    this->qFree =  SimpleVector(this->ndof);

    this->velocity =  SimpleVector(this->ndof);
    this->velocity0 =  SimpleVector(this->ndof);
    this->velocityFree =  SimpleVector(this->ndof);

    this->p =  SimpleVector(this->ndof);

    this->fExt =  SimpleVector(this->ndof);
    this->fInt =  SimpleVector(this->ndof);
    this->QNLInertia =  SimpleVector(this->ndof);


    // Mass
    if ((static_cast <LagrangianDSXML*>(this->dsxml))->isMPlugin())
    {
      plugin = lgptr->getMPlugin();
      this->setComputeMassFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    // \warning : VA:  It is a very good idea to take the constant Mass Matrix, but for the moment a constant
    //  Mass Matrix is only read by a LagrangianLinearTIDS
    else this->mass = lgptr->getMMatrix();

    this->q0 = lgptr->getQ0();
    if (lgptr->hasQ()) this->q = lgptr->getQ();
    else
    {
      this->q = this->q0;
      cout << "Warning : q is not defined in the XML \n q is initialized with q0" << endl;
    }

    if (lgptr->hasQMemory()) this->qMemory = SiconosMemory::SiconosMemory(lgptr->getQMemoryXML());    //lgptr->getQMemory();
    else cout << "Warning : qMemory is not defined in the XML " << endl;

    this->velocity0 = lgptr->getVelocity0();
    if (lgptr->hasVelocity()) this->velocity = lgptr->getVelocity();
    else
    {
      this->velocity = this->velocity0;
      cout << "Warning : velocity is not defined in the XML \n velocity is initialized with velocity0" << endl;
    }

    if (lgptr->hasVelocityMemory()) this->velocityMemory = SiconosMemory::SiconosMemory(lgptr->getVelocityMemoryXML());    //lgptr->getVelocityMemory();
    else cout << "Warning : velocityMemory is not defined in the XML " << endl;

    static_cast<CompositeVector*>(this->x)->add(this->q);
    static_cast<CompositeVector*>(this->x)->add(this->velocity);

    static_cast<CompositeVector*>(this->x0)->add(this->q0);
    static_cast<CompositeVector*>(this->x0)->add(this->velocity0);

    // QNLInertia
    if (this->DSType == LNLDS)
    {
      // FInt
      if (lgptr->hasFint())
      {
        if (lgptr->isFintPlugin())
        {
          plugin = lgptr->getFintPlugin();
          this->setComputeFIntFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
        }
        else this->fInt = lgptr->getFintVector();
      }
      else
      {
        this->fInt = /*SiconosVector*/SimpleVector::SimpleVector();
        cout << "Warning : Fint is not defined in this LagrangianDS ( " << this->DSType << " )." << endl;
      }

      // JacobianQFInt
      if (lgptr->isJacobianQFintPlugin())
      {
        plugin = lgptr->getJacobianQFintPlugin();
        this->setComputeJacobianQFIntFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
      }
      else this->jacobianQFInt = lgptr->getJacobianQFintMatrix();

      // JacobianVelocityFInt
      if (lgptr->isJacobianVelocityFintPlugin())
      {
        plugin = lgptr->getJacobianVelocityFintPlugin();
        this->setComputeJacobianVelocityFIntFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
      }
      else this->jacobianVelocityFInt = lgptr->getJacobianVelocityFintMatrix();

      // JacobianQQNLInertia
      if (lgptr->isJacobianQQNLInertiaPlugin())
      {
        plugin = lgptr->getJacobianQQNLInertiaPlugin();
        this->setComputeJacobianQQNLInertiaFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
      }
      else this->jacobianQQNLInertia = lgptr->getJacobianQQNLInertiaMatrix();

      // JacobianVelocityQNLInertiaFunction
      if (lgptr->isJacobianVelocityQNLInertiaPlugin())
      {
        plugin = lgptr->getJacobianVelocityQNLInertiaPlugin();
        this->setComputeJacobianVelocityQNLInertiaFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
      }
      else this->jacobianVelocityQNLInertia = lgptr->getJacobianVelocityQNLInertiaMatrix();

      if (lgptr->isQNLInertiaPlugin())
      {
        plugin = lgptr->getQNLInertiaPlugin();
        this->setComputeQNLInertiaFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
      }
      else this->QNLInertia = lgptr->getQNLInertiaVector();
    }

    // FExt
    if (lgptr->isFextPlugin())
    {
      plugin = lgptr->getFextPlugin();
      this->setComputeFExtFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));

    }
    else
    {
      this->fExt = lgptr->getFextVector();
    }

    /****************************/

    this->number = this->dsxml->getNumber();

    if (this->dsxml->hasId() == true) this->id = this->dsxml->getId();
    else cout << "Warning : Id is not defined in the XML " << endl;

    if (this->dsxml->hasX0() == true) *(this->x0) = this->dsxml->getX0();
    else cout << "Warning : x0 is not defined in the XML " << endl;

    if (this->dsxml->hasX() == true)
    {
      *(this->x) = this->dsxml->getX();
    }
    else cout << "Warning : x is not defined in the XML " << endl;

    if (this->dsxml->hasXDot() == true)(this->xDot) = this->dsxml->getXDot();
    else cout << "Warning : xDot is not defined in the XML " << endl;

    if (this->dsxml->hasXMemory() == true) this->xMemory = SiconosMemory::SiconosMemory(this->dsxml->getXMemoryXML()); //this->dsxml->getXMemory();
    else cout << "Warning : xMemory is not defined in the XML " << endl;

    if (this->dsxml->hasXDotMemory() == true) this->xDotMemory = SiconosMemory::SiconosMemory(this->dsxml->getXDotMemoryXML()); //this->dsxml->getXDotMemory();
    else cout << "Warning : xDotMemory is not defined in the XML " << endl;

    if (this->dsxml->hasStepsInMemory() == true) this->stepsInMemory = this->dsxml->getStepsInMemory();
    else cout << "Warning : stepsInMemory is not defined in the XML " << endl;
  }
  else RuntimeException::selfThrow("LagrangianDS::fillDSWithDSXML - object DSXML does not exist");
  OUT("LagrangianDS::fillDSWithDSXML\n");
}

void LagrangianDS::display() const
{
  IN("LagrangianDS::display\n");

  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the LagrangianDS " << endl;
  DynamicalSystem::display();
  cout << "| ndof : " << this->ndof << endl;
  cout << "| q " << endl;
  this->q.display();
  cout << "| q0 " << endl;
  this->q0.display();
  cout << "| qFree " << endl;
  this->qFree.display();
  cout << "| velocity " << endl;
  this->velocity.display();
  cout << "| velocity0 " << endl;
  this->velocity0.display();
  cout << "| velocityFree " << endl;
  this->velocityFree.display();
  cout << "| p " << endl;
  this->p.display();
  cout << "-----------------------------------------------------" << endl << endl;

  OUT("LagrangianDS::display\n");
}


void LagrangianDS::init()
{
  IN("LagrangianDS::init\n");

  this->r = SimpleVector::SimpleVector();
  this->BC = NULL;
  this->jacobianX = SiconosMatrix::SiconosMatrix();
  this->setVectorFieldFunction("BasicPlugin.so", "vectorField");
  this->setComputeJacobianXFunction("BasicPlugin.so", "computeJacobianX");

  velocityFree = SimpleVector::SimpleVector();
  qFree = SimpleVector::SimpleVector();

  this->x = new CompositeVector();
  this->x0 = new CompositeVector();
  this->xFree = new CompositeVector();

  this->setComputeMassFunction("BasicPlugin.so", "computeMass");
  this->setComputeFIntFunction("BasicPlugin.so", "computeFInt");
  this->setComputeFExtFunction("BasicPlugin.so", "computeFExt");
  /*
   * \WARNING bizarre que quelque chose propre au LTIDS soit ici
   */
  if (this->DSType == LNLDS)
  {
    this->setComputeQNLInertiaFunction("BasicPlugin.so", "computeQNLInertia");
    this->setComputeJacobianQFIntFunction("BasicPlugin.so", "computeJacobianQFInt");
    this->setComputeJacobianVelocityFIntFunction("BasicPlugin.so", "computeJacobianVelocityFInt");
    this->setComputeJacobianQQNLInertiaFunction("BasicPlugin.so", "computeJacobianQQNLInertia");
    this->setComputeJacobianVelocityQNLInertiaFunction("BasicPlugin.so", "computeJacobianVelocityQNLInertia");
  }

  OUT("LagrangianDS::init\n");
}

void LagrangianDS::saveDSToXML()
{
  IN("LagrangianDS::saveDSToXML\n");
  DynamicalSystem::saveDSToXML();

  if (this->dsxml != NULL)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(this->dsxml);
    lgptr->setNdof(this->ndof);

    if (this->DSType == LTIDS)
    {
      lgptr->setMMatrix(&(this->mass));
    }
    else if (this->DSType == LNLDS)
    {
      lgptr->setMPlugin(this->massFunctionName);
    }

    lgptr->setQ(&(this->q));
    lgptr->setQ0(&(this->q0));
    lgptr->setQMemory(&(this->qMemory));

    lgptr->setVelocity(&(this->velocity));
    lgptr->setVelocity0(&(this->velocity0));
    lgptr->setVelocityMemory(&(this->velocityMemory));

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
      {
        lgptr->setFextVector(&(this->fExt));
      }
    }
    else
    {
      lgptr->setFextPlugin(this->fExtFunctionName);
    }

    if (this->DSType != LTIDS)  // for a LagrangianLinearTIDS, these plugin must not be saved
    {
      // FInt
      if (lgptr->hasFint())
      {
        if (!lgptr->isFintPlugin())
        {
          if (this->fInt.size() > 0)
            lgptr->setFintVector(&(this->fInt));
          else cout << "Warning : Fint can't be saved, the Fint vector is not defined." << endl;
        }
      }
      else
      {
        lgptr->setFintPlugin(this->fIntFunctionName);
      }

      // JacobianQFInt
      if (lgptr->hasJacobianQFint())
      {
        if (!lgptr->isJacobianQFintPlugin())
        {
          lgptr->setJacobianQFintMatrix(&(this->jacobianQFInt));
        }
      }
      else
      {
        lgptr->setJacobianQFintPlugin(this->jacobianQFIntFunctionName);
      }

      // JacobianVelocityFInt
      if (lgptr->hasJacobianVelocityFint())
      {
        if (!lgptr->isJacobianVelocityFintPlugin())
        {
          lgptr->setJacobianVelocityFintMatrix(&(this->jacobianVelocityFInt));
        }
      }
      else
      {
        lgptr->setJacobianVelocityFintPlugin(this->jacobianVelocityFIntFunctionName);
      }

      // JacobianQQNLInertia
      if (lgptr->hasJacobianQQNLInertia())
      {
        if (!lgptr->isJacobianQQNLInertiaPlugin())
        {
          lgptr->setJacobianQQNLInertiaMatrix(&(this->jacobianQQNLInertia));
        }
      }
      else
      {
        lgptr->setJacobianQQNLInertiaPlugin(this->jacobianQQNLInertiaFunctionName);
      }

      // JacobianVelocityQNLInertiaFunction
      if (lgptr->hasJacobianVelocityQNLInertia())
      {
        if (!lgptr->isJacobianVelocityQNLInertiaPlugin())
        {
          lgptr->setJacobianVelocityQNLInertiaMatrix(&(this->jacobianVelocityQNLInertia));
        }
      }
      else
      {
        lgptr->setJacobianVelocityQNLInertiaPlugin(this->jacobianVelocityQNLInertiaFunctionName);
      }

      // QNLInertia
      if (lgptr->hasQNLInertia())
      {
        if (!lgptr->isQNLInertiaPlugin())
        {
          lgptr->setQNLInertiaVector(&(this->QNLInertia));
        }
      }
      else
      {
        lgptr->setQNLInertiaPlugin(this->QNLInertiaFunctionName);
      }
    }
  }
  else RuntimeException::selfThrow("LagrangianDS::saveDSToXML - object DSXML does not exist");
  OUT("LagrangianDS::saveDSToXML\n");
}

void LagrangianDS::SiconosVectorSizeInit()
{
  /*
   * initilaisation of the SiconosVectors size
   */
  this->r = SimpleVector(this->n);

  this->q = SimpleVector(this->ndof);
  this->q0 =  SimpleVector(this->ndof);
  this->qFree =  SimpleVector(this->ndof);

  this->velocity =  SimpleVector(this->ndof);
  this->velocity0 =  SimpleVector(this->ndof);
  this->velocityFree =  SimpleVector(this->ndof);

  this->p =  SimpleVector(this->ndof);

  this->fExt =  SimpleVector(this->ndof);
  this->fInt =  SimpleVector(this->ndof);
  this->QNLInertia =  SimpleVector(this->ndof);
}

void LagrangianDS::CompositeVectorInit()
{
  /*
   * initialisation of the CompositeVectors
   */
  cout << "LagrangianDS::CompositeVectorInit()" << endl;
  static_cast<CompositeVector*>(this->x)->add(this->q);
  static_cast<CompositeVector*>(this->x)->add(this->velocity);

  static_cast<CompositeVector*>(this->x0)->add(this->q0);
  static_cast<CompositeVector*>(this->x0)->add(this->velocity0);
  cout << "/LagrangianDS::CompositeVectorInit()" << endl;
}


LagrangianDS* LagrangianDS::convert(DynamicalSystem* ds)
{
  cout << "LagrangianDS::convert (DynamicalSystem* ds)" << endl;
  LagrangianDS* lnlds = dynamic_cast<LagrangianDS*>(ds);
  return lnlds;
}

