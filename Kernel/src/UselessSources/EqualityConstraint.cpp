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
#include "EqualityConstraint.h"
using namespace std;

EqualityConstraint::EqualityConstraint():
  type(NLINEAREC), number(0), id("none"), ecXML(NULL), computeInputName("none"),
  computeOutputName("none"), computeOutputPtr(NULL), computeInputPtr(NULL)
{}

EqualityConstraint::EqualityConstraint(EqualityConstraintXML* newEcXML):
  type(NLINEAREC), number(0), id("none"), ecXML(newEcXML), computeInputName("none"),
  computeOutputName("none"), computeOutputPtr(NULL), computeInputPtr(NULL)
{}

EqualityConstraint::~EqualityConstraint()
{
  if (G != NULL) delete G;
}

std::vector<DSInputOutput*> EqualityConstraint::getDSInputOutputs(void)
{
  return dsioVector;
}

DSInputOutput* EqualityConstraint::getDSInputOutput(const unsigned int& i)
{
  if (i >= dsioVector.size())
    RuntimeException::selfThrow("EqualityConstraint - getDSInputOutput : \'i\' is out of range");
  return dsioVector[i];
}

void EqualityConstraint::setDSInputOutputs(std::vector<DSInputOutput*> dsioVect)
{
  dsioVector = dsioVect;
}

void EqualityConstraint::addDSInputOutput(DSInputOutput* dsio)
{
  //  DSInputOutput* dsioTmp;
  //  dsioTmp = new DSInputOutput();
  //  *dsioTmp = *dsio;

  /*
   *  in EqualityConstraint class, we don't create new objects in the DSInputOutput vector
   *    => we only save a link (pointer) on the DSInputOutputs of the DynamicalSystems !!
   */
  dsioVector.push_back(dsio);
}

void EqualityConstraint::saveEqualityConstraintToXML()
{
  if (ecXML != NULL)
  {
    /*
     * these attributes are only required for LagrangianNonLinear DSInputOutput !
     */
    //    disoxml->setComputeInputPlugin( computeInputName );
    //    dsioxml->setComputeOutputPlugin( computeOutputName );
    ecXML->setG(G);
  }
  else RuntimeException::selfThrow("EqualityConstraint::saveEqualityConstraintToXML - object EqualityConstraintXML does not exist");
}

void EqualityConstraint::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the EqualityConstraint " << endl;
  cout << "| id : " << id << endl;
  cout << "| number : " << number << endl;
  cout << "| G : " << endl;
  G->display();
  cout << "-----------------------------------------------------" << endl << endl;
}

void EqualityConstraint::init()
{
  number = 0;
  id = "none";
  ecXML = NULL;
}

////////////////////////////////
void EqualityConstraint::computeOutput(double time)
{
  if (computeOutputPtr == NULL) RuntimeException::selfThrow("computeOutput() is not linked to a plugin function");

  //to do
  //computeOutputPtr(&x(0), &time, &lambdaPtr(0), &y(0));
  //  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
  //
  //  DynamicalSystem *ds1 ,*ds2;
  //  SiconosVector *y = interaction->getYPtr();
  //  SiconosVector *yDot = interaction->getYDotPtr();
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
  //        //*y = (h * q) + b;
  //
  //      BlockVector vel;
  //      vel.add(*(d1->getVelocityPtr()));
  //      vel.add(*(d2->getVelocityPtr()));
  //      *yDot = (h * vel);
  //
  //      computeOutputPtr(*q, 0.0, lambda, y);
  //      }
  //    else
  //    {
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //    }
}

void EqualityConstraint::computeInput(double time)
{
  if (computeInputPtr == NULL) RuntimeException::selfThrow("computeInput() is not linked to a plugin function");

  //to do
  //computeInputPtr(&x(0), &time, &lambdaPtr(0), &r(0));
}

void EqualityConstraint::setComputeOutputFunction(std::string pluginPath, std::string functionName)
{
  computeOutputPtr = NULL;
  cShared.setFunction(&computeOutputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeOutputName = plugin + ":" + functionName;
}

void EqualityConstraint::setComputeInputFunction(std::string pluginPath, std::string functionName)
{
  computeInputPtr = NULL;
  cShared.setFunction(&computeInputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeInputName = plugin + ":" + functionName;
}
///////////////////////////////

void EqualityConstraint::fillEqualityConstraintWithEqualityConstraintXML()
{
  if (ecXML != NULL)
  {
    string plugin;
    // computeInput
    if (ecXML->hasComputeInput())
    {
      cout << "EqualityConstraintPluginType == " << type << endl;
      plugin = (ecXML)->getComputeInputPlugin();
      setComputeInputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning - No computeInput method is defined in a EqualityConstraint " << getType() << endl;

    // computeOutput
    if (ecXML->hasComputeOutput())
    {
      cout << "EqualityConstraintPluginType == " << type << endl;
      plugin = (ecXML)->getComputeOutputPlugin();
      setComputeOutputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning - No computeOutput method is defined in a Relation " << getType() << endl;

    number = ecXML->getNumber();
    G = new SimpleMatrix(ecXML->getG());
  }
  //else RuntimeException::selfThrow("EqualityConstraint::fillEqualityConstraintWithEqualityConstraintXML - object EqualityConstraintXML does not exist");
}

void EqualityConstraint::createEqualityConstraint(EqualityConstraintXML *newEcXML,
    int newNumber, SiconosMatrix *newG,
    std::vector<DSInputOutput*> *newDsioVector)
{
  if (ecXML != NULL)
  {
    ecXML = newEcXML;
    type = NLINEAREC;
    fillEqualityConstraintWithEqualityConstraintXML();
  }
  else
  {
    ecXML = NULL;
    type = NLINEAREC;
    number = newNumber;
    G = new SimpleMatrix(*newG);
    dsioVector = *newDsioVector;
  }
}

