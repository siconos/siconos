/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "Relation.h"

using namespace std;

// Default constructor
Relation::Relation(RELATION::TYPES newType,
                   RELATION::SUBTYPES newSub):
  relationType(newType), subType(newSub), hName("unamed"), gName("unamed")
{
  zeroPlugin();
}

// xml constructor
Relation::Relation(SP::RelationXML relxml,
                   RELATION::TYPES newType,
                   RELATION::SUBTYPES newSub):
  relationType(newType), subType(newSub),
  relationxml(relxml), hName("unamed"), gName("unamed")
{
  zeroPlugin();
  if (! relationxml)
    RuntimeException::selfThrow("Relation::fillRelationWithRelationXML - object RelationXML does not exist");
}
void Relation::zeroPlugin()
{
  pluginH.reset(new PluggedObject());
  pluginjXH.reset(new PluggedObject());
  pluginjLH.reset(new PluggedObject());
  pluginG.reset(new PluggedObject());
  pluginjLG.reset(new PluggedObject());
  pluginf.reset(new PluggedObject());
  plugine.reset(new PluggedObject());
}

void Relation::initializeMemory()
{
  mResiduy.reset(new BlockVector());
  unsigned int nslawSize = getInteractionPtr()->getNonSmoothLawPtr()->getNsLawSize();
  unsigned int numberOfRelations = getInteractionPtr()->getNumberOfRelations();
  for (unsigned int j = 0; j < numberOfRelations ; ++j)
    mResiduy->insertPtr(SP::SimpleVector(new SimpleVector(nslawSize)));

  mH_alpha.reset(new BlockVector());
  for (unsigned int j = 0; j < numberOfRelations; ++j)
    mH_alpha->insertPtr(SP::SimpleVector(new SimpleVector(nslawSize)));

}


Relation::~Relation()
{
}

void Relation::display() const
{
  cout << "=====> Relation of type "
       << relationType
       << " and subtype "
       << subType << endl;
  if (interaction.lock())
    cout
        << "- Interaction id"
        << interaction.lock()->getId()
        << endl;
  else cout << "- Linked interaction -> NULL" << endl;
}

// --- ResiduY functions
void Relation::computeResiduY(double t)
{
  //Residu_y = y_alpha_k+1 - H_alpha;
  *mResiduy = *mH_alpha;
  scal(-1, *mResiduy, *mResiduy);

  //      cout<<"Relation::computeResiduY mH_alpha"<<endl;
  //      mH_alpha->display();
  //      cout<<"Relation::computeResiduY Y"<<endl;
  //      getInteractionPtr()->getYPtr(0)->display();

  (*mResiduy) += *(getInteractionPtr()->getYPtr(0));

  //      cout<<" Relation::computeResiduY residuY"<<endl;
  //      mResiduy->display();

}
void Relation::computeG(double t)
{
  computeInput(t);
}

void Relation::setComputeJacLHFunction(const std::string& pluginPath, const std::string& functionName)
{
  pluginjLH->setComputeFunction(pluginPath, functionName);
  //  Plugin::setFunction(&pluginjLH, pluginPath, functionName);
  //    SSL::buildPluginName(pluginNamejLOutput,pluginPath,functionName);
}
void Relation::setComputeJacXHFunction(const std::string& pluginPath, const std::string& functionName)
{
  pluginjXH->setComputeFunction(pluginPath, functionName);
  //    Plugin::setFunction(&pluginjXH, pluginPath, functionName);
}

/** To set a plug-in function to compute input function g
 *  \param string : the complete path to the plugin
 *  \param string : the function name to use in this plugin
 */
void Relation::setComputeGFunction(const std::string& pluginPath, const std::string& functionName)
{
  pluginG->setComputeFunction(pluginPath, functionName);
  //    SSL::buildPluginName(pluginNameInput,pluginPath,functionName);
}
void Relation::setComputeFFunction(const std::string& pluginPath, const std::string& functionName)
{
  //  Plugin::setFunction(&pluginf, pluginPath, functionName,gName);
  pluginf->setComputeFunction(pluginPath, functionName);
  //    SSL::buildPluginName(pluginNamefplugin,pluginPath,functionName);
}
void Relation::setComputeEFunction(const std::string& pluginPath, const std::string& functionName)
{
  plugine->setComputeFunction(pluginPath, functionName);
  //Plugin::setFunction(&plugine, pluginPath, functionName,gName);
  //    SSL::buildPluginName(pluginNameeplugin,pluginPath,functionName);
}

/** To set a plug-in function to compute the jacobian according to x of the input
 *  \param string : the complete path to the plugin
 *  \param string : the function name to use in this plugin
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */
void Relation::setComputeJacLGFunction(const std::string& pluginPath, const std::string& functionName)
{
  pluginjLG->setComputeFunction(pluginPath, functionName);
  //  Plugin::setFunction(&pluginjLG, pluginPath, functionName);
  //    SSL::buildPluginName(pluginNamejLOutput,pluginPath,functionName);
}

void Relation::setComputeHFunction(const std::string& pluginPath, const std::string& functionName)
{
  pluginH->setComputeFunction(pluginPath, functionName);
  //  Plugin::setFunction(&pluginH, pluginPath, functionName,hName);
  //    SSL::buildPluginName(pluginNameOutput,pluginPath,functionName);
}
void Relation::saveRelationToXML() const
{
  RuntimeException::selfThrow("Relation - saveRelationToXML: not yet implemented for relation of type " + getType());
}
