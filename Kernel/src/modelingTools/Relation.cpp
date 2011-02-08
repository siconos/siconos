/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "Relation.hpp"

using namespace std;

// Default constructor
Relation::Relation(RELATION::TYPES newType,
                   RELATION::SUBTYPES newSub):
  relationType(newType), subType(newSub)
{
  zeroPlugin();
}

// xml constructor
Relation::Relation(SP::RelationXML relxml,
                   RELATION::TYPES newType,
                   RELATION::SUBTYPES newSub):
  relationType(newType), subType(newSub),
  relationxml(relxml)
{
  zeroPlugin();
  if (! relationxml)
    RuntimeException::selfThrow("Relation::fillRelationWithRelationXML - object RelationXML does not exist");
}
void Relation::zeroPlugin()
{
  _pluginh.reset(new PluggedObject());
  _pluginJachx.reset(new PluggedObject());
  _pluginJachlambda.reset(new PluggedObject());
  _pluging.reset(new PluggedObject());
  _pluginJacLg.reset(new PluggedObject());
  _pluginf.reset(new PluggedObject());
  _plugine.reset(new PluggedObject());
}

void Relation::initializeMemory()
{
  _Residuy.reset(new BlockVector());
  unsigned int nslawSize = interaction()->nonSmoothLaw()->size();
  unsigned int numberOfRelations = interaction()->numberOfRelations();
  for (unsigned int j = 0; j < numberOfRelations ; ++j)
    _Residuy->insertPtr(SP::SimpleVector(new SimpleVector(nslawSize)));

  _h_alpha.reset(new BlockVector());
  for (unsigned int j = 0; j < numberOfRelations; ++j)
    _h_alpha->insertPtr(SP::SimpleVector(new SimpleVector(nslawSize)));

}




Relation::~Relation()
{
}


const std::string Relation::getJachxName() const
{
  if (_pluginJachx->fPtr)
    return _pluginJachx->getPluginName();
  return "unamed";
}

const std::string Relation::gethName() const
{
  if (_pluginh->fPtr)
    return _pluginh->getPluginName();
  return "unamed";

}
const std::string Relation::getgName() const
{
  if (_pluging->fPtr)
    return _pluging->getPluginName();
  return "unamed";

}

const std::string Relation::getJacgName(unsigned int) const
{
  return "unamed";
}


void Relation::display() const
{
  cout << "=====> Relation of type "
       << relationType
       << " and subtype "
       << subType << endl;
  if (_interaction.lock())
    cout
        << "- Interaction id"
        << _interaction.lock()->getId()
        << endl;
  else cout << "- Linked interaction -> NULL" << endl;
}

// --- ResiduY functions
void Relation::computeResiduY(double t)
{
  //Residu_y = y_alpha_k+1 - H_alpha;
  *_Residuy = *_h_alpha;
  scal(-1, *_Residuy, *_Residuy);

  //      cout<<"Relation::computeResiduY _h_alpha"<<endl;
  //      _h_alpha->display();
  //      cout<<"Relation::computeResiduY Y"<<endl;
  //      interaction()->y(0)->display();

  (*_Residuy) += *(interaction()->y(0));

  //      cout<<" Relation::computeResiduY residuY"<<endl;
  //      _Residuy->display();

}
void Relation::computeg(double t)
{
  unsigned int i = interaction()->getRelativeDegree();
  if (i)
    i--;
  computeInput(t, i);
}

void Relation::setComputeJachlambdaFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJachlambda->setComputeFunction(pluginPath, functionName);
  //  Plugin::setFunction(&_pluginJachlambda, pluginPath, functionName);
  //    SSL::buildPluginName(pluginNamejLOutput,pluginPath,functionName);
}
void Relation::setComputeJachxFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJachx->setComputeFunction(pluginPath, functionName);
  //    Plugin::setFunction(&_pluginJachx, pluginPath, functionName);
}

/** To set a plug-in function to compute input function g
 *  \param string : the complete path to the plugin
 *  \param string : the function name to use in this plugin
 */
void Relation::setComputegFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluging->setComputeFunction(pluginPath, functionName);
  //    SSL::buildPluginName(pluginNameInput,pluginPath,functionName);
}
void Relation::setComputeFFunction(const std::string& pluginPath, const std::string& functionName)
{
  //  Plugin::setFunction(&pluginf, pluginPath, functionName,gName);
  _pluginf->setComputeFunction(pluginPath, functionName);
  //    SSL::buildPluginName(pluginNamefplugin,pluginPath,functionName);
}
void Relation::setComputeEFunction(const std::string& pluginPath, const std::string& functionName)
{
  _plugine->setComputeFunction(pluginPath, functionName);
  //Plugin::setFunction(&_plugine, pluginPath, functionName,gName);
  //    SSL::buildPluginName(pluginNameeplugin,pluginPath,functionName);
}

/** To set a plug-in function to compute the jacobian according to x of the input
 *  \param string : the complete path to the plugin
 *  \param string : the function name to use in this plugin
 *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
 */
void Relation::setComputeJacglambdaFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginJacLg->setComputeFunction(pluginPath, functionName);
  //  Plugin::setFunction(&_pluginJacLg, pluginPath, functionName);
  //    SSL::buildPluginName(pluginNamejLOutput,pluginPath,functionName);
}

void Relation::setComputehFunction(const std::string& pluginPath, const std::string& functionName)
{
  _pluginh->setComputeFunction(pluginPath, functionName);
  //  Plugin::setFunction(&_pluginh, pluginPath, functionName,hName);
  //    SSL::buildPluginName(pluginNameOutput,pluginPath,functionName);
}
const  SP::SiconosVector Relation::residuR()
{
  assert(0);
  return data[0];
}
void Relation::saveRelationToXML() const
{
  RuntimeException::selfThrow("Relation - saveRelationToXML: not yet implemented for relation of type " + getType());
}
