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
  relationType(newType),
  subType(newSub), hPlugged(false), gPlugged(false),
  hName("unamed"), gName("unamed")
{}

// xml constructor
Relation::Relation(SP::RelationXML relxml,
                   RELATION::TYPES newType,
                   RELATION::SUBTYPES newSub):
  relationType(newType), subType(newSub),
  relationxml(relxml), hPlugged(false),
  gPlugged(false), hName("unamed"), gName("unamed")
{
  if (! relationxml)
    RuntimeException::selfThrow("Relation::fillRelationWithRelationXML - object RelationXML does not exist");
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

  //  cout<<"Relation::computeResiduY mH_alpha"<<endl;
  //  mH_alpha->display();
  //  cout<<"Relation::computeResiduY Y"<<endl;
  //  getInteractionPtr()->getYPtr(0)->display();

  (*mResiduy) += *(getInteractionPtr()->getYPtr(0));

  //  cout<<" Relation::computeResiduY residu Y"<<endl;
  //  mResiduy->display();

}
void Relation::computeG(double t)
{
  computeInput(t, 0);
}

SP::SiconosMatrix Relation::getCPtr()
{
  return getJacHPtr(0);
}
SP::SiconosMatrix Relation::getDPtr()
{
  return getJacHPtr(1);
}
SP::SiconosMatrix Relation::getBPtr()
{
  return getJacGPtr(0);
}


void Relation::saveRelationToXML() const
{
  RuntimeException::selfThrow("Relation - saveRelationToXML: not yet implemented for relation of type " + getType());
}
