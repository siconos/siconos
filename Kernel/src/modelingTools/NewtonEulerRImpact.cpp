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


#include "NewtonEulerRImpact.hpp"
#include <boost/math/quaternion.hpp>
using namespace std;
void NewtonEulerRImpact::computeJachq(double t)
{
  DSIterator itDS = interaction()->dynamicalSystemsBegin();
  SP::DynamicalSystem aux = *itDS;
  //assert (&(*aux)==&(*_ds1));
  itDS++;

  bool has2Bodies = false;
  if (itDS != interaction()->dynamicalSystemsEnd())
    has2Bodies = true;
  _jachq->setValue(0, 0, _Nc->getValue(0));
  _jachq->setValue(0, 1, _Nc->getValue(1));
  _jachq->setValue(0, 2, _Nc->getValue(2));
  if (has2Bodies)
  {
    _jachq->setValue(0, 7, -_Nc->getValue(0));
    _jachq->setValue(0, 8, -_Nc->getValue(1));
    _jachq->setValue(0, 9, -_Nc->getValue(2));
  }
  SP::BlockVector BlockX = boost::static_pointer_cast<BlockVector>((data[q0]));
  for (int iDS = 0; iDS < 2; iDS++)
  {
    if (!has2Bodies && iDS == 1)
      continue;
    double sign = 1.0;
    SP::SiconosVector q = (BlockX->getAllVect())[iDS];
    //      printf("ds%d->q :",iDS);q->display();
    ::boost::math::quaternion<float>    quatGP;
    if (iDS == 0)
    {
      ::boost::math::quaternion<float>    quatAux(0, _Pc1->getValue(0) - q->getValue(0), _Pc1->getValue(1) - q->getValue(1),
          _Pc1->getValue(2) - q->getValue(2));
      quatGP = quatAux;
    }
    else
    {
      sign = -1.0;
      ::boost::math::quaternion<float>    quatAux(0, _Pc2->getValue(0) - q->getValue(0), _Pc2->getValue(1) - q->getValue(1),
          _Pc2->getValue(2) - q->getValue(2));
      quatGP = quatAux;
    }
    //      printf("GP :%lf, %lf, %lf\n",quatGP.R_component_2(),quatGP.R_component_3(),quatGP.R_component_4());
    ::boost::math::quaternion<float>    quatQ(q->getValue(3), q->getValue(4), q->getValue(5), q->getValue(6));
    ::boost::math::quaternion<float>    quatcQ(q->getValue(3), -q->getValue(4), -q->getValue(5), -q->getValue(6));
    ::boost::math::quaternion<float>    quat0(1, 0, 0, 0);
    ::boost::math::quaternion<float>    quatBuff;
    quatBuff = (quatGP * quatQ) + (quatcQ * quatGP);
    _jachq->setValue(0, 7 * iDS + 3, sign * (quatBuff.R_component_2()*_Nc->getValue(0) +
                     quatBuff.R_component_3()*_Nc->getValue(1) + quatBuff.R_component_4()*_Nc->getValue(2)));
    for (int i = 1; i < 4; i++)
    {
      ::boost::math::quaternion<float>    quatei(0, (i == 1) ? 1 : 0, (i == 2) ? 1 : 0, (i == 3) ? 1 : 0);
      quatBuff = quatei * quatcQ * quatGP - quatGP * quatQ * quatei;
      _jachq->setValue(0, 7 * iDS + 3 + i, sign * (quatBuff.R_component_2()*_Nc->getValue(0) +
                       quatBuff.R_component_3()*_Nc->getValue(1) + quatBuff.R_component_4()*_Nc->getValue(2)));
    }
  }
  //    printf("computeJachq :");_jachq->display();
  //    printf("q1dot : ");sOCCContacts[_indexContact]._DS1->dotq()->display();
  //    printf("q2dot : ");sOCCContacts[_indexContact]._DS2->dotq()->display();

}
