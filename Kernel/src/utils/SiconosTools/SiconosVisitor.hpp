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

/*! \file SiconosVisitor.hpp
  \brief A general visitor interface for siconos objects.
*/


#ifndef SiconosVisitor_hpp
#define SiconosVisitor_hpp

#include "RuntimeException.h"

/* objects that may be visited (1) */
class DynamicalSystem;
class Sphere;
class Disk;
class Circle;

/** This define a visitor pattern

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) June 14, 2009

   User must define a derivation of Siconos::Visitor class "myvisitor"
   with the correct visit() functions.

   Then the visitor can be used as :

   A_visitable_Siconos_Object->accept(Siconos::Visitor myvisitor)

   The visitor itself may also be visited in order to compound visits on
   several objects.

*/

//namespace Siconos
//{

#define QUOTE(M) #M

#define FAIL(X) \
  { RuntimeException::selfThrow                                         \
      (QUOTE(you must define a visit function for SP :: X in a derived class of SiconosVisitor)); }


#define VISIT(X)                                                    \
  virtual void visit(boost::shared_ptr<X>) FAIL(X);                 \
  virtual void visit(X&) FAIL(X);

class SiconosVisitor
{
public:

  /* idem (1) */
  VISIT(DynamicalSystem);
  VISIT(Sphere);
  VISIT(Disk);
  VISIT(Circle);

};

//}

TYPEDEF_SPTR(SiconosVisitor);


#undef VISIT
#undef FAIL
#undef QUOTE

#endif /* SiconosVisitor_hpp */
