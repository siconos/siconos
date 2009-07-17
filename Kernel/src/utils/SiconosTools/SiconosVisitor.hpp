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

/** A visitor pattern.

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) June 14, 2009

   User have to instantiate a derivation of SiconosVisitor class :

   struct myvisitor : public SiconosVisitor

   with some wanted visit() functions.

   Then the visitor can be used as :

   A_visitable_Siconos_Object->accept(Siconos::Visitor myvisitor)

*/

#include "RuntimeException.h"

/* objects that may be visited (1) */
class DynamicalSystem;
class Sphere;
class Disk;
class Circle;
class CircleCircleR;
class DiskDiskR;
class DiskPlanR;
class Sphere;
class SphereSphereR;
class SpherePlanR;
class NonSmoothLaw;
class MixedComplementarityConditionNSL;
class ComplementarityConditionNSL;
class NewtonImpactNSL;
class NewtonImpactFrictionNSL;

//namespace Siconos
//{

#define SICONOS_VISITOR_QUOTE(M) #M

#define SICONOS_VISITOR_FAIL(X)                                                         \
  { RuntimeException::selfThrow                                         \
      ( SICONOS_VISITOR_QUOTE(you must define a visit function for SP :: X in a derived class of SiconosVisitor)); }


/** hook to be inserted in a virtual class definiton */
#define VIRTUAL_ACCEPT_VISITORS(FROMCLASS)                              \
  virtual void accept(SP::SiconosVisitor)                               \
  { RuntimeException::selfThrow                                         \
      ( SICONOS_VISITOR_QUOTE(this class derived from FROMCLASS does not accept a visitor for shared pointers)); }; \
  virtual void accept(SiconosVisitor&)                                  \
  { RuntimeException::selfThrow                                         \
      ( SICONOS_VISITOR_QUOTE(this class derived from FROMCLASS does not accept a visitor)); }

/** hooks to be inserted in class definition */
#define ACCEPT_STD_VISITORS()                                           \
  virtual void accept(SiconosVisitor& tourist) { tourist.visit(*this); } \
 
#define ACCEPT_SP_VISITORS()                                            \
  virtual void accept(SP::SiconosVisitor tourist) { tourist->visit(shared_from_this()); }

#define ACCEPT_VISITORS() \
  ACCEPT_SP_VISITORS();   \
  ACCEPT_STD_VISITORS()


#define VISIT(X)                                                    \
  virtual void visit(boost::shared_ptr<X>) SICONOS_VISITOR_FAIL(X);                 \
  virtual void visit(X&) SICONOS_VISITOR_FAIL(X);

class SiconosVisitor
{
public:

  /* idem (1) */
  VISIT(DynamicalSystem);
  VISIT(Disk);
  VISIT(Circle);
  VISIT(DiskPlanR);
  VISIT(CircleCircleR);
  VISIT(DiskDiskR);
  VISIT(Sphere);
  VISIT(SphereSphereR);
  VISIT(SpherePlanR);
  VISIT(NonSmoothLaw);
  VISIT(MixedComplementarityConditionNSL);
  VISIT(ComplementarityConditionNSL);
  VISIT(NewtonImpactNSL);
  VISIT(NewtonImpactFrictionNSL);
};

//}


TYPEDEF_SPTR(SiconosVisitor);

#undef VISIT

#endif /* SiconosVisitor_hpp */
