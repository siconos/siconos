/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

   SiconosVisitor also define a type visitor object under the
   namespace Type:: and some functions to access type of visitables
   classes:

   Type::value(c) : the type of the visitable object c as an enum.

   Type::name(c)  : the name of the Type::value as a std::string

*/

#include "SiconosPointers.hpp"
#include "RuntimeException.hpp"

/* all Siconos classes that may be visited are defined there */
#include "SiconosVisitables.hpp"

/* convenient macros */
#define SICONOS_VISITOR_QUOTE(M) #M

#define SICONOS_VISITOR_FAIL(X)                                         \
  { RuntimeException::selfThrow                                         \
      ( SICONOS_VISITOR_QUOTE(you must define a visit function for X in a derived class of SiconosVisitor)); }


/** hook to be inserted in a virtual class definiton */
#define VIRTUAL_ACCEPT_VISITORS(FROMCLASS)                              \
  template<typename Archive> friend class SiconosSerializer;            \
  virtual void acceptSP(SP::SiconosVisitor)                             \
  {                                                                     \
    RuntimeException::selfThrow                                         \
      ( SICONOS_VISITOR_QUOTE(this class derived from FROMCLASS does not accept a visitor for shared pointers)); \
  };                                                                    \
  virtual void accept(SiconosVisitor&) const                            \
  {                                                                     \
    RuntimeException::selfThrow                                         \
      ( "accept: no visitor defined");                     \
  };                                                                    \
  virtual void acceptSerializer(SiconosVisitor&)                        \
  {                                                                     \
    RuntimeException::selfThrow                                         \
      ( "acceptSerializer: no serializer defined");                                           \
  };                                                                    \
  virtual inline Type::Siconos acceptType(FindType& ft) const           \
  { RuntimeException::selfThrow                                         \
      ( SICONOS_VISITOR_QUOTE(this class derived from FROMCLASS does not accept a type visitor)); \
    return Type::void_type;                                             \
  }                                                                     \
 
/** hooks to be inserted in class definition */
#define ACCEPT_STD_VISITORS()                                           \
  template<typename Archive> friend class SiconosSerializer;            \
  virtual void accept(SiconosVisitor& tourist) const { tourist.visit(*this); } \
  virtual void acceptSerializer(SiconosVisitor& serializer) { serializer.visit(*this); } \
  virtual inline Type::Siconos acceptType(FindType& ft) const { return ft.visit(*this); } \
 
#define ACCEPT_SP_VISITORS()                                            \
  virtual void acceptSP(SP::SiconosVisitor tourist) { tourist->visit(shared_from_this()); }

#define ACCEPT_VISITORS()                       \
  ACCEPT_SP_VISITORS()                          \
  ACCEPT_STD_VISITORS()                         \
 




/* objects that may be visited (1) */
#undef REGISTER
#undef REGISTER_STRUCT

#define REGISTER(X) class X;
#define REGISTER_STRUCT(X) struct X;

#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y) REGISTER(X)
#define REGISTER_BASE_EXTERN(X,Y) REGISTER(X)

SICONOS_VISITABLES()

/* associated types */
#undef REGISTER
#undef REGISTER_STRUCT
#define REGISTER(X) X,

#define REGISTER_STRUCT(X) X,

#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y) REGISTER(X)
#define REGISTER_BASE_EXTERN(X,Y) REGISTER(X)
namespace Type
{
enum Siconos
{
  SICONOS_VISITABLES()
  void_type
};
}

//namespace Siconos
//{


/* the type visitor */
#undef REGISTER
#define REGISTER(X)                                                 \
  virtual Type::Siconos visit(const X&) const { return Type::X; };  \
 
#undef REGISTER_STRUCT
#define REGISTER_STRUCT(X) REGISTER(X)

#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y)                                         \
  virtual Type::Siconos visit(const X&) const { return Type::Y; }; \
 
#define REGISTER_BASE_EXTERN(X,Y) REGISTER_BASE(X,Y)

struct FindType
{
  SICONOS_VISITABLES()
};

/* the base visitor */
#undef REGISTER
#define REGISTER(X)             \
  virtual void visit(boost::shared_ptr<X>) SICONOS_VISITOR_FAIL(SP :: X); \
  virtual void visit(X&) SICONOS_VISITOR_FAIL(X);                         \
  virtual void visit(const X&) SICONOS_VISITOR_FAIL(X);

#undef REGISTER_STRUCT
#define REGISTER_STRUCT(X) REGISTER(X)

#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y) REGISTER(X)

#define REGISTER_BASE_EXTERN(X,Y) REGISTER_BASE(X,Y)

struct SiconosVisitor
{
  SICONOS_VISITABLES()
  virtual ~SiconosVisitor() {};
};

/* some functions in Type namespace */
namespace Type
{
static FindType find;

template <typename C>
inline Siconos value(const C& c)
{
  return c.acceptType(find);
}

#undef REGISTER
#define REGISTER(X) case Type:: X : r.reset(new std::string(#X)); break;

#undef REGISTER_STRUCT
#define REGISTER_STRUCT(X) REGISTER(X)
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y) REGISTER(X)

#define REGISTER_BASE_EXTERN(X,Y) REGISTER_BASE(X,Y)

namespace
{
boost::shared_ptr<std::string> str(const Siconos& X)
{
  boost::shared_ptr<std::string> r;

  switch (X)
  {
    SICONOS_VISITABLES()
  default:
    assert(0);
  }

  return(r);
}
}


template <class C>
std::string name(const C& c)
{
  return *(Type::str(Type::value(c)));
}

}




/* ask something to a class with a visitor */

/* example :
 *
 * struct ForMass : public Question<SP::SiconosMatrix>
 * {
 *    void visit(const LagrangianDS& ds)
 *    {
 *       Question::answer = ds.mass();
 *    }
 *
 * }
 * SP::DynamicalSystem ds
 * [...]
 *
 * SP::SiconosMatrix mass = ask<ForMass>(*ds);
 */


/** get some value from a visitable object with the help of a
    GeneralQuestion
    \param v a visitable object
 */
template <class GeneralQuestion, class Visitable>
typename GeneralQuestion::type ask(const Visitable& v)
{
  GeneralQuestion t;

  v.accept(t);

  return t.answer;

};

/** get some value from a visitable object with the help of a
    parameterized GeneralQuestion
    \param v a visitable object
    \param arg the GeneralQuestion argument
 */
template <class GeneralQuestion, class Visitable, class Argument>
typename GeneralQuestion::type ask(const Visitable& v, const Argument& arg)
{
  GeneralQuestion t(arg);

  v.accept(t);

  return t.answer;

};

/** apply a SiconosVisitor to a visitable object
 * \param v a visitable object
 */
template <class Visitor, class Visitable>
void apply(const Visitable& v)
{
  static Visitor t;

  v.accept(t);

};

/** apply a parameterized SiconosVisitor to a visitable object
 * \param v a visitable object
 * \param arg the SiconosVisitor argument
 */
template <class VisitorWithArgument, class Visitable, class Argument>
void apply(const Visitable& v, const Argument& arg)
{
  VisitorWithArgument t(arg);

  v.accept(t);

};

/** apply a parameterized SiconosVisitor to a visitable object
 * \param v a visitable object
 * \param arg1 the first SiconosVisitor argument
 * \param arg2 the second SiconosVisitor argument
 */
template <class VisitorWith2Arguments, class Visitable, class Argument1, class Argument2>
void apply(const Visitable& v, const Argument1& arg1, const Argument2& arg2)
{
  VisitorWith2Arguments t(arg1, arg2);

  v.accept(t);

};

/* use boost array for the initialization of non const reference */
#include <boost/type_traits.hpp>
#include <boost/array.hpp>

/** a generic return value visitor */
template <class AnswerType>
struct Question : public SiconosVisitor
{
  typedef AnswerType type;
  type answer;

  Question() : answer(boost::array<typename boost::remove_reference<AnswerType>::type, 1>()[0])
  {};
  Question(AnswerType ref) : answer(ref) {};

};

#define ANSWER(T,CODE) \
  void visit(const T& ds)                       \
  {                                             \
    answer = ds . CODE;                         \
  }

#define ANSWER_V(T,CODE) \
  void visit(const T& ds)                       \
  {                                             \
    answer = CODE;                              \
  }

#define ANSWER_F(T,CODE) \
  void visit(const T& ds)                       \
  {                                             \
    answer = CODE(ds);                          \
  }




TYPEDEF_SPTR(SiconosVisitor);

#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN

#endif /* SiconosVisitor_hpp */
