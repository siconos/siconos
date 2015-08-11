/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

/*! \file Question.hpp
  \brief ask something to a class with a visitor */

/** example :
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

/* use boost array for the initialization of non const reference */
//#include <boost/type_traits.hpp>

#ifndef Question_hpp
#define Question_hpp

#include "SiconosVisitor.hpp"

#if (__cplusplus >= 201103L) && !defined(USE_BOOST_FOR_CXX11)
#include <type_traits>
#include <array>
#else
#include <boost/type_traits/remove_reference.hpp>
#include <boost/array.hpp>
#endif

/** a generic return value visitor */
template <class AnswerType>
struct Question : public SiconosVisitor
{
  typedef AnswerType type;
  type answer;

  Question() : answer(std11::array<typename std11::remove_reference<AnswerType>::type, 1>()[0])
  {};
  Question(AnswerType ref) : answer(ref) {};

};


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

}

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

}

/** apply a SiconosVisitor to a visitable object
 * \param v a visitable object
 */
template <class Visitor, class Visitable>
void apply(const Visitable& v)
{
  static Visitor t;

  v.accept(t);

}

/** apply a parameterized SiconosVisitor to a visitable object
 * \param v a visitable object
 * \param arg the SiconosVisitor argument
 */
template <class VisitorWithArgument, class Visitable, class Argument>
void apply(const Visitable& v, const Argument& arg)
{
  VisitorWithArgument t(arg);

  v.accept(t);

}

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

}



#define ANSWER(T,CODE)                          \
  void visit(const T& ds)                       \
  {                                             \
    answer = ds . CODE;                         \
  }

#define ANSWER_V(T,CODE)                        \
  void visit(const T& ds)                       \
  {                                             \
    answer = CODE;                              \
  }

#define ANSWER_F(T,CODE)                        \
  void visit(const T& ds)                       \
  {                                             \
    answer = CODE(ds);                          \
  }


#endif
