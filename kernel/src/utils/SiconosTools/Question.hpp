/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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

#include <SiconosConfig.h>
#if defined(SICONOS_STD_ARRAY) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
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

#endif
