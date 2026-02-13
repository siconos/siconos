/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
 * struct ForMass : public Question<std::shared_ptr<SiconosMatrix>>
 * {
 *    void visit(const LagrangianDS& ds)
 *    {
 *       Question::answer = ds.mass();
 *    }
 *
 * }
 * std::shared_ptr<DynamicalSystem> ds;
 * [...]
 *
 * auto mass = ask<ForMass>(*ds);
 */

#ifndef Question_hpp
#define Question_hpp

#include <array>
#include <optional>

#include "NSLVisitor.hpp"

namespace siconos::modeling::nonsmooth_laws {

/** a generic return value visitor */
template <typename AnswerType>
struct Question : public Visitor {
  using type = AnswerType;

  //  type answer{std::a rray<typename std::remove_reference<AnswerType>::type, 1>{}[0]};
  std::optional<type> answer;
  Question() = default;

  Question(AnswerType ref) : answer(ref){};
};

/** get some value from a visitable object with the help of a
    GeneralQuestion
    \param v a visitable object
 */
template <class GeneralQuestion, class Visitable>
typename GeneralQuestion::type ask(const Visitable& v) {
  GeneralQuestion t;

  v.accept(t);

  return *t.answer;
}

/** get some value from a visitable object with the help of a
    parameterized GeneralQuestion
    \param v a visitable object
    \param arg the GeneralQuestion argument
 */
template <class GeneralQuestion, class Visitable, class Argument>
typename GeneralQuestion::type ask(const Visitable& v, const Argument& arg) {
  GeneralQuestion t(arg);

  v.accept(t);

  return *t.answer;
}

/** apply a siconos::modeling::nonsmooth_laws::Visitor to a visitable object
 * \param v a visitable object
 */
template <class Vtor, class Visitable>
void apply(const Visitable& v) {
  static Vtor t;

  v.accept(t);
}
}  // namespace siconos::modeling::nonsmooth_laws
#endif
