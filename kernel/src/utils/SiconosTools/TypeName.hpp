/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file TypeName.hpp
  \brief get a string name from visitable classes
*/

#ifndef TypeName_hpp
#define TypeName_hpp

#include <string>
#include <variant>
// #include "SiconosVisitables.hpp" // For SICONOS_VISITABLES macro

// #undef REGISTER
// #define REGISTER(X)                        \
//   case siconos::internal::type::Siconos:: X:	\
//     r = std::make_shared<std::string>(#X); \
//     break;

// #undef REGISTER_STRUCT
// #define REGISTER_STRUCT(X) REGISTER(X)

// namespace siconos::internal {

// struct SiconosVisitor;

// inline std::shared_ptr<std::string> str(const siconos::internal::type::Siconos& X)
// {
//   std::shared_ptr<std::string> r;

//   switch (X) {
//     SICONOS_VISITABLES()
//     default:
//       assert(0);
//   }

//   return (r);
// }

// template <class C>
// std::string name(const C& c)
// {
//   return *(siconos::internal::type::str(siconos::internal::type::value(c)));
// }



//   template <class T > inline std::shared_ptr<std::string> str(const siconos::internal::type::Siconos& X)
// {
//   std::shared_ptr<std::string> r;

//   switch (X) {
//     SICONOS_VISITABLES()
//     default:
//       assert(0);
//   }

//   return (r);
// }

  
// }  // namespace siconos::internal


namespace siconos::typesv{

  struct FindType{

    // auto operator()(const siconos::experimental::B& obj) const {return experimental::Type::B; }
    // auto operator()(const siconos::experimental::C& obj) const {return experimental::Type::C; }
    // auto operator()(std::shared_ptr<siconos::experimental::B> obj) const {return experimental::Type::B; }
    // auto operator()(std::shared_ptr<siconos::experimental::C> obj) const {return experimental::Type::C; }
  };
  
  // auto type_value = []( auto& obj) { return obj.name(); };

  static FindType find;

  template <typename C> auto type_value(const C &c) {return std::visit(find, c);}
  
  template <typename T> constexpr auto str(const T& X)
  {
  
  switch (X) {
  case T::B:
    //r = std::make_shared<std::string>("siconos::experimental::B");
    return "siconos::experimental::B";
    break;
  case T::C:
    //r = std::make_shared<std::string>("siconos::experimental::C");
    return "siconos::experimental::C";
    break;
  default:
    assert(0);
  }
  }

  template <class C> std::string type_name(const C& c)
  {
     return str(type_value(c));
  }

  
    
}


#endif
