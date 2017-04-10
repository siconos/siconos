/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*! \file SiconosVisitor.hpp
  \brief A general visitor interface for siconos objects.
*/


#ifndef SiconosVisitor_hpp
#define SiconosVisitor_hpp

/** A visitor pattern.

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) June 14, 2009

   User has to instantiate a derivation of SiconosVisitor class :

   struct myvisitor : public SiconosVisitor
   {
     void visit(const LagrangianDS& ds)
     {
       ...
     }
   }

   with some wanted visit() functions.

   Then the visitor may be used as :

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
      ( "accept: no visitor defined");                                  \
  };                                                                    \
  virtual void accept_modifier(SiconosVisitor&)                         \
  {                                                                     \
    RuntimeException::selfThrow                                         \
      ( "accept: no modifier defined");                                 \
  };                                                                    \
  virtual void acceptSerializer(SiconosVisitor&)                        \
  {                                                                     \
    RuntimeException::selfThrow                                         \
      ( "acceptSerializer: no serializer defined");                     \
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
  virtual void accept_modifier(SiconosVisitor& tourist) { tourist.visit(*this); } \
  virtual void acceptSerializer(SiconosVisitor& serializer) { serializer.visit(*this); } \
  virtual inline Type::Siconos acceptType(FindType& ft) const { return ft.visit(*this); } \

#define ACCEPT_NONVIRTUAL_VISITORS()                                    \
  template<typename Archive> friend class SiconosSerializer;            \
  void accept(SiconosVisitor& tourist) const { tourist.visit(*this); }  \
  void acceptSerializer(SiconosVisitor& serializer) { serializer.visit(*this); } \
  inline Type::Siconos acceptType(FindType& ft) const { return ft.visit(*this); } \

#define ACCEPT_SP_VISITORS()                                            \
  virtual void acceptSP(SP::SiconosVisitor tourist) { tourist->visit(shared_from_this()); }

#define ACCEPT_VISITORS()                       \
  ACCEPT_SP_VISITORS()                          \
  ACCEPT_STD_VISITORS()                         \

/** hooks to be inserted in class definition */
#define ACCEPT_BASE_STD_VISITORS(BASE)                                  \
  template<typename Archive> friend class SiconosSerializer;            \
  virtual void acceptBase(SiconosVisitor& tourist) const { tourist.visit(*static_cast<const BASE *>(this)); } \
  virtual void accept(SiconosVisitor& tourist) const { tourist.visit(*this); } \
  virtual void acceptSerializerBase(SiconosVisitor& serializer) { serializer.visit(*static_cast<const BASE *>(this)); } \
  virtual void acceptSerializer(SiconosVisitor& serializer) { serializer.visit(*this); } \
  virtual inline Type::Siconos acceptType(FindType& ft) const { return ft.visit(*static_cast<const BASE *>(this)); } \

#define ACCEPT_BASE_NONVIRTUAL_VISITORS(BASE)                           \
  template<typename Archive> friend class SiconosSerializer;            \
  void acceptBase(SiconosVisitor& tourist) const { tourist.visit(*static_cast<const BASE *>(this)); } \
  void accept(SiconosVisitor& tourist) const { tourist.visit(*this); } \
  void acceptSerializerBase(SiconosVisitor& serializer) { serializer.visit(*static_cast<const BASE *>(this)); } \
  void acceptSerializer(SiconosVisitor& serializer) { serializer.visit(*this); } \
  inline Type::Siconos acceptType(FindType& ft) const { return ft.visit(*static_cast<const BASE *>(this)); } \

#define ACCEPT_BASE_SP_VISITORS(BASE)                                   \
  virtual void acceptSPBase(SP::SiconosVisitor tourist) { tourist->visit(std11::static_pointer_cast<BASE>(shared_from_this())); }\
  virtual void acceptSP(SP::SiconosVisitor tourist) { tourist->visit(shared_from_this()); }

#define ACCEPT_BASE_VISITORS(BASE)                           \
  ACCEPT_BASE_SP_VISITORS(BASE)                              \
  ACCEPT_BASE_STD_VISITORS(BASE)                             \

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
  Type::Siconos visit(const X&) const { return Type::X; };  \

#undef REGISTER_STRUCT
#define REGISTER_STRUCT(X) REGISTER(X)

#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y)                                         \
  Type::Siconos visit(const X&) const { return Type::Y; }; \

#define REGISTER_BASE_EXTERN(X,Y) REGISTER_BASE(X,Y)

struct FindType
{
  SICONOS_VISITABLES()
};

#define SICONOS_VISITOR_FAIL(X)                                         \
  RuntimeException::selfThrow                                           \
      ( SICONOS_VISITOR_QUOTE(you must define a visit function for X in a derived class of SiconosVisitor))

/* the base visitor */
#undef REGISTER
#define REGISTER(X)             \
  virtual void visit(std11::shared_ptr<X>) { SICONOS_VISITOR_FAIL(SP :: X);}; \
  virtual void visit(X&) { SICONOS_VISITOR_FAIL(X);};                   \
  virtual void visit(const X&) { SICONOS_VISITOR_FAIL(X);};             \

#undef REGISTER_STRUCT
#define REGISTER_STRUCT(X) REGISTER(X)

#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y) REGISTER(X)

#define REGISTER_BASE_EXTERN(X,Y) REGISTER_BASE(X,Y)

struct SiconosVisitor
{
  typedef struct{} arguments_type;
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
}


TYPEDEF_SPTR(SiconosVisitor)

#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN


#endif /* SiconosVisitor_hpp */
