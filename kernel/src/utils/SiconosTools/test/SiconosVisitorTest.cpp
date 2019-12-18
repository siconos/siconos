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
#include "SiconosVisitorTest.hpp"
class ObjectA;
class ObjectB;

#include "../SiconosVisitables.hpp"
#undef SICONOS_VISITABLES
#define SICONOS_VISITABLES()                    \
  KERNEL_CLASSES()                              \
  REGISTER(ObjectA)                             \
  REGISTER(ObjectB)

#include "TypeName.hpp"
#include "../SiconosVisitor.hpp"



// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosVisitorTest);


void SiconosVisitorTest::setUp()
{
}

void SiconosVisitorTest::tearDown()
{
}

class DynamicalSystem
{
public:
  VIRTUAL_ACCEPT_VISITORS(DynamicalSystem);
  virtual ~DynamicalSystem() {}
};

class LagrangianDS : public DynamicalSystem
{
public:
  ACCEPT_STD_VISITORS();
};

/* TypeOf */
void SiconosVisitorTest::t1()
{
  DynamicalSystem *ds = new LagrangianDS();

  CPPUNIT_ASSERT(Type::value(*ds) == Type::LagrangianDS);

  delete(ds);

}

/* standard visitor */
void SiconosVisitorTest::t2()
{

  struct MyVisitor : public SiconosVisitor
  {
#if !defined(_MSC_VER)
    using SiconosVisitor::visit;
#endif

    void visit(const LagrangianDS& ds)
    {
      ;
    }

  };

  DynamicalSystem *ds = new LagrangianDS();

  try
  {
    MyVisitor myvisitor;

    ds->accept(myvisitor);

    delete(ds);
  }

  catch(...)
  {
    CPPUNIT_ASSERT(false);
    delete(ds);
  }


}

void SiconosVisitorTest::t3()
{

  DynamicalSystem *ds = new LagrangianDS();

  CPPUNIT_ASSERT(Type::name(*ds) == "LagrangianDS");

  delete ds;
}

struct Object
{

  VIRTUAL_ACCEPT_VISITORS();

  virtual ~Object() {};

};

struct ObjectA : public Object
{
  int id;
  int dummya;

  ACCEPT_STD_VISITORS();
  virtual ~ObjectA() {};

};

struct ObjectB : public Object
{
  int id;
  int dummyb;

  ACCEPT_STD_VISITORS();
  virtual ~ObjectB() {};

};

#define VISITOR_CLASSES()\
  REGISTER(ObjectA)\
  REGISTER(ObjectB)

#include "VisitorMaker.hpp"
using namespace Experimental;

struct GetId : public SiconosVisitor
{

  int result;

  template<typename T>
  void operator()(const T& obj)
  {
    result = obj.id;
  }
};


void SiconosVisitorTest::t4()
{
  Object *o;

  ObjectA ooa;
  ObjectB oob;

  ooa.id = 0;
  oob.id = 1;


  Visitor < Classes < ObjectA, ObjectB >,
          GetId >::Make visitor;

  o = & ooa;

  o->accept(visitor);

  CPPUNIT_ASSERT(visitor.result == 0);

  o = & oob;

  o->accept(visitor);

  CPPUNIT_ASSERT(visitor.result == 1);

};
