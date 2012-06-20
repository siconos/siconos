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
// === Test file for LinearDS class ===
// -> call successively all the constructors and display new object.

#include "Model.hpp"
#include "LagrangianDS.hpp"
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"
#include <sys/time.h>
#include <iostream>
#include <math.h>

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // problem size
    unsigned int size = 3;
    // initial state
    SP::SiconosVector q0(new SiconosVector(3));
    (*q0)(0) = 1;
    (*q0)(1) = 2;
    (*q0)(2) = 3;
    SP::SiconosVector v0(new SiconosVector(3));
    (*v0)(0) = 4;
    (*v0)(1) = 5;
    (*v0)(2) = 6;

    // constructor from min set of data
    cout << "======== Test 1 ============= " << endl;
    SP::LagrangianDS lds1(new LagrangianDS(1, size, *q0, *v0));
    lds1->display();

    cout << "======== Test 2 ============= " << endl;
    // constructor from a set of data, M from a plugin.
    lds1.reset(new LagrangianDS(1, size, *q0, *v0, "LagPlugin:computeMass"));
    lds1->computeMass(2);
    lds1->mass()->display();
    lds1->computeMass(2, v0);
    lds1->mass()->display();

    cout << "======== Test 3 ============= " << endl;
    // constructor from a set of data, Mass a given matrix.
    SP::SiconosMatrix A(new SiconosMatrix("mat.dat", true));
    lds1.reset(new LagrangianDS(1, size, *q0, *v0, *A));
    lds1->mass()->display();

    cout << "======== Test 4 ============= " << endl;
    // getter/setter
    // by "value" :
    lds1.reset(new LagrangianDS(1, size, *q0, *v0));
    lds1->setMass(*A);
    lds1->getMass().display();
    lds1->setFInt(*q0);
    lds1->getFInt().display();
    lds1->setFExt(*v0);
    lds1->getFExt().display();
    lds1->setQNLInertia(*q0);
    lds1->getQNLInertia().display();
    lds1->setJacobianFIntq(*A);
    lds1->getJacobianFIntq().display();
    lds1->setJacobianFintVelocity(*A);
    lds1->getJacobianFintVelocity().display();
    lds1->setJacobianQNLInertiaq(*A);
    lds1->getJacobianQNLInertiaq().display();
    lds1->setJacobianQNLInertiaVelocity(*A);
    lds1->getJacobianQNLInertiaVelocity().display();

    cout << "======== Test 5 ============= " << endl;
    // getter/setter
    // by pointer :
    SP::SiconosMatrix B(new SiconosMatrix("matB.dat", true));
    lds1->setMassPtr(B);
    lds1->mass()->display();
    lds1->setFIntPtr(v0);
    lds1->fInt()->display();
    lds1->setFExtPtr(q0);
    lds1->fExt()->display();
    lds1->setQNLInertiaPtr(v0);
    lds1->qNLInertia()->display();
    lds1->setJacobianFIntqPtr(B);
    lds1->jacobianFIntq()->display();
    lds1->setJacobianFintVelocityPtr(B);
    lds1->jacobianFIntVelocity()->display();
    lds1->setJacobianQNLInertiaqPtr(B);
    lds1->jacobianQNLInertiaq()->display();
    lds1->setJacobianQNLInertiaVelocityPtr(B);
    lds1->jacobianQNLInertiaVelocity()->display();

    B->zero();
    lds1->mass()->display();
    lds1->fInt()->display();
    lds1->fExt()->display();
    lds1->qNLInertia()->display();
    lds1->jacobianFIntq()->display();
    lds1->jacobianFIntVelocity()->display();
    lds1->jacobianQNLInertiaq()->display();
    lds1->jacobianQNLInertiaVelocity()->display();

    cout << "======== Test 6 ============= " << endl;
    // Plugin
    lds1.reset(new LagrangianDS(1, size, *q0, *v0));
    lds1->setComputeMassFunction("LagPlugin.so", "computeMass");
    lds1->computeMass(2);
    lds1->mass()->display();
    lds1->setComputeFIntFunction("LagPlugin.so", "computeFInt");
    lds1->computeFInt(2);
    lds1->fInt()->display();
    lds1->computeFInt(3, v0, q0);
    lds1->fInt()->display();
    lds1->setComputeFExtFunction("LagPlugin.so", "computeFExt");
    lds1->computeFExt(2);
    lds1->fExt()->display();
    lds1->setComputeQNLInertiaFunction("LagPlugin.so", "computeQNLInertia");
    lds1->computeQNLInertia();
    lds1->qNLInertia()->display();
    lds1->computeQNLInertia(v0, q0);
    lds1->qNLInertia()->display();
    lds1->setComputeJacobianFIntqFunction("LagPlugin.so", "computeJacobianFIntq");
    lds1->computeJacobianFIntq(2);
    lds1->jacobianFIntq()->display();
    lds1->computeJacobianFIntq(2, v0, q0);
    lds1->jacobianFIntq()->display();
    lds1->setComputeJacobianFintVelocityFunction("LagPlugin.so", "computeJacobianFintVelocity");
    lds1->computeJacobianFintVelocity(2);
    lds1->jacobianFIntVelocity()->display();
    lds1->computeJacobianFintVelocity(2, v0, q0);
    lds1->jacobianFIntVelocity()->display();
    lds1->setComputeJacobianQNLInertiaqFunction("LagPlugin.so", "computeJacobianQNLInertiaq");
    lds1->computeJacobianQNLInertiaq();
    lds1->jacobianQNLInertiaq()->display();
    lds1->computeJacobianQNLInertiaq(v0, q0);
    lds1->jacobianQNLInertiaq()->display();
    lds1->setComputeJacobianQNLInertiaVelocityFunction("LagPlugin.so", "computeJacobianQNLInertiaVelocity");
    lds1->computeJacobianQNLInertiaVelocity();
    lds1->jacobianQNLInertiaVelocity()->display();
    lds1->computeJacobianQNLInertiaVelocity(v0, q0);
    lds1->jacobianQNLInertiaVelocity()->display();

    // set f, u, B
    SP::SiconosVector u(new SiconosVector(2));
    (*u)(0) = 1.8;
    (*u)(1) = 1.4;
    SP::SiconosVector f(new SiconosVector(size));
    (*f)(0) = 2.5;
    (*f)(1) = 4;
    (*f)(2) = 9;


    SP::SiconosVector u2(new SiconosVector(1));
    (*u2)(0) = 34;

    cout << "======== Test 6 ============= " << endl;
    // xml constructor
    // parse xml file:
    xmlDocPtr doc;
    xmlNodePtr cur;
    doc = xmlParseFile("LagrangianDS_test.xml");
    if (!doc)
      XMLException::selfThrow("Document not parsed successfully");
    cur = xmlDocGetRootElement(doc);
    if (!cur)
    {
      XMLException::selfThrow("empty document");
      xmlFreeDoc(doc);
    }

    // get rootNode
    xmlNode * node;

    if (xmlStrcmp(cur->name, (const xmlChar *) "SiconosModel"))
    {
      XMLException::selfThrow("document of the wrong type, root node !=SiconosModel");
      xmlFreeDoc(doc);
    }

    // look for LagrangianDS node
    node = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
    cout << node->name << endl;
    node = SiconosDOMTreeTools::findNodeChild(node, "DS_Definition");
    cout << node->name << endl;
    node = SiconosDOMTreeTools::findNodeChild(node, "LagrangianDS");
    cout << node->name << endl;

    // xml constructor
    SP::LagrangianDSXML tmpxml(new LagrangianDSXML(node, false));

    lds1.reset(new LagrangianDS(tmpxml));
    lds1->computeMass(2);
    lds1->mass()->display();
    lds1->computeFInt(2);
    lds1->fInt()->display();
    lds1->computeFInt(3, v0, q0);
    lds1->fInt()->display();
    lds1->computeFExt(2);
    lds1->fExt()->display();
    lds1->computeQNLInertia();
    lds1->qNLInertia()->display();
    lds1->computeQNLInertia(v0, q0);
    lds1->qNLInertia()->display();
    lds1->computeJacobianFIntq(2);
    lds1->jacobianFIntq()->display();
    lds1->computeJacobianFIntq(2, v0, q0);
    lds1->jacobianFIntq()->display();
    lds1->computeJacobianFintVelocity(2);
    lds1->jacobianFIntVelocity()->display();
    lds1->computeJacobianFintVelocity(2, v0, q0);
    lds1->jacobianFIntVelocity()->display();
    lds1->computeJacobianQNLInertiaq();
    lds1->jacobianQNLInertiaq()->display();
    lds1->computeJacobianQNLInertiaq(v0, q0);
    lds1->jacobianQNLInertiaq()->display();
    lds1->computeJacobianQNLInertiaVelocity();
    lds1->jacobianQNLInertiaVelocity()->display();
    lds1->computeJacobianQNLInertiaVelocity(v0, q0);
    lds1->jacobianQNLInertiaVelocity()->display();
  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/BouncingBall\'" << endl;
  }
}
