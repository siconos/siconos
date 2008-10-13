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
// === Test file for LinearDS class ===
// -> call successively all the constructors and display new object.

#include "Model.h"
#include "LagrangianDS.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
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
    SimpleVector * q0 = new SimpleVector(3);
    (*q0)(0) = 1;
    (*q0)(1) = 2;
    (*q0)(2) = 3;
    SimpleVector * v0 = new SimpleVector(3);
    (*v0)(0) = 4;
    (*v0)(1) = 5;
    (*v0)(2) = 6;

    // constructor from min set of data
    cout << "======== Test 1 ============= " << endl;
    LagrangianDS * lds1 = new LagrangianDS(1, size, *q0, *v0);
    lds1->display();
    delete lds1;

    cout << "======== Test 2 ============= " << endl;
    // constructor from a set of data, M from a plugin.
    lds1 = new LagrangianDS(1, size, *q0, *v0, "LagPlugin:computeMass");
    lds1->computeMass(2);
    lds1->getMassPtr()->display();
    lds1->computeMass(2, v0);
    lds1->getMassPtr()->display();
    delete lds1;

    cout << "======== Test 3 ============= " << endl;
    // constructor from a set of data, Mass a given matrix.
    SiconosMatrix * A = new SiconosMatrix("mat.dat", true);
    lds1 = new LagrangianDS(1, size, *q0, *v0, *A);
    lds1->getMassPtr()->display();
    delete lds1;

    cout << "======== Test 4 ============= " << endl;
    // getter/setter
    // by "value" :
    lds1 = new LagrangianDS(1, size, *q0, *v0);
    lds1->setMass(*A);
    lds1->getMass().display();
    lds1->setFInt(*q0);
    lds1->getFInt().display();
    lds1->setFExt(*v0);
    lds1->getFExt().display();
    lds1->setQNLInertia(*q0);
    lds1->getQNLInertia().display();
    lds1->setJacobianQFInt(*A);
    lds1->getJacobianQFInt().display();
    lds1->setJacobianVelocityFInt(*A);
    lds1->getJacobianVelocityFInt().display();
    lds1->setJacobianQQNLInertia(*A);
    lds1->getJacobianQQNLInertia().display();
    lds1->setJacobianVelocityQNLInertia(*A);
    lds1->getJacobianVelocityQNLInertia().display();

    cout << "======== Test 5 ============= " << endl;
    // getter/setter
    // by pointer :
    SiconosMatrix * B = new SiconosMatrix("matB.dat", true);
    lds1->setMassPtr(B);
    lds1->getMassPtr()->display();
    lds1->setFIntPtr(v0);
    lds1->getFIntPtr()->display();
    lds1->setFExtPtr(q0);
    lds1->getFExtPtr()->display();
    lds1->setQNLInertiaPtr(v0);
    lds1->getQNLInertiaPtr()->display();
    lds1->setJacobianQFIntPtr(B);
    lds1->getJacobianQFIntPtr()->display();
    lds1->setJacobianVelocityFIntPtr(B);
    lds1->getJacobianVelocityFIntPtr()->display();
    lds1->setJacobianQQNLInertiaPtr(B);
    lds1->getJacobianQQNLInertiaPtr()->display();
    lds1->setJacobianVelocityQNLInertiaPtr(B);
    lds1->getJacobianVelocityQNLInertiaPtr()->display();

    B->zero();
    lds1->getMassPtr()->display();
    lds1->getFIntPtr()->display();
    lds1->getFExtPtr()->display();
    lds1->getQNLInertiaPtr()->display();
    lds1->getJacobianQFIntPtr()->display();
    lds1->getJacobianVelocityFIntPtr()->display();
    lds1->getJacobianQQNLInertiaPtr()->display();
    lds1->getJacobianVelocityQNLInertiaPtr()->display();

    cout << "======== Test 6 ============= " << endl;
    // Plugin
    delete lds1;
    lds1 = new LagrangianDS(1, size, *q0, *v0);
    lds1->setComputeMassFunction("LagPlugin.so", "computeMass");
    lds1->computeMass(2);
    lds1->getMassPtr()->display();
    lds1->setComputeFIntFunction("LagPlugin.so", "computeFInt");
    lds1->computeFInt(2);
    lds1->getFIntPtr()->display();
    lds1->computeFInt(3, v0, q0);
    lds1->getFIntPtr()->display();
    lds1->setComputeFExtFunction("LagPlugin.so", "computeFExt");
    lds1->computeFExt(2);
    lds1->getFExtPtr()->display();
    lds1->setComputeQNLInertiaFunction("LagPlugin.so", "computeQNLInertia");
    lds1->computeQNLInertia();
    lds1->getQNLInertiaPtr()->display();
    lds1->computeQNLInertia(v0, q0);
    lds1->getQNLInertiaPtr()->display();
    lds1->setComputeJacobianQFIntFunction("LagPlugin.so", "computeJacobianQFInt");
    lds1->computeJacobianQFInt(2);
    lds1->getJacobianQFIntPtr()->display();
    lds1->computeJacobianQFInt(2, v0, q0);
    lds1->getJacobianQFIntPtr()->display();
    lds1->setComputeJacobianVelocityFIntFunction("LagPlugin.so", "computeJacobianVelocityFInt");
    lds1->computeJacobianVelocityFInt(2);
    lds1->getJacobianVelocityFIntPtr()->display();
    lds1->computeJacobianVelocityFInt(2, v0, q0);
    lds1->getJacobianVelocityFIntPtr()->display();
    lds1->setComputeJacobianQQNLInertiaFunction("LagPlugin.so", "computeJacobianQQNLInertia");
    lds1->computeJacobianQQNLInertia();
    lds1->getJacobianQQNLInertiaPtr()->display();
    lds1->computeJacobianQQNLInertia(v0, q0);
    lds1->getJacobianQQNLInertiaPtr()->display();
    lds1->setComputeJacobianVelocityQNLInertiaFunction("LagPlugin.so", "computeJacobianVelocityQNLInertia");
    lds1->computeJacobianVelocityQNLInertia();
    lds1->getJacobianVelocityQNLInertiaPtr()->display();
    lds1->computeJacobianVelocityQNLInertia(v0, q0);
    lds1->getJacobianVelocityQNLInertiaPtr()->display();

    delete lds1;

    // set f, u, B
    SimpleVector *u = new SimpleVector(2);
    (*u)(0) = 1.8;
    (*u)(1) = 1.4;
    SimpleVector *f = new SimpleVector(size);
    (*f)(0) = 2.5;
    (*f)(1) = 4;
    (*f)(2) = 9;


    SimpleVector *u2 = new SimpleVector(1);
    (*u2)(0) = 34;
    delete u2;


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
    LagrangianDSXML * tmpxml = new LagrangianDSXML(node, false);

    lds1 = new LagrangianDS(tmpxml);
    lds1->computeMass(2);
    lds1->getMassPtr()->display();
    lds1->computeFInt(2);
    lds1->getFIntPtr()->display();
    lds1->computeFInt(3, v0, q0);
    lds1->getFIntPtr()->display();
    lds1->computeFExt(2);
    lds1->getFExtPtr()->display();
    lds1->computeQNLInertia();
    lds1->getQNLInertiaPtr()->display();
    lds1->computeQNLInertia(v0, q0);
    lds1->getQNLInertiaPtr()->display();
    lds1->computeJacobianQFInt(2);
    lds1->getJacobianQFIntPtr()->display();
    lds1->computeJacobianQFInt(2, v0, q0);
    lds1->getJacobianQFIntPtr()->display();
    lds1->computeJacobianVelocityFInt(2);
    lds1->getJacobianVelocityFIntPtr()->display();
    lds1->computeJacobianVelocityFInt(2, v0, q0);
    lds1->getJacobianVelocityFIntPtr()->display();
    lds1->computeJacobianQQNLInertia();
    lds1->getJacobianQQNLInertiaPtr()->display();
    lds1->computeJacobianQQNLInertia(v0, q0);
    lds1->getJacobianQQNLInertiaPtr()->display();
    lds1->computeJacobianVelocityQNLInertia();
    lds1->getJacobianVelocityQNLInertiaPtr()->display();
    lds1->computeJacobianVelocityQNLInertia(v0, q0);
    lds1->getJacobianVelocityQNLInertiaPtr()->display();

    delete tmpxml;

    delete f;
    delete u;
    delete B;
    delete A;
    delete lds1;
    delete q0;
    delete v0;

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
