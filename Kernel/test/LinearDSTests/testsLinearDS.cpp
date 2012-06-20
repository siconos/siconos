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
#include "FirstOrderLinearDS.hpp"
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
    SP::SiconosVector x0(new SiconosVector(3));
    (*x0)(0) = 1;
    (*x0)(1) = 2;
    (*x0)(2) = 3;

    // constructor from min set of data
    SP::FirstOrderLinearDS lds1(new FirstOrderLinearDS(1, size, *x0));
    lds1->display();

    // constructor from a set of data, A from a plugin.
    lds1.reset(new FirstOrderLinearDS(1, size, *x0, "LDSPlugin.so", "computeA"));
    lds1->computeA(2);
    lds1->A()->display();

    // constructor from a set of data, A a given matrix.
    SP::SiconosMatrix A(new SiconosMatrix("mat.dat", true));
    lds1.reset(new FirstOrderLinearDS(1, *x0, *A));
    lds1->display();

    // set b, u, E
    SP::SiconosMatrix E(new SiconosMatrix("matE.dat", true));
    SP::SiconosVector u(new SiconosVector(2));
    (*u)(0) = 1.8;
    (*u)(1) = 1.4;
    SP::SiconosVector b(new SiconosVector(size));
    (*b)(0) = 2.5;
    (*b)(1) = 4;
    (*b)(2) = 9;

    // set by values
    lds1->setB(*b);
    lds1->setU(*u);
    lds1->setE(*E);
    lds1->display();

    // set with pointers
    lds1->setBPtr(b);
    lds1->setUPtr(u);
    lds1->setEPtr(E);
    lds1->display();

    // change set different size for u
    SP::SiconosVector u2(new SiconosVector(1));
    (*u2)(0) = 34;
    // lds1->setU(*u2); error -> should change uSize before
    lds1->setUSize(1);
    lds1->setU(*u2);
    //    lds1->setE(*E)
    lds1->display();
    // The same in reverse order
    lds1.reset(new FirstOrderLinearDS(1, *x0, *A));
    // set with pointers
    lds1->setBPtr(b);
    lds1->setUPtr(u);
    lds1->setEPtr(E);
    lds1->display();
    // set by values
    lds1->setB(*b);
    lds1->setU(*u);
    lds1->setE(*E);
    lds1->display();

    // set with plugins
    lds1->setComputeBFunction("LDSPlugin.so", "computeb");
    lds1->setComputeUFunction("LDSPlugin.so", "computeU");
    lds1->setComputeEFunction("LDSPlugin.so", "computeE");

    lds1->computeb(3);
    lds1->computeU(2);
    lds1->computeE(2);

    lds1->display();
    // xml constructor

    // parse xml file:
    xmlDocPtr doc;
    xmlNodePtr cur;
    doc = xmlParseFile("FirstOrderLinearDS_test.xml");
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

    // look for FirstOrderLinearDS node
    node = SiconosDOMTreeTools::findNodeChild(cur, "NSDS");
    cout << node->name << endl;
    node = SiconosDOMTreeTools::findNodeChild(node, "DS_Definition");
    cout << node->name << endl;
    node = SiconosDOMTreeTools::findNodeChild(node, "FirstOrderLinearDS");
    cout << node->name << endl;

    // xml constructor
    SP::FirstOrderLinearDSXML tmpxml(new FirstOrderLinearDSXML(node, false));

    lds1.reset(new FirstOrderLinearDS(tmpxml));
    lds1->computeE(3);
    lds1->computeb(2);

    lds1->display();

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
