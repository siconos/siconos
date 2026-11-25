// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <cppunit/TestPath.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include <functional>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  // The object to run tests
  CppUnit::TextUi::TestRunner runner;

  // Get test classes that have been registered
  CppUnit::TestFactoryRegistry& registry = CppUnit::TestFactoryRegistry::getRegistry();

  // Put tests into the runner
  CppUnit::Test* all = registry.makeTest();
  runner.addTest(all);

  // the test named by the argument, if any: a suite or a test case
  // (Suite::testMethod), its full name or its end after a "::" (a suite in
  // an anonymous namespace is "(anonymous namespace)::Suite")
  std::string path;
  if (argc > 1) {
    const std::string name = argv[1];
    const std::string suffix = "::" + name;
    std::function<CppUnit::Test*(CppUnit::Test*)> find =
        [&](CppUnit::Test* test) -> CppUnit::Test* {
      const std::string& n = test->getName();
      if (n == name || (n.size() > suffix.size() &&
                        n.compare(n.size() - suffix.size(), suffix.size(), suffix) == 0))
        return test;
      for (int i = 0; i < test->getChildTestCount(); ++i)
        if (CppUnit::Test* found = find(test->getChildTestAt(i))) return found;
      return nullptr;
    };
    CppUnit::Test* test = find(all);
    if (!test) {
      std::cerr << "no test " << name << std::endl;
      return 1;
    }
    CppUnit::TestPath test_path;
    all->findTestPath(test, test_path);
    // relative to the suite of the runner, which holds all
    path = test_path.toString().substr(1);
  }

  // Run tests; a failure is also the exit code, for callers that do not
  // read the output
  return runner.run(path, false, true, false) ? 0 : 1;
}
