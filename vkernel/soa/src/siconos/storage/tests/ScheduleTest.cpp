#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <cstdio>
#include <string>
#include <vector>

#include "siconos/simul/moreau_jean_schedule.hpp"
#include "siconos/storage/schedule.hpp"

namespace {

using namespace siconos::storage::schedule;

// a toy computation: its stages write their names
struct toy_context {
  std::vector<std::string> log;
  bool more = false;
};

struct first : stage<"first"> {
  static void run(toy_context& c) {
    c.log.emplace_back("first");
    c.more = true;
  }
};

struct second : stage<"second"> {
  static void run(toy_context& c) { c.log.emplace_back("second"); }
};

struct never : stage<"never"> {
  static void run(toy_context& c) { c.log.emplace_back("never"); }
};

struct more : condition<"more"> {
  static bool holds(const toy_context& c) { return c.more; }
};

struct less : condition<"less"> {
  static bool holds(const toy_context& c) { return !c.more; }
};

using toy = sequence<"toy", first, when<more, second>, when<less, never>>;

class ScheduleTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ScheduleTest);
  CPPUNIT_TEST(testRunInOrder);
  CPPUNIT_TEST(testObserver);
  CPPUNIT_TEST(testDescribe);
  CPPUNIT_TEST(testMoreauJean);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testRunInOrder() {
    toy_context c;
    run<toy>(c);
    CPPUNIT_ASSERT((c.log == std::vector<std::string>{"first", "second"}));
  }

  void testObserver() {
    toy_context c;
    std::vector<std::string> seen;
    run<toy>(c, [&seen](std::string_view name, auto&& action) {
      seen.emplace_back(name);
      action();
    });
    CPPUNIT_ASSERT((seen == std::vector<std::string>{"first", "second"}));
    CPPUNIT_ASSERT((c.log == seen));
  }

  void testDescribe() {
    CPPUNIT_ASSERT_EQUAL(std::string("toy:\n"
                                     "  first\n"
                                     "  if more:\n"
                                     "    second\n"
                                     "  if less:\n"
                                     "    never\n"),
                         describe<toy>());
  }

  // the Moreau-Jean time step, as the simulations run it: a change of the
  // algorithm changes this text
  void testMoreauJean() {
    const auto text = describe<siconos::simul::moreau_jean::time_step_schedule>();
    std::fputs(text.c_str(), stdout);
    CPPUNIT_ASSERT_EQUAL(
        std::string("Moreau-Jean time step:\n"
                    "  for each element:\n"
                    "    the relation sets hold every interaction\n"
                    "    M^-1 f, if f or K change\n"
                    "    the jacobians H, if the relations change\n"
                    "    the free state: v, q at step + 1 without the contact impulses\n"
                    "    the outputs y, ydot at step\n"
                    "    the free outputs y, ydot at step + 1\n"
                    "  the active set: the interactions whose gap may close, the systems they "
                    "involve\n"
                    "  if there are contacts:\n"
                    "    the sizes and offsets of the assembled problem\n"
                    "    for each element:\n"
                    "      the jacobians H, if the relations do not change\n"
                    "    H of the active interactions\n"
                    "    for each element:\n"
                    "      M of the involved systems\n"
                    "      K of the involved FEM systems\n"
                    "      lambda, ydot, mu of the active interactions\n"
                    "      the Newton impact law on the free ydot at step + 1\n"
                    "    q = H v_free, the free output derivatives\n"
                    "    W = H M^-1 H^T, or H (M + h^2 theta^2 K)^-1 H^T\n"
                    "    the frictional contact problem (W, q, mu): lambda\n"
                    "    the velocity correction M^-1 H^T lambda\n"
                    "    for each element:\n"
                    "      lambda of the active interactions from the solution\n"
                    "      v, q at step + 1 of the involved systems\n"),
        text);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ScheduleTest);

}  // namespace
