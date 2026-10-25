/* One rigid body type carrying disks and polytopes through per-kind sparse
   sets (collision/shape_of.hpp), on the registry backend.

   Checks: an item_ref to a shape can live in a dynamic property; has_shape,
   shape_of and set_shape answer per body; a body has exactly one shape kind;
   bodies_with iterates one kind only; storage::remove keeps the sets
   consistent (the last body takes the removed index with its shape). The
   static "shape" declaration gives the same answers at compile time.

   The per-kind answers are asserted both at compile time (file/block scope
   static_assert, verified when this TU builds) and at run time
   (CPPUNIT_ASSERT_MESSAGE below), so a failure shows up in ctest. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include <cmath>
#include <stdexcept>
#include <string>

#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/shape_of.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
using storage::pattern::symbol;

using disk_shape = collision::shape::disk;
using polytope_shape = collision::shape::polytope;

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

// per-kind world: one rigid body type, the shape kind per body
struct body : model::lagrangian_ds {};
struct data_t
    : storage::make_runtime_database<
          env, body, disk_shape, polytope_shape,
          storage::with_properties<
              storage::wrapped<body, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<polytope_shape, some::unbounded_collection>,
              storage::dynamic_storage<
                  storage::attached<body, symbol<"disk_shape">, some::item_ref<disk_shape>>>,
              storage::dynamic_storage<storage::attached<body, symbol<"polytope_shape">,
                                                         some::item_ref<polytope_shape>>>>> {};

// static world: the shape kind is part of the ds type
struct disk_ds : model::lagrangian_ds {};
struct sdata_t
    : storage::make_runtime_database<
          env, disk_ds, disk_shape,
          storage::with_properties<
              storage::wrapped<disk_ds, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::attached<disk_ds, symbol<"shape">, some::item_ref<disk_shape>>>> {};
}  // namespace siconos::config

namespace {
template <typename F>
bool throws(F&& f) {
  try {
    f();
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}
}  // namespace

class RigidShapesTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(RigidShapesTest);
  CPPUNIT_TEST(testRigidShapes);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp() override {}
  void tearDown() override {}

  void testRigidShapes() {
    using namespace siconos;
    using collision::shape::disk;
    using collision::shape::polytope;

    {
      config::data_t data;
      using body_h = decltype(storage::add<config::body>(data));
      static_assert(collision::per_kind_shapes_v<body_h>);
      static_assert(collision::may_have_shape<disk, body_h>());
      static_assert(collision::may_have_shape<polytope, body_h>());
      static_assert(!collision::always_shape<disk, body_h>());

      auto dsh = storage::add<config::disk_shape>(data);
      dsh.radius() = 0.5;
      auto psh = storage::add<config::polytope_shape>(data);
      psh.nvertices() = 4;
      psh.vertices()[0] = {1, -1};
      psh.vertices()[1] = {1, 1};
      psh.vertices()[2] = {-1, 1};
      psh.vertices()[3] = {-1, -1};

      auto b0 = storage::add<config::body>(data);
      auto b1 = storage::add<config::body>(data);
      auto b2 = storage::add<config::body>(data);
      b0.q() = {0., 0., 0.};
      b1.q() = {1., 0., 0.};
      b2.q() = {2., 0., 0.};

      collision::set_shape(b0, dsh);
      collision::set_shape(b1, psh);
      collision::set_shape(b2, psh);

      CPPUNIT_ASSERT_MESSAGE("body 0 is a disk only", collision::has_shape<disk>(b0) &&
                                                          !collision::has_shape<polytope>(b0));
      CPPUNIT_ASSERT_MESSAGE("body 1 is a polytope only", collision::has_shape<polytope>(b1) &&
                                                              !collision::has_shape<disk>(b1));
      CPPUNIT_ASSERT_MESSAGE("shape_of<disk> reads the radius",
                             collision::shape_of<disk>(b0).radius() == 0.5);
      CPPUNIT_ASSERT_MESSAGE("shape_of<polytope> reads the vertices",
                             collision::shape_of<polytope>(b1).nvertices() == 4);
      CPPUNIT_ASSERT_MESSAGE("shape_of of a missing kind throws",
                             throws([&] { collision::shape_of<disk>(b1); }));
      CPPUNIT_ASSERT_MESSAGE("a second shape kind throws",
                             throws([&] { collision::set_shape(b1, dsh); }));
      CPPUNIT_ASSERT_MESSAGE("the refused kind left no entry",
                             !collision::has_shape<disk>(b1));

      auto& polys = collision::bodies_with<config::body, polytope>(data);
      auto& disks = collision::bodies_with<config::body, disk>(data);
      CPPUNIT_ASSERT_MESSAGE("two polytope bodies, one disk body",
                             polys.size() == 2 && disks.size() == 1);

      // body 2 (polytope) takes index 0 when body 0 (disk) is removed
      storage::remove(data, b0);
      auto& qs = storage::attr_values<config::body, "q">(data, 0);
      CPPUNIT_ASSERT_MESSAGE("body 2 moved to index 0", std::size(qs) == 2 && qs[0](0) == 2.);
      CPPUNIT_ASSERT_MESSAGE("the disk entry went with body 0",
                             polys.size() == 2 && disks.size() == 0);
      CPPUNIT_ASSERT_MESSAGE("the polytope entry of body 2 followed it to index 0",
                             polys.contains(0) && polys.contains(1) && !polys.contains(2));
      auto nb0 = storage::make_handle(data, storage::index<config::body, std::size_t>(0));
      CPPUNIT_ASSERT_MESSAGE("the moved body still reads its polytope",
                             collision::has_shape<polytope>(nb0) &&
                                 collision::shape_of<polytope>(nb0).nvertices() == 4);
    }

    {
      config::sdata_t data;
      using ds_h = decltype(storage::add<config::disk_ds>(data));
      static_assert(!collision::per_kind_shapes_v<ds_h>);
      static_assert(collision::always_shape<disk, ds_h>());
      static_assert(!collision::may_have_shape<polytope, ds_h>());

      auto dsh = storage::add<config::disk_shape>(data);
      dsh.radius() = 0.25;
      auto d = storage::add<config::disk_ds>(data);
      collision::set_shape(d, dsh);
      CPPUNIT_ASSERT_MESSAGE(
          "static declaration: a disk ds is a disk",
          collision::has_shape<disk>(d) && !collision::has_shape<polytope>(d));
      CPPUNIT_ASSERT_MESSAGE("static declaration: shape_of<disk>",
                             collision::shape_of<disk>(d).radius() == 0.25);
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(RigidShapesTest);

int main() {
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(RigidShapesTest::suite());
  return runner.run() ? 0 : 1;
}
