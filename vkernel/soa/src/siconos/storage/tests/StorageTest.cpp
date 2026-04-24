
#include "StorageTest.hpp"

#include <vector>

#include "siconos/storage/add.hpp"
#include "siconos/storage/default_interface.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/storage/make.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/traits/traits.hpp"

namespace siconos::config {

using namespace siconos::storage;
using namespace siconos::storage::pattern;

struct item_a : item {
  using attributes = gather<attribute<"value", some::scalar>>;
};

struct item_b : item {
  using attributes = gather<attribute<"ref", some::item_ref<item_a>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) ref() { return attr<"ref">(*self()); };

    decltype(auto) href()
    {
      return make_ref_handle(self()->data(), attr<"ref">(*self()));
    };
  };
};

template <typename Item>
struct env {

  using boolean = bool;
  using scalar = float;
  using indice = std::size_t;
  using integer = int;

  template <typename T>
  using unbounded_collection = std::vector<T>;

  template <typename T>
  using item_ref = index<T, indice>;

  template <typename T>
  using default_storage = std::array<T, 1>;
};

struct params {};

struct make
    : storage::make<
          env, item_a, item_b,
          with_properties<wrapped<item_a, some::unbounded_collection>,
                          wrapped<item_b, some::unbounded_collection>>> {};

}  // namespace siconos::config

namespace config = siconos::config;
namespace store = siconos::storage;
namespace ct = store::mp;

CPPUNIT_TEST_SUITE_REGISTRATION(StorageTest);

void StorageTest::setUp() {}

void StorageTest::tearDown() {}

// Tests on index and handle
void StorageTest::testStorage0()
{
  auto data = siconos::config::make();

  auto a1 = store::add<config::item_a>(data);
  auto a2 = store::add<config::item_a>(data);

  // With unbounded storage index is incremented at each addition.
  CPPUNIT_ASSERT(a2.index().value() == a1.index().value() + 1);

  auto b1 = store::add<config::item_b>(data);
  auto b2 = store::add<config::item_b>(data);

  b1.ref() = a1.index();
  b2.ref() = a2.index();

  // Check index as lvalue reference is set from handle.
  CPPUNIT_ASSERT(b1.ref().value() == a1.index().value());
  CPPUNIT_ASSERT(b2.ref().value() == a2.index().value());

  auto a3 = store::add<config::item_a>(data);
  auto a4 = store::add<config::item_a>(data);

  CPPUNIT_ASSERT(a3.index().value() == 2);
  CPPUNIT_ASSERT(a4.index().value() == 3);

  b1.href() = a3;
  b2.href() = a4;

  CPPUNIT_ASSERT(b1.ref().value() == a3.index().value());
  CPPUNIT_ASSERT(b2.ref().value() == a4.index().value());
};
