/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright INRIA.
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

#include "StorageTest.hpp"

#include <vector>

#include "siconos/vkernel/storage/add.hpp"
#include "siconos/vkernel/storage/default_interface.hpp"
#include "siconos/vkernel/storage/handle.hpp"
#include "siconos/vkernel/storage/make.hpp"
#include "siconos/vkernel/storage/memory.hpp"
#include "siconos/vkernel/storage/mp/mp.hpp"
#include "siconos/vkernel/storage/pattern/pattern.hpp"
#include "siconos/vkernel/storage/properties.hpp"
#include "siconos/vkernel/storage/registry.hpp"
#include "siconos/vkernel/storage/some/some.hpp"
#include "siconos/vkernel/storage/traits/traits.hpp"

namespace siconos::vkernel::config {

using namespace siconos::vkernel::storage;
using namespace siconos::vkernel::storage::pattern;

struct item_a : item {
  using attributes = gather<attribute<"value", some::scalar>>;
};

struct item_b : item {
  using attributes = gather<attribute<"ref", some::item_ref<item_a>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) ref() { return attr<"ref">(*self()); };

    decltype(auto) href() { return make_ref_handle(self()->data(), attr<"ref">(*self())); };
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

struct make : storage::make_runtime_database<
                  env, item_a, item_b,
                  with_properties<wrapped<item_a, some::unbounded_collection>,
                                  wrapped<item_b, some::unbounded_collection>>> {};

}  // namespace siconos::vkernel::config

namespace config = siconos::vkernel::config;
namespace store = siconos::vkernel::storage;
namespace ct = store::mp;

CPPUNIT_TEST_SUITE_REGISTRATION(StorageTest);

void StorageTest::setUp() {}

void StorageTest::tearDown() {}

// Tests on index and handle
void StorageTest::testStorage0() {
  auto data = siconos::vkernel::config::make();

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
