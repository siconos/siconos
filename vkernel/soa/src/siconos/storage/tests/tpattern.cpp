
#include <siconos/storage/pattern/pattern.hpp>

using namespace siconos::storage::pattern;
using namespace siconos::storage::some;

static_assert(std::is_same_v<decltype(cons(int{}, gather<char, float>{})),
                             gather<int, char, float>>);

static_assert(std::is_same_v<decltype(append(gather<int>{}, gather<float>{})),
              gather<int, gather<float>>>);

static_assert(std::is_same_v<decltype(concat(gather<int>{}, gather<float>{})),
                             gather<int, float>>);


// static_assert(
//     attribute_name(siconos::storage::pattern::attribute<
//                    "hello", siconos::storage::some::scalar>{}) ==
//     "hello");

int main(){};
