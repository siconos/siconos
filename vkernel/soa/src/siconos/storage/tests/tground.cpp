#include <siconos/storage/mp/mp.hpp>

using namespace siconos::storage;
using namespace siconos::storage::mp;

struct A {
  using a_t = void;
};
struct B : A {};
struct C : B {};
struct D {};
template <typename T>
struct E : T {};

template <typename T>
concept match_a_t = requires { typename T::a_t; };

static_assert(std::is_same_v<decltype(filter(std::tuple<A, B, C, D, E<B>>{},
                                             derive_from<A>)),
                             std::tuple<A, B, C, E<B>>>);

static_assert(match_a_t<A>);
static_assert(match_a_t<B>);

static_assert(
    std::is_same_v<decltype(filter(std::tuple<A, B, C, D, E<B>>{},
                                   is_a_model<[]<typename T>() consteval {
                                     return match_a_t<T>;
                                   }>)),
                   std::tuple<A, B, C, E<B>>>);

static_assert(contains(std::tuple<int, char, double>{}, int{}));

static_assert(!any_of(std::tuple<float, float, double>{2.0, 3.0, 2.0},
                      []<typename T>(T) { return std::is_integral<T>(); }));

static_assert(any_of(std::tuple<float, char, double>{2.0, 'a', 2.0},
                     []<typename T>(T) { return std::is_integral<T>(); }));

static_assert(
    std::is_same_v<decltype(filter(std::make_tuple(4, 5, 6.0, 'a', "hello"),
                                   []<typename T>(T) {
                                     return std::is_floating_point<T>();
                                   })),
                   std::tuple<double>>);

static_assert(transform(make_tuple(1, 2, 3), [](auto x) { return x + 1; }) ==
              make_tuple(2, 3, 4));

static_assert(filter(all_type_c(make_tuple(A{}, B{}, C{}, D{})),
                     mp::derive_from<A>) ==
              all_type_c(make_tuple(A{}, B{}, C{})));

static_assert((make_tuple(make_tuple(type_c<A>, type_c<A>, type_c<B>,
                                     type_c<C>)) |
               unique) == make_tuple(type_c<A>, type_c<B>, type_c<C>));
static_assert(to_set(all_type_c(tuple_unique(make_tuple(A{}, B{}, A{},
                                                        C{})))) ==
              to_set(all_type_c(make_tuple(A{}, B{}, C{}))));
int main() {}
