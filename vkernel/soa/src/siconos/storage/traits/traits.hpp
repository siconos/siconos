#pragma once

#include <concepts>
#include <tuple>

#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"

namespace siconos::storage::traits {

using siconos::storage::pattern::nth_t;
using siconos::storage::pattern::param;
using siconos::storage::pattern::rec;
namespace match = siconos::storage::pattern::match;

static auto translate = rec([]<typename E, typename T>(auto&& translate, E,
                                                       T) {
  if constexpr (match::type_t<T>) {
    // return the embedded type
    return typename T::xtype{};
  }
  else if constexpr (!match::attribute<T>) {
    // return the type itself
    return T{};
  }
  else if constexpr (std::derived_from<T, some::undefined_polymorphic_type>) {
    return mp::unpack((typename T::types{}), []<typename... Ts>(Ts...) {
      return
          typename E::template variant<decltype(translate(E{}, Ts{}))...>{};
    });
  }
  else if constexpr (std::derived_from<T, some::undefined_tuple>) {
    return mp::unpack((typename T::types{}), []<typename... Ts>(Ts...) {
      return typename E::template tuple<decltype(translate(E{}, Ts{}))...>{};
    });
  }
  else if constexpr (std::derived_from<T, some::boolean>) {
    return typename E::boolean{};
  }
  else if constexpr (std::derived_from<T, some::scalar>) {
    return typename E::scalar{};
  }
  else if constexpr (std::derived_from<T, some::indice>) {
    return typename E::indice{};
  }
  else if constexpr (std::derived_from<T, some::integer>) {
    return typename E::integer{};
  }
  else if constexpr (std::derived_from<T, some::undefined_indice_parameter>) {
    return mp::get_m<param<T::name>>(typename E::params{});
  }
  else if constexpr (std::derived_from<T, some::undefined_indice_value>) {
    return T{};
  }
  else if constexpr (std::derived_from<T, some::undefined_type_parameter>) {
    return mp::get_m<param<T::name>>(typename E::params{});
  }

  else if constexpr (std::derived_from<
                         T, some::undefined_unbounded_collection>) {
    return typename E::template unbounded_collection<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T,
                                       some::undefined_bounded_collection>) {
    return typename E::template bounded_collection<
        decltype(translate(E{}, typename T::type{})),
        decltype(translate(E{}, nth_t<0, typename T::sizes>{}))::value>{};
  }

  else if constexpr (std::derived_from<T,
                                       some::structure_for_assembled_data>) {
    if constexpr (std::derived_from<T, some::undefined_diagonal_matrix> &&
                  std::derived_from<T, some::unbounded_storage>) {
      return
          typename E::template assembled_diagonal_matrix<decltype(translate(
              E{}, typename T::type{}))>{};
    }
    else if constexpr (std::derived_from<T,
                                         some::undefined_unbounded_matrix>) {
      return typename E::template assembled_matrix<decltype(translate(
          E{}, typename T::type{}))>{};
    }
    else if constexpr (std::derived_from<T,
                                         some::undefined_unbounded_vector>) {
      return typename E::template assembled_vector<decltype(translate(
          E{}, typename T::type{}))>{};
    }
  }
  else if constexpr (std::derived_from<T, some::undefined_vector>) {
    return typename E::template vector<
        decltype(translate(E{}, typename T::type{})),
        decltype(translate(E{}, nth_t<0, typename T::sizes>{}))::value>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_array>) {
    return typename E::template array<
        decltype(translate(E{}, typename T::type{})),
        decltype(translate(E{}, nth_t<0, typename T::sizes>{}))::value>{};
  }

  else if constexpr (std::derived_from<T, some::undefined_diagonal_matrix> &&
                     std::derived_from<T, some::unbounded_storage>) {
    return typename E::template unbounded_diagonal_matrix<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_diagonal_matrix>) {
    return typename E::template diagonal_matrix<
        decltype(translate(E{}, typename T::type{})),
        decltype(translate(E{}, nth_t<0, typename T::sizes>{}))::value>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_matrix>) {
    using nrows = decltype(translate(E{}, nth_t<0, typename T::sizes>{}));
    using ncols = decltype(translate(E{}, nth_t<1, typename T::sizes>{}));
    return typename E::template matrix<decltype(translate(
                                           E{}, typename T::type{})),
                                       nrows::value, ncols::value>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_unbounded_matrix>) {
    return typename E::template unbounded_matrix<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_unbounded_vector>) {
    return typename E::template unbounded_vector<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_map>) {
    return typename E::template map<
        decltype(translate(E{}, (typename T::types{})[0_c])),
        decltype(translate(E{}, (typename T::types{})[1_c]))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_graph>) {
    return typename E::template graph<
        decltype(translate(E{}, (typename T::types{}[0_c]))),
        decltype(translate(E{}, (typename T::types{}[1_c])))>{};
  }
  else if constexpr (std::derived_from<T, some::undefined_vdescriptor>) {
    return typename E::template vdescriptor<decltype(translate(
        E{}, typename T::type{}))>{};
  }
  else if constexpr (std::derived_from<T, some::item_ref<typename T::type>>) {
    return typename E::template item_ref<typename T::type>{};
  }
  else {
    []<typename Attr = T, bool flag = false>() {
      static_assert(flag, "cannot translate attribute");
    }();
  }
});

template <typename T, typename E>
concept translatable = requires(T t) { translate(E{}, t); };

template <typename E>
struct config {
  template <translatable<E> A>
  struct convert {
    using type = decltype(translate(E{}, A{}));
  };
};

template <typename T>
using value_type = typename T::value_type;  // eigen ok.

}  // namespace siconos::storage::traits
