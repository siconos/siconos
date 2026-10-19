/* CXD is a program dedicated to convex collision detection
 * (currently part of the Siconos software).
 *
 * Copyright 2026 INRIA.
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
#ifndef CXD_BROAD_SWEEP_HPP
#define CXD_BROAD_SWEEP_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "cxd_broad.hpp"

namespace cxd
{
  template<std::size_t Dim, typename Shapes, typename Padding>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase_sweep(const Shapes& shapes, const Padding padding)
  {
    using shape = std::decay_t<decltype(shapes[0])>;
    using scalar = typename decltype(traits::env(std::declval<shape>()))::scalar;
    const scalar padding_value = static_cast<scalar>(padding);
    if (padding_value < scalar{})
      throw std::invalid_argument("broad phase padding must be non-negative");

    struct proxy
    {
      std::size_t index;
      aabb<scalar, Dim> bounds;
    };

    std::vector<proxy> proxies;
    proxies.reserve(shapes.size());
    for (std::size_t index = 0; index < shapes.size(); ++index)
    {
      aabb<scalar, Dim> bounds = bounding_box<Dim>(shapes[index]);
      for (std::size_t d = 0; d < Dim; ++d)
      {
        bounds.minimum[d] -= padding_value;
        bounds.maximum[d] += padding_value;
      }
      proxies.push_back({index, bounds});
    }

    std::sort(proxies.begin(), proxies.end(),
              [](const proxy& first, const proxy& second)
              {
                if (first.bounds.minimum[0] != second.bounds.minimum[0])
                  return first.bounds.minimum[0] < second.bounds.minimum[0];
                return first.index < second.index;
              });

    std::vector<proxy> active;
    std::vector<std::pair<std::size_t, std::size_t>> candidates;
    for (const proxy& current : proxies)
    {
      std::erase_if(active,
                    [&current](const proxy& value)
                    {
                      return value.bounds.maximum[0] < current.bounds.minimum[0];
                    });

      for (const proxy& value : active)
      {
        bool overlap = true;
        for (std::size_t d = 0; d < Dim; ++d)
        {
          if (value.bounds.maximum[d] < current.bounds.minimum[d] ||
              current.bounds.maximum[d] < value.bounds.minimum[d])
          {
            overlap = false;
            break;
          }
        }
        if (overlap)
        {
          candidates.emplace_back(std::min(value.index, current.index),
                                  std::max(value.index, current.index));
        }
      }
      active.push_back(current);
    }
    return candidates;
  }

  template<typename Shapes, typename Padding>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase_sweep(const Shapes& shapes, const Padding padding)
  {
    return broad_phase_sweep<2>(shapes, padding);
  }

  template<typename Shapes>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase_sweep(const Shapes& shapes)
  {
    using shape = std::decay_t<decltype(shapes[0])>;
    using scalar = typename decltype(traits::env(std::declval<shape>()))::scalar;
    return broad_phase_sweep<2>(shapes, scalar{});
  }

}

#endif
