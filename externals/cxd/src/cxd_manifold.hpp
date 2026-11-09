/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
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
#ifndef CXD_MANIFOLD_HPP
#define CXD_MANIFOLD_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <unordered_map>
#include <vector>

namespace cxd
{
  template<std::size_t Dim, std::size_t MaxContactsPerPair>
  class contact_manifold
  {
  public:
    struct key_type
    {
      std::size_t first;
      std::size_t second;

      auto operator<=>(const key_type&) const = default;
    };

    struct contact
    {
      std::array<float, Dim> normal{};
      float penetration = 0.F;
      std::array<float, Dim> location{};
    };

    void update(key_type key, std::span<const contact> contacts)
    {
      if (contacts.empty())
      {
        erase(key);
        return;
      }

      auto& entry = map_[key];
      entry.count = std::min(contacts.size(), MaxContactsPerPair);
      std::copy_n(contacts.begin(), entry.count, entry.contacts.begin());
      touch(key);
    }

    std::span<const contact> get(key_type key) const
    {
      if (auto it = map_.find(key); it != map_.end())
        return {it->second.contacts.data(), it->second.count};
      return {};
    }

    bool contains(key_type key) const
    {
      return map_.contains(key);
    }

    void erase(key_type key)
    {
      if (map_.erase(key))
        touch(key);
    }

    void clear()
    {
      map_.clear();
      dirty_.clear();
    }

    std::span<const key_type> changed_keys() const
    {
      return {dirty_.data(), dirty_.size()};
    }

    void mark_synced()
    {
      dirty_.clear();
    }

    static constexpr std::size_t capacity_per_pair()
    {
      return MaxContactsPerPair;
    }

    std::size_t size() const
    {
      return map_.size();
    }

    bool empty() const
    {
      return map_.empty();
    }

    template<typename F>
    requires std::invocable<F, key_type, std::span<const contact>>
    void for_each(F&& f) const
    {
      for (const auto& [key, entry] : map_)
        f(key, std::span{entry.contacts.data(), entry.count});
    }

    std::vector<key_type> keys() const
    {
      std::vector<key_type> result;
      result.reserve(map_.size());
      for (const auto& [key, _] : map_)
        result.push_back(key);
      return result;
    }

  private:
    struct pair_data
    {
      std::array<contact, MaxContactsPerPair> contacts;
      std::size_t count = 0;
    };

    struct key_hash
    {
      std::size_t operator()(const key_type& key) const noexcept
      {
        std::size_t h1 = std::hash<std::size_t>{}(key.first);
        std::size_t h2 = std::hash<std::size_t>{}(key.second);
        return h1 ^ (h2 << 1);
      }
    };

    void touch(key_type key)
    {
      if (std::find(dirty_.begin(), dirty_.end(), key) == dirty_.end())
        dirty_.push_back(key);
    }

    std::unordered_map<key_type, pair_data, key_hash> map_;
    std::vector<key_type> dirty_;
  };
}

#endif
