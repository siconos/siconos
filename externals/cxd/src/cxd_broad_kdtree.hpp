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
#ifndef CXD_BROAD_KDTREE_HPP
#define CXD_BROAD_KDTREE_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "cxd_aabb.hpp"

namespace cxd
{
  namespace kdtree
  {
    template<typename Scalar, std::size_t Dim>
    struct node
    {
      aabb<Scalar, Dim> bounds;
      std::vector<std::size_t> item_indices;
      std::unique_ptr<node> left;
      std::unique_ptr<node> right;
      bool is_leaf = false;
      std::size_t axis = 0;

      constexpr bool is_leaf_node() const { return is_leaf; }
    };

    template<typename Scalar, std::size_t Dim>
    struct kdtree
    {
      std::unique_ptr<node<Scalar, Dim>> root;
      std::size_t leaf_count = 0;
      std::size_t max_depth = 0;
      std::size_t max_leaf_size = 8;

      kdtree() = default;
      explicit kdtree(std::size_t max_leaf_size_) : max_leaf_size(max_leaf_size_) {}

      void build(std::vector<std::pair<std::size_t, aabb<Scalar, Dim>>>& items, std::size_t depth = 0)
      {
        if (items.empty())
          return;

        root = build_recursive(items, 0, items.size(), depth);
      }

      void refit(const std::vector<std::pair<std::size_t, aabb<Scalar, Dim>>>& items)
      {
        if (!root)
          return;
        refit_recursive(root.get(), items);
      }

    private:
      std::unique_ptr<node<Scalar, Dim>> build_recursive(
        std::vector<std::pair<std::size_t, aabb<Scalar, Dim>>>& items,
        std::size_t begin, std::size_t end, std::size_t depth)
      {
        max_depth = std::max(max_depth, depth);

        const std::size_t count = end - begin;
        if (count <= max_leaf_size)
        {
          auto leaf = std::make_unique<node<Scalar, Dim>>();
          leaf->is_leaf = true;
          leaf->bounds = items[begin].second;
          for (std::size_t i = begin + 1; i < end; ++i)
          {
            for (std::size_t d = 0; d < Dim; ++d)
            {
              leaf->bounds.minimum[d] = std::min(leaf->bounds.minimum[d], items[i].second.minimum[d]);
              leaf->bounds.maximum[d] = std::max(leaf->bounds.maximum[d], items[i].second.maximum[d]);
            }
          }
          leaf->item_indices.reserve(count);
          for (std::size_t i = begin; i < end; ++i)
            leaf->item_indices.push_back(i);
          std::sort(leaf->item_indices.begin(), leaf->item_indices.end(),
                    [&](std::size_t a, std::size_t b) { return items[a].first < items[b].first; });
          leaf_count = 1;
          return leaf;
        }

        const std::size_t axis = depth % Dim;
        auto mid = begin + count / 2;
        std::nth_element(items.begin() + begin, items.begin() + mid, items.begin() + end,
                         [axis](const auto& a, const auto& b)
                         {
                           return (a.second.minimum[axis] + a.second.maximum[axis]) <
                                  (b.second.minimum[axis] + b.second.maximum[axis]);
                         });

        using node_type = node<Scalar, Dim>;
        auto node = std::make_unique<node_type>();
        node->axis = axis;

        kdtree<Scalar, Dim> left_tree(max_leaf_size);
        kdtree<Scalar, Dim> right_tree(max_leaf_size);
        left_tree.root = left_tree.build_recursive(items, begin, mid, depth + 1);
        right_tree.root = right_tree.build_recursive(items, mid, end, depth + 1);

        node->left = std::move(left_tree.root);
        node->right = std::move(right_tree.root);
        
        for (std::size_t d = 0; d < Dim; ++d)
        {
          node->bounds.minimum[d] = std::min(node->left->bounds.minimum[d], node->right->bounds.minimum[d]);
          node->bounds.maximum[d] = std::max(node->left->bounds.maximum[d], node->right->bounds.maximum[d]);
        }
        
        leaf_count = left_tree.leaf_count + right_tree.leaf_count;
        max_depth = std::max(max_depth, std::max(left_tree.max_depth, right_tree.max_depth));
        return node;
      }

      void refit_recursive(node<Scalar, Dim>* n,
                           const std::vector<std::pair<std::size_t, aabb<Scalar, Dim>>>& items)
      {
        if (!n)
          return;

        if (n->is_leaf)
        {
          n->bounds = items[n->item_indices[0]].second;
          for (std::size_t i = 1; i < n->item_indices.size(); ++i)
          {
            for (std::size_t d = 0; d < Dim; ++d)
            {
              n->bounds.minimum[d] = std::min(n->bounds.minimum[d], items[n->item_indices[i]].second.minimum[d]);
              n->bounds.maximum[d] = std::max(n->bounds.maximum[d], items[n->item_indices[i]].second.maximum[d]);
            }
          }
          return;
        }

        refit_recursive(n->left.get(), items);
        refit_recursive(n->right.get(), items);

        for (std::size_t d = 0; d < Dim; ++d)
        {
          n->bounds.minimum[d] = std::min(n->left->bounds.minimum[d], n->right->bounds.minimum[d]);
          n->bounds.maximum[d] = std::max(n->left->bounds.maximum[d], n->right->bounds.maximum[d]);
        }
      }

    public:
      void query(const aabb<Scalar, Dim>& query_bounds,
                 std::vector<std::pair<std::size_t, std::size_t>>& results,
                 const std::vector<std::pair<std::size_t, aabb<Scalar, Dim>>>& all_items) const
      {
        query_recursive(root.get(), query_bounds, results, all_items);
      }

    private:
      void query_recursive(const node<Scalar, Dim>* n,
                           const aabb<Scalar, Dim>& query_bounds,
                           std::vector<std::pair<std::size_t, std::size_t>>& results,
                           const std::vector<std::pair<std::size_t, aabb<Scalar, Dim>>>& all_items) const
      {
        if (!n || !n->bounds.overlaps(query_bounds))
          return;

        if (n->is_leaf)
        {
          for (std::size_t idx : n->item_indices)
          {
            if (all_items[idx].second.overlaps(query_bounds))
              results.emplace_back(all_items[idx].first, all_items[idx].first);
          }
          return;
        }

        query_recursive(n->left.get(), query_bounds, results, all_items);
        query_recursive(n->right.get(), query_bounds, results, all_items);
      }
    };

    template<typename Scalar, std::size_t Dim>
    void collect_pairs_recursive(const node<Scalar, Dim>* a,
                                 const node<Scalar, Dim>* b,
                                 std::vector<std::pair<std::size_t, std::size_t>>& results,
                                 const std::vector<std::pair<std::size_t, aabb<Scalar, Dim>>>& all_items)
    {
      if (!a || !b || !a->bounds.overlaps(b->bounds))
        return;

      if (a->is_leaf && b->is_leaf)
      {
        if (a == b)
        {
          for (std::size_t ia = 0; ia < a->item_indices.size(); ++ia)
          {
            for (std::size_t ib = ia + 1; ib < b->item_indices.size(); ++ib)
            {
              const std::size_t pos_a = a->item_indices[ia];
              const std::size_t pos_b = b->item_indices[ib];
              const std::size_t idx_a = all_items[pos_a].first;
              const std::size_t idx_b = all_items[pos_b].first;
              if (idx_a >= idx_b)
                continue;
              if (all_items[pos_a].second.overlaps(all_items[pos_b].second))
                results.emplace_back(idx_a, idx_b);
            }
          }
        }
        else
        {
          for (std::size_t ia = 0; ia < a->item_indices.size(); ++ia)
          {
            for (std::size_t ib = 0; ib < b->item_indices.size(); ++ib)
            {
              const std::size_t pos_a = a->item_indices[ia];
              const std::size_t pos_b = b->item_indices[ib];
              const std::size_t idx_a = all_items[pos_a].first;
              const std::size_t idx_b = all_items[pos_b].first;
              if (idx_a >= idx_b)
                continue;
              if (all_items[pos_a].second.overlaps(all_items[pos_b].second))
                results.emplace_back(idx_a, idx_b);
            }
          }
        }
        return;
      }

      if (a->is_leaf)
      {
        collect_pairs_recursive(a, b->left.get(), results, all_items);
        collect_pairs_recursive(a, b->right.get(), results, all_items);
      }
      else if (b->is_leaf)
      {
        collect_pairs_recursive(a->left.get(), b, results, all_items);
        collect_pairs_recursive(a->right.get(), b, results, all_items);
      }
      else
      {
        collect_pairs_recursive(a->left.get(), b->left.get(), results, all_items);
        collect_pairs_recursive(a->left.get(), b->right.get(), results, all_items);
        collect_pairs_recursive(a->right.get(), b->left.get(), results, all_items);
        collect_pairs_recursive(a->right.get(), b->right.get(), results, all_items);
      }
    }
  }

  template<std::size_t Dim, typename Shapes, typename Padding>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase_kdtree_impl(const Shapes& shapes, const Padding padding)
  {
    using shape = std::decay_t<decltype(shapes[0])>;
    using scalar = typename decltype(traits::env(std::declval<shape>()))::scalar;
    const scalar padding_value = static_cast<scalar>(padding);
    if (padding_value < scalar{})
      throw std::invalid_argument("broad phase padding must be non-negative");

    std::vector<std::pair<std::size_t, aabb<scalar, Dim>>> items;
    items.reserve(shapes.size());
    for (std::size_t index = 0; index < shapes.size(); ++index)
    {
      aabb<scalar, Dim> bounds = bounding_box<Dim>(shapes[index]);
      for (std::size_t d = 0; d < Dim; ++d)
      {
        bounds.minimum[d] -= padding_value;
        bounds.maximum[d] += padding_value;
      }
      items.emplace_back(index, bounds);
    }

    kdtree::kdtree<scalar, Dim> tree(8);
    tree.build(items);
    tree.refit(items);

    std::vector<std::pair<std::size_t, std::size_t>> candidates;
    kdtree::collect_pairs_recursive(tree.root.get(), tree.root.get(), candidates, items);

    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

    return candidates;
  }

  template<std::size_t Dim, typename Shapes>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase_kdtree_impl(const Shapes& shapes)
  {
    using shape = std::decay_t<decltype(shapes[0])>;
    using scalar = typename decltype(traits::env(std::declval<shape>()))::scalar;
    return broad_phase_kdtree_impl<Dim>(shapes, scalar{});
  }

  template<typename Shapes, typename Padding>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase_kdtree(const Shapes& shapes, const Padding padding)
  {
    using point = std::decay_t<decltype(shapes[0][0])>;
    if constexpr (point{}.size() == 2)
      return broad_phase_kdtree_impl<2>(shapes, padding);
    else if constexpr (point{}.size() == 3)
      return broad_phase_kdtree_impl<3>(shapes, padding);
    else
      static_assert(std::false_type{}, "broad_phase_kdtree only supports 2D and 3D shapes");
  }

  template<typename Shapes>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase_kdtree(const Shapes& shapes)
  {
    using point = std::decay_t<decltype(shapes[0][0])>;
    if constexpr (point{}.size() == 2)
      return broad_phase_kdtree_impl<2>(shapes);
    else if constexpr (point{}.size() == 3)
      return broad_phase_kdtree_impl<3>(shapes);
    else
      static_assert(std::false_type{}, "broad_phase_kdtree only supports 2D and 3D shapes");
  }

  template<typename Shapes, typename Padding>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase(const Shapes& shapes, const Padding padding)
  {
    return broad_phase_kdtree(shapes, padding);
  }

  template<typename Shapes>
  requires requires (const Shapes& shapes)
  {
    { shapes.size() } -> std::convertible_to<std::size_t>;
    requires any_vectors<std::decay_t<decltype(shapes[0])>>;
  }
  auto broad_phase(const Shapes& shapes)
  {
    return broad_phase_kdtree(shapes);
  }
}

#endif
