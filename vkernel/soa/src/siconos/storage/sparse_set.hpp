#pragma once

#include <type_traits>
#include <unordered_map>
#include <vector>

namespace siconos::storage {

template <typename Key, typename T>
struct sparse_set {
  using value_type = T;
  using key_type = Key;
  using size_type = std::size_t;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  static constexpr size_type invalid_index =
      static_cast<size_type>(-1);

  std::vector<T> dense_;
  std::vector<Key> dense_keys_;
  std::conditional_t<std::is_integral_v<Key>,
                     std::vector<size_type>,
                     std::unordered_map<Key, size_type>> sparse_;

  sparse_set() = default;

  bool contains(Key key) const
  {
    if constexpr (std::is_integral_v<Key>) {
      auto k = static_cast<size_type>(key);
      return k < sparse_.size() && sparse_[k] != invalid_index;
    } else {
      return sparse_.find(key) != sparse_.end();
    }
  }

  T& at(Key key)
  {
    if constexpr (std::is_integral_v<Key>) {
      return dense_.at(sparse_.at(static_cast<size_type>(key)));
    } else {
      return dense_.at(sparse_.at(key));
    }
  }

  const T& at(Key key) const
  {
    if constexpr (std::is_integral_v<Key>) {
      return dense_.at(sparse_.at(static_cast<size_type>(key)));
    } else {
      return dense_.at(sparse_.at(key));
    }
  }

  T* try_get(Key key)
  {
    if constexpr (std::is_integral_v<Key>) {
      auto k = static_cast<size_type>(key);
      if (k < sparse_.size() && sparse_[k] != invalid_index) {
        return &dense_[sparse_[k]];
      }
      return nullptr;
    } else {
      auto it = sparse_.find(key);
      return it != sparse_.end() ? &dense_[it->second] : nullptr;
    }
  }

  const T* try_get(Key key) const
  {
    if constexpr (std::is_integral_v<Key>) {
      auto k = static_cast<size_type>(key);
      if (k < sparse_.size() && sparse_[k] != invalid_index) {
        return &dense_[sparse_[k]];
      }
      return nullptr;
    } else {
      auto it = sparse_.find(key);
      return it != sparse_.end() ? &dense_[it->second] : nullptr;
    }
  }

  T& get_or_create(Key key)
  {
    if (auto* p = try_get(key)) return *p;
    return emplace(key);
  }

  // Non-mutating, always-valid lookup: returns a shared empty T when key is
  // absent instead of nullptr/throwing, so callers can index unconditionally
  // (e.g. `for (auto x : sparse[i])`) without a branch, matching how a
  // static per-item array already gives a possibly-empty T for every index.
  const T& operator[](Key key) const
  {
    if (auto* p = try_get(key)) return *p;
    static const T empty{};
    return empty;
  }

  template <typename... Args>
  T& emplace(Key key, Args&&... args)
  {
    if constexpr (std::is_integral_v<Key>) {
      auto k = static_cast<size_type>(key);
      if (k >= sparse_.size()) {
        sparse_.resize(k + 1, invalid_index);
      }
      auto idx = sparse_[k];
      if (idx != invalid_index) {
        dense_[idx] = T(std::forward<Args>(args)...);
        return dense_[idx];
      }
      idx = dense_.size();
      dense_.emplace_back(std::forward<Args>(args)...);
      dense_keys_.push_back(key);
      sparse_[k] = idx;
      return dense_.back();
    } else {
      auto it = sparse_.find(key);
      if (it != sparse_.end()) {
        dense_[it->second] = T(std::forward<Args>(args)...);
        return dense_[it->second];
      }
      auto idx = dense_.size();
      dense_.emplace_back(std::forward<Args>(args)...);
      dense_keys_.push_back(key);
      sparse_[key] = idx;
      return dense_.back();
    }
  }

  bool erase(Key key)
  {
    if constexpr (std::is_integral_v<Key>) {
      auto k = static_cast<size_type>(key);
      if (k >= sparse_.size() || sparse_[k] == invalid_index) {
        return false;
      }
      auto idx = sparse_[k];
      auto last_key = dense_keys_.back();
      dense_[idx] = std::move(dense_.back());
      dense_keys_[idx] = last_key;
      sparse_[static_cast<size_type>(last_key)] = idx;
      dense_.pop_back();
      dense_keys_.pop_back();
      sparse_[k] = invalid_index;
      return true;
    } else {
      auto it = sparse_.find(key);
      if (it == sparse_.end()) return false;
      auto idx = it->second;
      auto last_key = dense_keys_.back();
      dense_[idx] = std::move(dense_.back());
      dense_keys_[idx] = last_key;
      sparse_[last_key] = idx;
      dense_.pop_back();
      dense_keys_.pop_back();
      sparse_.erase(it);
      return true;
    }
  }

  iterator begin() { return dense_.begin(); }
  iterator end() { return dense_.end(); }
  const_iterator begin() const { return dense_.begin(); }
  const_iterator end() const { return dense_.end(); }
  const_iterator cbegin() const { return dense_.cbegin(); }
  const_iterator cend() const { return dense_.cend(); }

  struct kv_iterator {
    std::vector<T>* dense;
    std::vector<Key>* keys;
    size_type idx;
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::pair<const Key&, T&>;
    using difference_type = std::ptrdiff_t;
    value_type operator*() const { return {(*keys)[idx], (*dense)[idx]}; }
    kv_iterator& operator++()
    {
      ++idx;
      return *this;
    }
    bool operator==(const kv_iterator& o) const { return idx == o.idx; }
    bool operator!=(const kv_iterator& o) const { return idx != o.idx; }
  };

  kv_iterator kv_begin() { return {&dense_, &dense_keys_, 0}; }
  kv_iterator kv_end() { return {&dense_, &dense_keys_, dense_.size()}; }

  size_type size() const { return dense_.size(); }
  bool empty() const { return dense_.empty(); }
  void clear()
  {
    dense_.clear();
    dense_keys_.clear();
    sparse_.clear();
  }
  void reserve(size_type n)
  {
    dense_.reserve(n);
    dense_keys_.reserve(n);
    sparse_.reserve(n);
  }
};

}  // namespace siconos::storage

