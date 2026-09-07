#pragma once

#include <any>
#include <unordered_map>
#include <stdexcept>

namespace siconos::storage {

template <typename Key>
struct dynamic_properties {
  using key_t = Key;
  std::unordered_map<Key, const std::type_info*> _type_registry;
  using map_t = std::unordered_map<Key, std::any>;
  map_t _data;

  template <typename Value>
  Value& get_or_create(const Key& key) {
    auto it = _type_registry.find(key);
    if (it != _type_registry.end() && *it->second != typeid(Value)) {
      throw std::runtime_error("Type mismatch for key: " + key);
    }
    _type_registry[key] = &typeid(Value);
    auto& any = _data[key];
    if (!any.has_value()) any = Value{};
    return std::any_cast<Value&>(any);
  }

  template <typename Value>
  struct proxy {
    dynamic_properties& storage;
    Key key;

    // Implicit conversion TO Value& (allows: Value& v = storage["key"])
    operator Value&() { return storage.get_or_create<Value>(key); }

    // Assignment from Value
    proxy& operator=(const Value& v)
    {
      storage.get_or_create<Value>(key) = v;
      return *this;
    }

    // Member access
    Value* operator->() { return &storage.get_or_create<Value>(key); }
    Value& operator*() { return storage.get_or_create<Value>(key); }
  };

  // operator[] returns a GENERIC proxy that works for any Value
  // Use a type-erased proxy that converts on use
  struct generic_proxy {
    dynamic_properties& storage;
    Key key;

    // Allow: auto& v = storage["key"]; (deduces from LHS)
    template <typename Value>
    operator Value&()
    {
      return storage.get_or_create<Value>(key);
    }

    // Allow: storage["key"] = value; (deduces from RHS)
    template <typename Value>
    generic_proxy& operator=(const Value& v)
    {
      storage.get_or_create<Value>(key) = v;
      return *this;
    }

    // Allow: storage["key"].method(); (deduces from member access)
    template <typename Value>
    Value* operator->()
    {
      return &storage.get_or_create<Value>(key);
    }
  };

  generic_proxy operator[](const Key& key) { return {*this, key}; }


  template <typename Value>
  Value& get(const Key& key)
  {
    return get_or_create<Value>(key);
  }

  template <typename Value>
  const Value& get(const Key& key) const
  {
    auto it = _data.find(key);
    if (it == _data.end()) {
      throw std::out_of_range("Key not found in dynamic_properties");
    }
    return std::any_cast<const Value&>(it->second);
  }

  template <typename Value>
  Value* try_get(const Key& key)
  {
    auto it = _data.find(key);
    return it != _data.end() ? std::any_cast<Value>(&it->second) : nullptr;
  }

  template <typename Value>
  const Value* try_get(const Key& key) const
  {
    auto it = _data.find(key);
    return it != _data.end() ? std::any_cast<const Value>(&it->second)
                             : nullptr;
  }

  bool contains(const Key& key) const { return _data.contains(key); }
  bool erase(const Key& key) { return _data.erase(key); }
  auto begin() { return _data.begin(); }
  auto end() { return _data.end(); }
  size_t size() const { return _data.size(); }
  void clear() { _data.clear(); }
};
}  // namespace siconos::storage
