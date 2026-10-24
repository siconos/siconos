#pragma once

#include <array>
#include <cstddef>
#include <functional>

// std::hash for the std::array keys of the contact maps (ds pair, fixed disk
// centre, segment coefficients, with or without a manifold point index)
namespace std {

// std::array<T, 2>
template <typename T>
struct hash<std::array<T, 2>> {
  size_t operator()(const std::array<T, 2>& arr) const noexcept {
    size_t h1 = std::hash<T>{}(arr[0]);
    size_t h2 = std::hash<T>{}(arr[1]);
    // boost::hash_combine
    return h1 ^ (h2 + 0x9e3779b97f4a7c15 + (h1 << 6) + (h1 >> 2));
  }
};

// std::array<T, 4>
template <typename T>
struct hash<std::array<T, 4>> {
  size_t operator()(const std::array<T, 4>& arr) const noexcept {
    size_t seed = 0;
    for (size_t i = 0; i < 4; ++i) {
      size_t h = std::hash<T>{}(arr[i]);
      seed ^= h + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

template <typename T>
struct hash<std::array<T, 3>> {
  size_t operator()(const std::array<T, 3>& arr) const noexcept {
    size_t seed = 0;
    for (size_t i = 0; i < 3; ++i) {
      size_t h = std::hash<T>{}(arr[i]);
      seed ^= h + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

template <typename T>
struct hash<std::array<T, 5>> {
  size_t operator()(const std::array<T, 5>& arr) const noexcept {
    size_t seed = 0;
    for (size_t i = 0; i < 5; ++i) {
      size_t h = std::hash<T>{}(arr[i]);
      seed ^= h + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

}  // namespace std
