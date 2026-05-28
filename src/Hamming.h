#pragma once

#include <cassert>

#include "matrix/CSCMatrix.h"
#include "utils/ZipRows.h"

namespace similarity {

template <typename Field>
size_t hamming_fixed_ratio(const SparseVector<Field>& x,
                           const SparseVector<Field>& y, Field ratio) {
  size_t diff = 0;

  for (const auto [i, x, y] : SparseZipRange{x, y}) {
    if (FieldTraits<Field>::is_nonzero(x - y * ratio)) {
      ++diff;
    }
  }

  return diff;
}

template <typename Field>
bool hamming_fixed_ratio_leq(const SparseVector<Field>& x,
                             const SparseVector<Field>& y, Field ratio,
                             size_t max_distance) {
  size_t diff = 0;

  for (const auto [i, x, y] : SparseZipRange{x, y}) {
    if (FieldTraits<Field>::is_nonzero(x - y * ratio)) {
      ++diff;

      if (diff > max_distance) {
        return false;
      }
    }
  }

  return true;
}

// Returns (distance, ratio) if x - ratio * y has the least number of nonzero
// variables (this number is equal to distance).
template <typename Field>
std::pair<size_t, Field> hamming(const SparseVector<Field>& x,
                                 const SparseVector<Field>& y) {
  size_t distance = x.size() + y.size();
  Field min_distance_ratio = 1.;

  for (auto [alpha_index, x_value] : x) {
    Field y_value;
    bool found = false;

    for (auto [index, value] : y) {
      if (index == alpha_index) {
        y_value = value;
        found = true;
        break;
      }
    }

    if (!found) {
      continue;
    }

    const Field ratio = x_value / y_value;
    const size_t new_distance = hamming_fixed_ratio(x, y, ratio);

    if (new_distance < distance) {
      distance = new_distance;
      min_distance_ratio = ratio;
    }
  }

  return {distance, min_distance_ratio};
}

// If real hamming(x, y) = h, then:
// 1. if h * 2 < union(x, y), then the method returns h
// 2. otherwise something strange happens, TODO: think about it
template <typename Field>
size_t fast_hamming(const SparseVector<Field>& x,
                    const SparseVector<Field>& y) {
  size_t xor_size = 0;

  int balance = 0;
  Field current_ratio = 0;

  for (auto [i, x, y] : SparseZipRange{x, y}) {
    if (!FieldTraits<Field>::is_nonzero(x) ||
        !FieldTraits<Field>::is_nonzero(y)) {
      ++xor_size;
      continue;
    }

    const Field ratio = x / y;

    if (!FieldTraits<Field>::is_nonzero(ratio - current_ratio)) {
      ++balance;
    } else if (balance == 0) {
      current_ratio = ratio;
    } else {
      --balance;
    }
  }

  size_t diff = 0;

  for (auto [i, x, y] : SparseZipRange{x, y}) {
    if (!FieldTraits<Field>::is_nonzero(x) ||
        !FieldTraits<Field>::is_nonzero(y)) {
      continue;
    }

    if (FieldTraits<Field>::is_nonzero(x / y - current_ratio)) {
      ++diff;
    }
  }

  return diff + xor_size;
}

// Checks that d(a, b) <= max_distance. If true, returns \alpha such that
// | a - \alpha b | <= max_distance. Otherwise, returns std::nullopt.
// max_distance must be <= 8
template <typename Field>
std::optional<Field> hamming_leq(const SparseVector<Field>& x,
                                 const SparseVector<Field>& y,
                                 size_t max_distance) {
  constexpr size_t kRatiosCapacity = 8;
  assert(max_distance < kRatiosCapacity);

  std::array<std::pair<Field, size_t>, kRatiosCapacity> ratios;
  ratios.fill({0, 0});
  size_t xor_size = 0;
  size_t intersection_size = 0;

  for (auto [i, x, y] : SparseZipRange{x, y}) {
    if (!FieldTraits<Field>::is_nonzero(x) ||
        !FieldTraits<Field>::is_nonzero(y)) {
      ++xor_size;

      if (xor_size > max_distance) {
        return std::nullopt;
      }
    } else {
      ++intersection_size;

      const Field ratio = x / y;
      bool found_slot = false;

      for (size_t j = 0; j <= max_distance - xor_size; ++j) {
        if (!FieldTraits<Field>::is_nonzero(ratios[j].first)) {
          ratios[j].first = ratio;
          ratios[j].second = 1;

          found_slot = true;
          break;
        }
        if (!FieldTraits<Field>::is_nonzero(ratio - ratios[j].first)) {
          ++ratios[j].second;

          found_slot = true;
          break;
        }
      }

      if (!found_slot) {
        return std::nullopt;
      }
    }
  }

  size_t max_count = 0;
  Field max_count_ratio = 1;

  for (size_t i = 0; i <= max_distance - xor_size; ++i) {
    if (max_count < ratios[i].second) {
      max_count_ratio = ratios[i].first;
      max_count = ratios[i].second;
    }
  }

  if (xor_size + intersection_size - max_count > max_distance) {
    return std::nullopt;
  }

  return max_count_ratio;
}

}  // namespace similarity
