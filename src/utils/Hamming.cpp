#include "utils/Hamming.h"

#include <algorithm>
#include <cassert>

#include "ZipRows.h"

size_t similarity::hamming_fixed_ratio(const SparseVector<double>& x,
                                       const SparseVector<double>& y,
                                       double ratio) {
  size_t diff = 0;

  for (auto [i, x, y] : SparseZipRange{x, y}) {
    if (FieldTraits<double>::is_nonzero(x - y * ratio)) {
      ++diff;
    }
  }

  return diff;
}

bool similarity::hamming_fixed_ratio_leq(const SparseVector<double>& x,
                                         const SparseVector<double>& y,
                                         double ratio, size_t max_distance) {
  size_t diff = 0;

  for (auto [i, x, y] : SparseZipRange{x, y}) {
    if (FieldTraits<double>::is_nonzero(x - y * ratio)) {
      ++diff;

      if (diff > max_distance) {
        return false;
      }
    }
  }

  return true;
}

std::pair<size_t, double> similarity::hamming(const SparseVector<double>& x,
                                              const SparseVector<double>& y) {
  size_t distance = x.size() + y.size();
  double min_distance_ratio = 1.;

  for (auto [alpha_index, x_value] : x) {
    double y_value;
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

    const double ratio = x_value / y_value;
    const size_t new_distance = hamming_fixed_ratio(x, y, ratio);

    if (new_distance < distance) {
      distance = new_distance;
      min_distance_ratio = ratio;
    }
  }

  return {distance, min_distance_ratio};
}

size_t similarity::fast_hamming(const SparseVector<double>& x,
                                const SparseVector<double>& y) {
  const double inf = 1e20;

  int balance = 0;
  double current_ratio = 0;

  for (auto [i, x, y] : SparseZipRange{x, y}) {
    double ratio = y != 0 ? x / y : inf;

    if (!FieldTraits<double>::is_nonzero(ratio - current_ratio)) {
      ++balance;
    } else if (balance == 0) {
      current_ratio = ratio;
    } else {
      --balance;
    }
  }

  size_t diff = 0;

  for (auto [i, x, y] : SparseZipRange{x, y}) {
    const double ratio = y != 0 ? x / y : inf;

    if (FieldTraits<double>::is_nonzero(ratio - current_ratio)) {
      ++diff;
    }
  }

  return diff;
}

std::optional<double> similarity::hamming_leq(const SparseVector<double>& x,
                                              const SparseVector<double>& y,
                                              size_t max_distance) {
  constexpr size_t kRatiosCapacity = 8;
  assert(max_distance < kRatiosCapacity);

  std::array<std::pair<double, size_t>, kRatiosCapacity> ratios;
  ratios.fill({0, 0});
  size_t xor_size = 0;
  size_t intersection_size = 0;

  for (auto [i, x, y] : SparseZipRange{x, y}) {
    if (!FieldTraits<double>::is_nonzero(x) ||
        !FieldTraits<double>::is_nonzero(y)) {
      ++xor_size;

      if (xor_size > max_distance) {
        return std::nullopt;
      }
    } else {
      ++intersection_size;

      const double ratio = x / y;
      bool found_slot = false;

      for (size_t j = 0; j <= max_distance - xor_size; ++j) {
        if (!FieldTraits<double>::is_nonzero(ratios[j].first)) {
          ratios[j].first = ratio;
          ratios[j].second = 1;

          found_slot = true;
          break;
        }
        if (!FieldTraits<double>::is_nonzero(ratio - ratios[j].first)) {
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
  double max_count_ratio = 1;

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
