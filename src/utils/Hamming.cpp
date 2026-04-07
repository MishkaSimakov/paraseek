#include "utils/Hamming.h"

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

    if (!FieldTraits<double>::is_nonzero(ratio - current_ratio)) {
      ++diff;
    }
  }

  return diff;
}
