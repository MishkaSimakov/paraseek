#include "Hamming.h"

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

size_t similarity::hamming(const SparseVector<double>& x,
                           const SparseVector<double>& y) {
  size_t distance = x.size() + y.size();

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

    distance = std::min(distance, hamming_fixed_ratio(x, y, x_value / y_value));
  }

  return distance;
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
    double ratio = y != 0 ? x / y : inf;

    if (ratio != current_ratio) {
      ++diff;
    }
  }

  return diff;
}
