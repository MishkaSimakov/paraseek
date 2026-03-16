#include "Hamming.h"

size_t calculate_diff(const SparseVector<double>& x,
                      const SparseVector<double>& y, double ratio) {
  size_t y_index = 0;

  size_t diff = 0;

  for (size_t x_index = 0; x_index < x.size(); ++x_index) {
    while (y_index < y.size() && y[y_index].first < x[x_index].first) {
      ++y_index;
      ++diff;
    }

    if (y_index < y.size() && y[y_index].first == x[x_index].first) {
      if (FieldTraits<double>::is_nonzero(y[y_index].second * ratio -
                                          x[x_index].second)) {
        ++diff;
      }

      ++y_index;
    } else {
      ++diff;
    }
  }

  diff += y.size() - y_index;

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

    distance = std::min(distance, calculate_diff(x, y, x_value / y_value));
  }

  return distance;
}
