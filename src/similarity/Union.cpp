#include "Union.h"

size_t similarity::union_size(const SparseVector<double>& x,
                              const SparseVector<double>& y) {
  size_t y_index = 0;
  size_t result = 0;

  for (size_t x_index = 0; x_index < x.size(); ++x_index) {
    while (y_index < y.size() && y[y_index].first < x[x_index].first) {
      ++y_index;
      ++result;
    }

    if (y_index < y.size() && y[y_index].first == x[x_index].first) {
      ++y_index;
    }

    ++result;
  }

  result += y.size() - y_index;

  return result;
}
