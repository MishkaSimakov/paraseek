#include "Jaccard.h"

#include "Hamming.h"
#include "Union.h"

double similarity::jaccard(const SparseVector<double>& x,
                           const SparseVector<double>& y) {
  return 1. - static_cast<double>(hamming(x, y)) /
                  static_cast<double>(union_size(x, y));
}
