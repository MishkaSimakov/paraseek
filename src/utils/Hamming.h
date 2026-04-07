#pragma once

#include "matrix/CSCMatrix.h"

namespace similarity {

size_t hamming_fixed_ratio(const SparseVector<double>& x,
                           const SparseVector<double>& y, double ratio);

bool hamming_fixed_ratio_leq(const SparseVector<double>& x,
                               const SparseVector<double>& y, double ratio,
                               size_t max_distance);

// Returns (distance, ratio) if x - ratio * y has the least number of nonzero
// variables (this number is equal to distance).
std::pair<size_t, double> hamming(const SparseVector<double>& x,
                                  const SparseVector<double>& y);

// If real hamming(x, y) = h, then:
// 1. if h * 2 < union(x, y), then the method returns h
// 2. otherwise something strange happens, TODO: think about it
size_t fast_hamming(const SparseVector<double>& x,
                    const SparseVector<double>& y);

// Checks that d(a, b) <= max_distance. If true, returns \alpha such that
// | a - \alpha b | <= max_distance. Otherwise, returns std::nullopt.
// max_distance must be <= 8
std::optional<double> hamming_leq(const SparseVector<double>& x,
                                  const SparseVector<double>& y,
                                  size_t max_distance);

}  // namespace similarity
