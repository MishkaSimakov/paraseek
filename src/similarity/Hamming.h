#pragma once

#include "matrix/CSCMatrix.h"

namespace similarity {

size_t hamming_fixed_ratio(const SparseVector<double>& x,
                           const SparseVector<double>& y, double ratio);

size_t hamming(const SparseVector<double>& x, const SparseVector<double>& y);

// If real hamming(x, y) = h, then:
// 1. if h * 2 < union(x, y), then the method returns h
// 2. otherwise something strange happens, TODO: think about it
size_t fast_hamming(const SparseVector<double>& x,
                    const SparseVector<double>& y);

}  // namespace similarity
