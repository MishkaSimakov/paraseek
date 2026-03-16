#pragma once

#include "matrix/CSCMatrix.h"

namespace similarity {

size_t hamming(const SparseVector<double>& x, const SparseVector<double>& y);

}
