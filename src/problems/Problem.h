#pragma once

#include <cassert>

#include "Bound.h"
#include "matrix/CSCMatrix.h"
#include "matrix/Matrix.h"

struct Problem {
  CSCMatrix<double> A;
  std::vector<double> b;
  std::vector<double> c;

  std::vector<Bound<double>> bounds;

  // Objective value shift. Objective is: min c^T x + shift
  double shift;

  Problem(CSCMatrix<double> A, std::vector<double> b, std::vector<double> c,
          std::vector<Bound<double>> bounds, double shift = 0)
      : A(std::move(A)),
        b(std::move(b)),
        c(std::move(c)),
        bounds(std::move(bounds)),
        shift(shift) {
    auto [n, d] = this->A.shape();

    assert(this->b.size() == n);
    assert(this->c.size() == d);
    assert(this->bounds.size() == d);
  }
};
