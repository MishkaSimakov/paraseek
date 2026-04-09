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

  // During preprocessing passes problem may be proven to be unfeasible.
  bool proven_unfeasible = false;

  static Problem with_size(size_t n, size_t d) {
    CSCMatrix<double> A(n);
    std::vector<double> b(n);
    std::vector<double> c(d, 0);
    std::vector<Bound<double>> bounds(d);

    return Problem(std::move(A), std::move(b), std::move(c), std::move(bounds));
  }

  Problem(CSCMatrix<double> A, std::vector<double> b, std::vector<double> c,
          std::vector<Bound<double>> bounds, double shift = 0)
      : A(std::move(A)),
        b(std::move(b)),
        c(std::move(c)),
        bounds(std::move(bounds)),
        shift(shift) {
    // auto [n, d] = this->A.shape();
    //
    // assert(this->b.size() == n);
    // assert(this->c.size() == d);
    // assert(this->bounds.size() == d);
  }
};
