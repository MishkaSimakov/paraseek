#pragma once

#include <cassert>

#include "Bound.h"
#include "matrix/CSCMatrix.h"
#include "matrix/Matrix.h"

// Problem is represented in the following format:
// min c^T x + shift, such that
// rhs_bounds.lower <= Ax <= rhs_bounds.upper
// bounds.lower <= x <= bounds.upper
template <typename Field>
struct Problem {
  CSCMatrix<Field> A;
  std::vector<Field> b;
  std::vector<Field> c;

  std::vector<Bound<Field>> bounds;

  // Objective value shift. Objective is: min c^T x + shift
  Field shift;

  // During preprocessing passes problem may be proven to be unfeasible.
  bool proven_unfeasible = false;

  // unsafe constructor, doesn't check dimensions
  Problem(std::monostate, CSCMatrix<Field> A, std::vector<Field> b,
          std::vector<Field> c, std::vector<Bound<Field>> bounds,
          Field shift = 0)
      : A(std::move(A)),
        b(std::move(b)),
        c(std::move(c)),
        bounds(std::move(bounds)),
        shift(shift) {}

 public:
  static Problem with_size(size_t n, size_t d) {
    CSCMatrix<Field> A(n);
    std::vector<Field> b(n);
    std::vector<Field> c(d, 0);
    std::vector<Bound<Field>> bounds(d);

    return unsafe(std::move(A), std::move(b), std::move(c), std::move(bounds));
  }

  // Doesn't check dimensions
  static Problem unsafe(CSCMatrix<Field> A, std::vector<Field> b,
                        std::vector<Field> c, std::vector<Bound<Field>> bounds,
                        Field shift = 0) {
    return Problem(std::monostate{}, std::move(A), std::move(b), std::move(c),
                   std::move(bounds), shift);
  }

  Problem(CSCMatrix<Field> A, std::vector<Field> b, std::vector<Field> c,
          std::vector<Bound<Field>> bounds, Field shift = 0)
      : Problem(std::monostate{}, std::move(A), std::move(b), std::move(c),
                std::move(bounds), shift) {
    auto [n, d] = this->A.shape();

    assert(this->b.size() == n);
    assert(this->c.size() == d);
    assert(this->bounds.size() == d);
  }
};
