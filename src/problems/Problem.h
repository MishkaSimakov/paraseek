#pragma once

#include <cassert>

#include "Bound.h"
#include "matrix/CSCMatrix.h"
#include "matrix/Matrix.h"

struct UnsafeTag {};

// Problem is represented in the following format:
// min c^T x + shift, such that
// rhs_bounds.lower <= Ax <= rhs_bounds.upper
// bounds.lower <= x <= bounds.upper
// x_i is integer iff is_integer[i] = true
template <typename Field>
struct Problem {
  CSCMatrix<Field> A;
  std::vector<Bound<Field>> rhs_bounds;
  std::vector<Field> c;
  std::vector<Bound<Field>> bounds;

  std::vector<bool> is_integer;

  // Objective value shift. Objective is: min c^T x + shift
  Field shift;

  // During preprocessing passes problem may be proven to be unfeasible.
  bool proven_unfeasible = false;

  // unsafe constructor, doesn't check dimensions
  Problem(UnsafeTag /* tag */, CSCMatrix<Field> A,
          std::vector<Bound<Field>> rhs_bounds, std::vector<Field> c,
          std::vector<Bound<Field>> bounds, std::vector<bool> is_integer,
          Field shift = 0)
      : A(std::move(A)),
        rhs_bounds(std::move(rhs_bounds)),
        c(std::move(c)),
        bounds(std::move(bounds)),
        is_integer(std::move(is_integer)),
        shift(shift) {}

  static Problem with_size(size_t n, size_t d) {
    CSCMatrix<Field> A(n);
    std::vector<Bound<Field>> rhs_bounds(n);
    std::vector<Field> c(d, 0);
    std::vector<Bound<Field>> bounds(d);
    std::vector<bool> is_integer(d, false);

    return Problem(UnsafeTag{}, std::move(A), std::move(rhs_bounds),
                   std::move(c), std::move(bounds), std::move(is_integer));
  }

  Problem(CSCMatrix<Field> A, std::vector<Bound<Field>> rhs_bounds,
          std::vector<Field> c, std::vector<Bound<Field>> bounds,
          std::vector<bool> is_integer, Field shift = 0)
      : Problem(UnsafeTag{}, std::move(A), std::move(rhs_bounds), std::move(c),
                std::move(bounds), std::move(is_integer), shift) {
    auto [n, d] = this->A.shape();

    assert(this->rhs_bounds.size() == n);
    assert(this->c.size() == d);
    assert(this->bounds.size() == d);
    assert(this->is_integer.size() == d);
  }
};
