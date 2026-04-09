#pragma once

#include <random>

#include "problems/Problem.h"
#include "utils/Random.h"

// If feasible = false, then returned problem may not be feasible.
template <typename Gen>
  requires std::uniform_random_bit_generator<std::remove_reference_t<Gen>>
Problem generate_random_problem(size_t n, size_t d, Gen&& generator,
                                bool feasible = true) {
  std::uniform_real_distribution<double> alpha_dist(0.1, 10);
  std::uniform_real_distribution<double> value_dist(-5.0, 5.0);
  std::uniform_real_distribution<double> prob(0.0, 1.0);

  // --- Generate random feasible solution ---
  std::vector<double> x_true(d);
  for (size_t j = 0; j < d; ++j) {
    x_true[j] = value_dist(generator);
  }

  // --- Generate sparse A ---
  Matrix<double> dense(n, d);

  size_t base_rows = std::max(1uz, n / 2);

  // 1. Generate base rows
  for (size_t i = 0; i < base_rows; ++i) {
    for (size_t j = 0; j < d; ++j) {
      if (prob(generator) < 0.4) {
        dense[i, j] = value_dist(generator);
      }
    }
  }

  // 2. Generate dependent rows
  constexpr size_t max_diff = 2;

  for (size_t i = base_rows; i < n; ++i) {
    const size_t ref = generator() % base_rows;
    double alpha = alpha_dist(generator);
    if (prob(generator) < 0.5) {
      alpha *= -1;
    }

    for (size_t j = 0; j < d; ++j) {
      dense[i, j] = alpha * dense[ref, j];
    }

    for (size_t j = 0; j < max_diff; ++j) {
      const size_t col = generator() % d;

      dense[i, col] += value_dist(generator);
    }
  }

  CSCMatrix<double> A(dense);

  // --- Compute b = A * x_true ---
  std::vector<double> b(n, 0.0);

  for (size_t col = 0; col < d; ++col) {
    for (auto [row, value] : A.get_column(col)) {
      b[row] += value * x_true[col];
    }
  }

  if (!feasible) {
    for (size_t row = 0; row < n; ++row) {
      if (prob(generator) < 0.5) {
        b[row] += value_dist(generator);
      }
    }
  }

  // --- Random objective ---
  std::vector<double> c(d);
  for (size_t j = 0; j < d; ++j) {
    c[j] = value_dist(generator);
  }

  // --- Bounds (loose but valid) ---
  std::vector<Bound<double>> bounds(d);

  std::uniform_real_distribution<double> bounds_dist(-10., 10.);

  for (size_t j = 0; j < d; ++j) {
    double lower = bounds_dist(generator);
    double upper = bounds_dist(generator);

    if (lower > upper) {
      std::swap(lower, upper);
    }

    bounds[j].lower = lower;
    bounds[j].upper = upper;
  }

  return {A, b, c, bounds};
}
