#include <gtest/gtest.h>

#include <random>

#include "helpers/CheckFeasible.h"
#include "helpers/FixedSplitter.h"
#include "helpers/Highs.h"
#include "helpers/Rational.h"
#include "nonzeros/ReduceNonzeros.h"
#include "problems/Random.h"

template <typename Field, typename FieldHasher>
void assert_reducer_correct(const Problem<Field>& problem) {
  const auto [n, d] = problem.A.shape();

  auto reduced = ReduceNonzeros<Field, FieldHasher>(6, 3).apply(problem);

  const auto solution_old = solve(problem);
  const auto solution_new = solve(reduced);

  // Status must match
  ASSERT_EQ(solution_old.status, solution_new.status);

  // If infeasible/unbounded -> nothing more to check
  if (solution_old.status != HighsModelStatus::kOptimal) {
    return;
  }

  // Objective must match
  ASSERT_NEAR(solution_old.objective, solution_new.objective, 1e-6);

  check_feasible(problem, solution_new.x);
}

template <typename Field>
Problem<Field> construct_problem(CSCMatrix<Field> matrix) {
  const auto [n, d] = matrix.shape();

  std::vector<Field> c(d, 0);
  std::vector<Bound<Field>> rhs_bounds(n, Bound<Field>{0, 0});
  std::vector<Bound<Field>> bounds(d, Bound<Field>{0, 0});
  std::vector<bool> is_integer(d, false);

  return Problem<Field>(std::move(matrix), std::move(rhs_bounds), std::move(c),
                        std::move(bounds), std::move(is_integer), 0);
}

using Field = Rational;
using FieldHasher = RationalHasher;

// using Field = double;
// using FieldHasher = seekers::DoubleHasher;

TEST(ReduceNonzerosTests, small_test_1) {
  CSCMatrix<Field> matrix = {
      {1, 1, 1, 1, 0},
      {1, 1, 1, 0, 0},
  };

  auto problem = construct_problem(matrix);
  auto split = FixedSplitter({{0, 1, 2}, {3, 4}});

  problem = ReduceNonzeros<Field, FieldHasher, FixedSplitter>(2, 1, split)
                .apply(problem);

  CSCMatrix<Field> expected = {
      {0, 0, 0, 1, 0},
      {1, 1, 1, 0, 0},
  };

  ASSERT_EQ(problem.A, expected);
}

TEST(ReduceNonzerosTests, small_test_2) {
  CSCMatrix<Field> matrix = {
      {2, 2, 2, 1, 0},
      {1, 1, 1, 0, 0},
  };

  auto problem = construct_problem(matrix);
  auto split = FixedSplitter({{0, 1, 2}, {3, 4}});

  problem = ReduceNonzeros<Field, FieldHasher, FixedSplitter>(2, 1, split)
                .apply(problem);

  CSCMatrix<Field> expected = {
      {0, 0, 0, 1, 0},
      {1, 1, 1, 0, 0},
  };

  ASSERT_EQ(problem.A, expected);
}

TEST(ReduceNonzerosTests, small_test_3) {
  const CSCMatrix<Field> matrix = {
      {2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0},
      {1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1},
  };

  auto problem = construct_problem(matrix);
  auto split = FixedSplitter({{0, 1, 2}, {3, 4, 5, 6, 7, 8, 9, 10}});

  problem = ReduceNonzeros<Field, FieldHasher, FixedSplitter>(2, 1, split)
                .apply(problem);

  ASSERT_EQ(problem.A, matrix);
}

TEST(ReduceNonzerosTests, test_random_small) {
  std::default_random_engine random;

  constexpr size_t rows_count = 3;

  for (size_t i = 0; i < 10'000; ++i) {
    auto problem = generate_random_problem<double>(rows_count, 5, random);

    ASSERT_NO_FATAL_FAILURE(
        (assert_reducer_correct<double, DoubleHasher>(problem)))
        << std::format("iteration = {}", i);
  }
}

TEST(ReduceNonzerosTests, test_random_big) {
  std::default_random_engine random;

  constexpr size_t rows_count = 100;

  for (size_t i = 0; i < 1000; ++i) {
    auto problem = generate_random_problem<double>(rows_count, 500, random);

    ASSERT_NO_FATAL_FAILURE(
        (assert_reducer_correct<double, DoubleHasher>(problem)))
        << std::format("iteration = {}", i);
  }
}
