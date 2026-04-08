#include <Highs.h>
#include <gtest/gtest.h>

#include <random>

#include "Reducer.h"
#include "helpers/Highs.h"
#include "helpers/LiftSolution.h"
#include "helpers/RandomProblem.h"
#include "problems/Problem.h"

void check_feasible(const Problem& problem, const std::vector<double>& x) {
  const auto [n, d] = problem.A.shape();

  ASSERT_EQ(x.size(), d);

  std::vector<double> real(n, 0);

  for (size_t col = 0; col < d; ++col) {
    for (const auto [row, value] : problem.A.get_column(col)) {
      real[row] += value * x[col];
    }
  }

  for (size_t row = 0; row < n; ++row) {
    ASSERT_NEAR(real[row], problem.b[row], 1e-6);
  }
}

void assert_reducer_correct(
    const Problem& problem,
    const std::vector<std::pair<size_t, size_t>>& rows) {
  const auto [n, d] = problem.A.shape();

  auto [reduced, mapping] = Reducer().apply(problem, rows);

  ASSERT_EQ(mapping.size(), d);

  // std::cout << problem.A << std::endl;
  // std::cout << reduced.A << std::endl;
  //
  // for (auto bound : problem.bounds) {
  //   std::cout << bound << " ";
  // }
  // std::cout << std::endl;
  //
  // for (auto bound : reduced.bounds) {
  //   std::cout << bound << " ";
  // }
  // std::cout << std::endl;
  //
  // for (auto b : problem.b) {
  //   std::cout << b << " ";
  // }
  // std::cout << std::endl;
  //
  // for (auto b : reduced.b) {
  //   std::cout << b << " ";
  // }
  // std::cout << std::endl;

  const auto solution_old = solve(problem);
  const auto solution_new = solve(reduced);

  // for (double value : solution_old.x) {
  //   std::cout << value << " ";
  // }
  // std::cout << std::endl;
  //
  // for (double value : solution_new.x) {
  //   std::cout << value << " ";
  // }
  // std::cout << std::endl;

  // Status must match
  ASSERT_EQ(solution_old.status, solution_new.status);

  // If infeasible/unbounded -> nothing more to check
  if (solution_old.status != HighsModelStatus::kOptimal) {
    return;
  }

  // Objective must match
  ASSERT_NEAR(solution_old.objective, solution_new.objective, 1e-6);

  // Lift new solution -> old space
  auto lifted = lift_solution(solution_new.x, mapping);

  check_feasible(problem, lifted);
}

TEST(ReducerTests, test_small_problem_1) {
  CSCMatrix<double> A = {{1, 0, 2, 3}, {0, 1, 1, 3}};

  std::vector<double> b = {0, 0};
  std::vector<double> c = {1, 2, 3, 4};

  std::vector<Bound<double>> bounds(4);

  const Problem problem(A, b, c, bounds);

  assert_reducer_correct(problem, {{0, 1}});
}

TEST(ReducerTests, test_small_problem_2) {
  CSCMatrix<double> A = {
      {1, 2, 0, 0},
      {2, 4, 1, 2},
      {0, 0, 3, 3},
  };

  std::vector<double> b = {1, 2, 3};
  std::vector<double> c = {1, 2, 3, 4};

  std::vector<Bound<double>> bounds(4);

  const Problem problem(A, b, c, bounds);

  ASSERT_NO_FATAL_FAILURE(assert_reducer_correct(problem, {{0, 1}}));
}

TEST(ReducerTests, test_random_small) {
  std::default_random_engine random;

  constexpr size_t rows_count = 2;
  std::vector<std::pair<size_t, size_t>> rows;

  // add all pairs into rows
  for (size_t i = 0; i < rows_count; ++i) {
    for (size_t j = i + 1; j < rows_count; ++j) {
      rows.emplace_back(i, j);
    }
  }

  for (size_t i = 0; i < 10'000; ++i) {
    auto problem = generate_random_problem(rows_count, 5, random);

    ASSERT_NO_FATAL_FAILURE(assert_reducer_correct(problem, rows));
  }
}

TEST(ReducerTests, test_random_big) {
  std::default_random_engine random;

  constexpr size_t rows_count = 25;
  std::vector<std::pair<size_t, size_t>> rows;

  // add all pairs into rows
  for (size_t i = 0; i < rows_count; ++i) {
    for (size_t j = i + 1; j < rows_count; ++j) {
      rows.emplace_back(i, j);
    }
  }

  for (size_t i = 0; i < 100; ++i) {
    auto problem = generate_random_problem(rows_count, 1000, random);

    ASSERT_NO_FATAL_FAILURE(assert_reducer_correct(problem, rows));
  }
}
