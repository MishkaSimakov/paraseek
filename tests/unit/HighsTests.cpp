#include <Highs.h>
#include <gtest/gtest.h>

#include "helpers/Highs.h"
#include "problems/Problem.h"

TEST(ReducerTests, CheckHighs) {
  // Constraints:
  // x + 2y <= 4
  // 3x + y <= 5
  CSCMatrix<double> A = {{1, 2, 1, 0}, {3, 1, 0, 1}};

  std::vector<Bound<double>> rhs_bounds = {
      {4, 4},
      {5, 5},
  };

  // Objective: maximize x + y
  std::vector<double> c = {-1, -1, 0, 0};

  // Bounds: x, y >= 0
  std::vector<Bound<double>> bounds = {
      {0, std::nullopt},
      {0, std::nullopt},
      {0, std::nullopt},
      {0, std::nullopt},
  };

  Problem problem(A, rhs_bounds, c, bounds);

  auto solution = solve(problem);

  ASSERT_DOUBLE_EQ(solution.x[0], 6. / 5.);
  ASSERT_DOUBLE_EQ(solution.x[1], 7. / 5.);

  ASSERT_DOUBLE_EQ(-solution.objective, 13. / 5.);

  ASSERT_EQ(solution.status, HighsModelStatus::kOptimal);
}
