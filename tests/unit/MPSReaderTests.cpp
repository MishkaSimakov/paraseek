#include <gtest/gtest.h>

#include "helpers/CheckFeasible.h"
#include "helpers/Highs.h"
#include "problems/ArchivedProblem.h"

// Solves problem without paraseek problem representation
Solution full_highs_solve(const std::string& name) {
  Highs highs;
  highs.setOptionValue("output_flag", true);

  if (highs.readModel(paths::problem(name + ".mps")) != HighsStatus::kOk) {
    throw std::runtime_error("HiGHS failed to read model.");
  }

  const HighsLp& lp = highs.getLp();

  highs.passModel(lp);
  highs.run();

  Solution solution;
  solution.status = highs.getModelStatus();

  const auto& highs_sol = highs.getSolution();
  solution.x = highs_sol.col_value;
  solution.objective = highs.getInfo().objective_function_value;

  return solution;
}

TEST(MPSReaderTests, ReadAndSolve) {
  const std::vector<std::string> small_problems = {"markshare_4_0"};

  for (const auto& name : small_problems) {
    SCOPED_TRACE(std::format("problem name: {}", name));

    const auto problem = get_archived<double>(name);
    const auto solution = solve(problem);

    // full highs pipeline without my problem class
    const auto expected = full_highs_solve(name);

    // check solutions match
    ASSERT_EQ(solution.status, expected.status);

    // If infeasible/unbounded -> nothing more to check
    if (solution.status != HighsModelStatus::kOptimal) {
      return;
    }

    // Objective must match
    ASSERT_NEAR(solution.objective, expected.objective, 1e-6);

    check_feasible(problem, expected.x);
  }
}
