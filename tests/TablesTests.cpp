#include <gtest/gtest.h>

#include "../src/problems/ProblemMatrix.h"
#include "../src/problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"

TEST(TablesTests, CompareWithBruteForce) {
  for (size_t i = 0; i < 50 && i < problems_names.size(); ++i) {
    auto matrix = get_problem_matrix(problems_names[i]);

    auto bf_result = seekers::BruteForce(2).seek(matrix);
    auto tbls_result = seekers::Tables(2).seek(matrix);

    auto tbls_normalized =
        seekers::normalize_tables_result(tbls_result.first, tbls_result.second);

    for (auto p : bf_result) {
      auto itr = tbls_normalized.find(p);

      ASSERT_FALSE(itr == tbls_normalized.end())
          << std::format("failed for problem: {}", problems_names[i]);
      tbls_normalized.erase(itr);
    }

    ASSERT_TRUE(tbls_normalized.empty())
        << std::format("failed for problem: {}", problems_names[i]);
  }
}
