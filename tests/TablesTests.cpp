#include <gtest/gtest.h>

#include <random>

#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

TEST(TablesTests, CompareWithBruteForce) {
  for (size_t i = 0; i < 100 && i < problems_names.size(); ++i) {
    auto matrix = get_problem_matrix(problems_names[i]);

    auto bf_result = seekers::BruteForce(2).seek(matrix);
    auto tbls_result = seekers::Tables(2).seek(matrix);

    auto tbls_normalized =
        seekers::normalize_tables_result(tbls_result.first, tbls_result.second);
    std::unordered_set tbls_normalized_set(tbls_normalized.begin(),
                                           tbls_normalized.end());

    for (auto p : bf_result) {
      auto itr = tbls_normalized_set.find(p);

      ASSERT_FALSE(itr == tbls_normalized_set.end())
          << std::format("failed for problem: {}", problems_names[i]);
      tbls_normalized_set.erase(itr);
    }

    ASSERT_TRUE(tbls_normalized_set.empty())
        << std::format("failed for problem: {}", problems_names[i]);
  }
}

TEST(TablesTests, SmallTest1) {
  CSCMatrix<double> matrix = {
      {1, 1, 1, 1},
      {1, 1, 1, 1},
  };

  seekers::TablesParameters params{
      .max_small_row_size = 100,
      .groups_count = 4,
  };

  auto result = seekers::Tables(2, params).seek(matrix);
  auto normalized =
      seekers::normalize_tables_result(result.first, result.second);

  ASSERT_EQ(normalized.size(), 1);
  ASSERT_EQ(normalized[0], (std::pair{0, 1}));
}

TEST(TablesTests, SmallTest2) {
  CSCMatrix<double> matrix = {
      {1, 2, 1, 1},
      {1, 1, 3, 4},
  };

  auto result = seekers::Tables(2).seek(matrix);
  auto normalized =
      seekers::normalize_tables_result(result.first, result.second);

  ASSERT_EQ(normalized.size(), 0);
}

TEST(TablesTests, RandomizedSmallTest) {
  constexpr size_t rows_size = 4;

  std::default_random_engine engine;
  std::uniform_int_distribution<int> elements_distribution(-5, 5);

  for (size_t i = 0; i < 100; ++i) {
    CSCMatrix<double> matrix(2);

    for (size_t j = 0; j < rows_size; ++j) {
      matrix.add_column();

      matrix.push_to_last_column(0, elements_distribution(engine));
      matrix.push_to_last_column(1, elements_distribution(engine));
    }

    seekers::TablesParameters params{
        .max_small_row_size = 100,
        .groups_count = 4,
    };

    auto result = seekers::Tables(2, params).seek(matrix);
    auto normalized =
        seekers::normalize_tables_result(result.first, result.second);

    if (similarity::hamming(matrix.get_row(0), matrix.get_row(1)) <= 2) {
      ASSERT_EQ(normalized.size(), 1);
      ASSERT_EQ(normalized[0], (std::pair{0, 1}));
    } else {
      ASSERT_EQ(normalized.size(), 0);
    }
  }
}
