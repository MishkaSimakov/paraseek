#include <gtest/gtest.h>

#include <random>
#include <unordered_set>

#include "helpers/RandomProblem.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Result.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

void compare_tables_with_brute_force(const CSCMatrix<double>& matrix,
                                     seekers::TablesParameters tables_params) {
  auto bf_result = seekers::BruteForce(2).seek(matrix);
  auto tbls_result = seekers::Tables(2, tables_params).seek(matrix);

  auto tbls_normalized = seekers::normalize_result(tbls_result).as_set();

  for (const auto p : bf_result) {
    auto itr = tbls_normalized.find(p);

    ASSERT_FALSE(itr == tbls_normalized.end());
    tbls_normalized.erase(itr);
  }

  ASSERT_TRUE(tbls_normalized.empty());
}

TEST(TablesTests, CompareWithBruteForce) {
  for (size_t i = 0; i < 100 && i < problems_names.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", problems_names[i]));

    auto matrix = get_problem_matrices(problems_names[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(
        matrix, {.groups_count = 4, .max_small_row_size = 4}));
  }
}

TEST(TablesTests, CompareWithBruteForceAllRowsSmall) {
  for (size_t i = 0; i < 100 && i < problems_names.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", problems_names[i]));

    auto matrix = get_problem_matrices(problems_names[i]).A;

    const seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = matrix.shape().second,
    };

    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(matrix, params));
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
  auto normalized = seekers::normalize_result(result);

  ASSERT_EQ(normalized.singular.size(), 1);
  ASSERT_EQ(normalized.singular[0], (std::pair{0, 1}));
}

TEST(TablesTests, SmallTest2) {
  CSCMatrix<double> matrix = {
      {1, 2, 1, 1},
      {1, 1, 3, 4},
  };

  auto result = seekers::Tables(2).seek(matrix);
  auto normalized = seekers::normalize_result(result);

  ASSERT_EQ(normalized.singular.size(), 0);
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
    auto normalized = seekers::normalize_result(result);

    if (similarity::hamming(matrix.get_row(0), matrix.get_row(1)).first <= 2) {
      ASSERT_EQ(normalized.singular.size(), 1);
      ASSERT_EQ(normalized.singular[0], (std::pair{0, 1}));
    } else {
      ASSERT_EQ(normalized.singular.size(), 0);
    }
  }
}

TEST(TablesTests, RandomizedTest2) {
  std::default_random_engine random;

  for (size_t i = 0; i < 1'000; ++i) {
    auto problem = generate_random_problem(50, 1000, random);

    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(
        problem.A, {.groups_count = 4, .max_small_row_size = 4}))
        << problem.A;
  }
}
