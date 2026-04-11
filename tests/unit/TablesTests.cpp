#include <gtest/gtest.h>

#include <random>

#include "helpers/RandomProblem.h"
#include "helpers/Rational.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Result.h"
#include "seekers/Tables.h"
#include "utils/Printing.h"

void compare_tables_with_brute_force(const CSCMatrix<Rational>& matrix,
                                     size_t max_diff,
                                     seekers::TablesParameters tables_params) {
  auto bf_result = seekers::BruteForce<Rational>(max_diff).seek(matrix);
  auto tbls_result =
      seekers::Tables<Rational, RationalHasher>(max_diff, tables_params)
          .seek(matrix);

  auto tbls_set = seekers::normalize_result(tbls_result).as_set();
  auto bf_set = std::set(bf_result.begin(), bf_result.end());

  // for (auto [i, j] : tbls_set) {
  //   if (similarity::hamming(matrix.get_row(i), matrix.get_row(j)).first >
  //   max_diff) {
  //     std::println("error!");
  //     similarity::hamming_leq(matrix.get_row(i), matrix.get_row(j), 3);
  //     similarity::hamming(matrix.get_row(i), matrix.get_row(j));
  //   }
  // }

  ASSERT_EQ(tbls_set, bf_set);
}

TEST(TablesTests, CompareWithBruteForce1) {
  for (size_t i = 0; i < 50 && i < problems_names.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", problems_names[i]));

    auto matrix = get_problem<Rational>(problems_names[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(
        matrix, 1, {.groups_count = 2, .max_small_row_size = 4}));
  }
}

TEST(TablesTests, CompareWithBruteForce2) {
  for (size_t i = 0; i < 50 && i < problems_names.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", problems_names[i]));

    auto matrix = get_problem<Rational>(problems_names[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(
        matrix, 2, {.groups_count = 4, .max_small_row_size = 4}));
  }
}

TEST(TablesTests, CompareWithBruteForce3) {
  for (size_t i = 0; i < 50 && i < problems_names.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", problems_names[i]));

    auto matrix = get_problem<Rational>(problems_names[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(
        matrix, 3, {.groups_count = 6, .max_small_row_size = 6}));
  }
}

TEST(TablesTests, CompareWithBruteForceAllRowsSmall) {
  for (size_t i = 0; i < 50 && i < problems_names.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", problems_names[i]));

    auto matrix = get_problem<Rational>(problems_names[i]).A;

    const seekers::TablesParameters params{
        .groups_count = 4,
        .max_small_row_size = matrix.shape().second,
    };

    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(matrix, 2, params));
  }
}

TEST(TablesTests, SmallTest1) {
  CSCMatrix<Rational> matrix = {
      {1, 1, 1, 1},
      {1, 1, 1, 1},
  };

  seekers::TablesParameters params{
      .max_small_row_size = 100,
      .groups_count = 4,
  };

  auto result =
      seekers::Tables<Rational, RationalHasher>(2, params).seek(matrix);
  auto normalized = seekers::normalize_result(result);

  ASSERT_EQ(normalized.singular.size(), 1);
  ASSERT_EQ(normalized.singular[0], (std::pair{0, 1}));
}

TEST(TablesTests, SmallTest2) {
  CSCMatrix<Rational> matrix = {
      {1, 2, 1, 1},
      {1, 1, 3, 4},
  };

  auto result = seekers::Tables<Rational, RationalHasher>(2).seek(matrix);
  auto normalized = seekers::normalize_result(result);

  ASSERT_EQ(normalized.singular.size(), 0);
}

TEST(TablesTests, RandomizedTest1) {
  constexpr size_t rows_size = 4;

  std::default_random_engine engine;
  std::uniform_int_distribution<int> elements_distribution(-5, 5);

  for (size_t i = 0; i < 100; ++i) {
    CSCMatrix<Rational> matrix(2);

    for (size_t j = 0; j < rows_size; ++j) {
      matrix.add_column();

      matrix.push_to_last_column(0, elements_distribution(engine));
      matrix.push_to_last_column(1, elements_distribution(engine));
    }

    seekers::TablesParameters params{
        .max_small_row_size = 100,
        .groups_count = 4,
    };

    auto result =
        seekers::Tables<Rational, RationalHasher>(2, params).seek(matrix);
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

  for (size_t i = 0; i < 100'000; ++i) {
    auto problem = generate_random_problem<Rational>(5, 10, random);

    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(
        problem.A, 2, {.groups_count = 4, .max_small_row_size = 4}))
        << linalg::to_dense(problem.A);
  }
}

TEST(TablesTests, RandomizedTest3) {
  std::default_random_engine random;

  for (size_t i = 0; i < 1'000; ++i) {
    auto problem = generate_random_problem<Rational>(50, 1000, random);

    ASSERT_NO_FATAL_FAILURE(compare_tables_with_brute_force(
        problem.A, 2, {.groups_count = 4, .max_small_row_size = 4}))
        << linalg::to_dense(problem.A) << " " << i;
  }
}
