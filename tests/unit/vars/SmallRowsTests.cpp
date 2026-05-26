#include <gtest/gtest.h>

#include <random>
#include <utility>

#include "helpers/RandomProblem.h"
#include "helpers/Rational.h"
#include "problems/ArchivedProblem.h"
#include "problems/ProblemsNames.h"
#include "utils/Printing.h"
#include "vars/BruteForce.h"
#include "vars/Result.h"
#include "vars/SmallRows.h"

void compare_with_brute_force(const CSCMatrix<Rational>& matrix,
                              size_t max_diff) {
  seekers::SmallRowsParameters params{
      .max_diff = max_diff,
      .max_row_size = matrix.shape().second,
      .entries_reduction = true,
  };

  const auto expected = seekers::BruteForce<Rational>(max_diff).seek(matrix);
  const auto result =
      seekers::SmallRows<Rational, RationalHasher>::seek(matrix, params).first;

  ASSERT_EQ(seekers::normalize_result(result).as_set(),
            std::set(expected.begin(), expected.end()));
}

TEST(SmallRowsTests, CompareWithBruteForce1) {
  for (size_t i = 0; i < 50 && i < benchmark_set.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", benchmark_set[i]));

    auto matrix = get_archived<Rational>(benchmark_set[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_with_brute_force(matrix, 1));
  }
}

TEST(SmallRowsTests, CompareWithBruteForce2) {
  for (size_t i = 0; i < 50 && i < benchmark_set.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", benchmark_set[i]));

    auto matrix = get_archived<Rational>(benchmark_set[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_with_brute_force(matrix, 2));
  }
}

TEST(SmallRowsTests, CompareWithBruteForce3) {
  for (size_t i = 0; i < 50 && i < benchmark_set.size(); ++i) {
    SCOPED_TRACE(std::format("problem name: {}", benchmark_set[i]));

    auto matrix = get_archived<Rational>(benchmark_set[i]).A;
    ASSERT_NO_FATAL_FAILURE(compare_with_brute_force(matrix, 3));
  }
}

TEST(SmallRowsTests, RandomizedTest1) {
  constexpr size_t rows_size = 4;

  std::default_random_engine engine;
  std::uniform_int_distribution<int> elements_distribution(-5, 5);

  for (size_t i = 0; i < 10000; ++i) {
    CSCMatrix<Rational> matrix(2);

    for (size_t j = 0; j < rows_size; ++j) {
      matrix.add_column();

      matrix.push_to_last_column(0, elements_distribution(engine));
      matrix.push_to_last_column(1, elements_distribution(engine));
    }

    seekers::SmallRowsParameters params{
        .max_diff = 2,
        .max_row_size = rows_size,
        .entries_reduction = true,
    };

    auto result =
        seekers::SmallRows<Rational, RationalHasher>::seek(matrix, params)
            .first;
    auto normalized = seekers::normalize_result(result);

    if (similarity::hamming(matrix.get_row(0), matrix.get_row(1)).first <= 2) {
      ASSERT_EQ(normalized.singular.size(), 1);
      ASSERT_EQ(normalized.singular[0], (std::pair{0, 1}));
    } else {
      ASSERT_EQ(normalized.singular.size(), 0);
    }
  }
}

TEST(SmallRowsTests, RandomizedTest2) {
  std::default_random_engine random;

  for (size_t i = 0; i < 100'000; ++i) {
    SCOPED_TRACE(std::format("i = {}", i));

    auto problem = generate_random_problem<Rational>(5, 10, random);

    ASSERT_NO_FATAL_FAILURE(compare_with_brute_force(problem.A, 2));
  }
}

TEST(SmallRowsTests, RandomizedTest3) {
  std::default_random_engine random;

  for (size_t i = 0; i < 1'000; ++i) {
    SCOPED_TRACE(std::format("i = {}", i));

    auto problem = generate_random_problem<Rational>(50, 1000, random);

    ASSERT_NO_FATAL_FAILURE(compare_with_brute_force(problem.A, 2));
  }
}
