#pragma once

#include <functional>

#include "BigRows.h"
#include "Result.h"
#include "SmallRows.h"
#include "splitters/Random.h"
#include "utils/Time.h"

namespace seekers {

struct AllRowsParameters {
  size_t max_diff;
  size_t threshold;

  // big rows params
  size_t groups_count;

  // small rows params
  bool entries_reduction = true;
  size_t small_column_limit = 2;

  std::function<void(std::vector<size_t>)> entries_count_logger;
};

struct AllRowsStatistics {
  size_t pairs_considered;
  size_t entries_considered;

  timing::Duration small_rows_duration;
  timing::Duration big_rows_duration;
};

template <typename Field, typename Hasher,
          typename Splitter = splitters::Random<Field>>
class AllRows {
 public:
  static std::pair<Result, AllRowsStatistics> seek(
      const CSCMatrix<Field>& matrix, const AllRowsParameters& params) {
    // seek big rows
    BigRowsParameters big_params{
        .min_row_size = params.threshold + 1,
        .groups_count = params.groups_count,
        .max_diff = params.max_diff,
    };

    auto [big_result, big_stats] =
        BigRows<Field, Hasher, Splitter>::seek(matrix, big_params);

    // seek small rows
    SmallRowsParameters small_params{
        .max_row_size = params.threshold + params.max_diff,
        .max_diff = params.max_diff,

        .entries_reduction = params.entries_reduction,
        .small_column_limit = params.small_column_limit,

        .entries_count_logger = params.entries_count_logger,
    };

    auto [small_result, small_stats] =
        SmallRows<Field, Hasher>::seek(matrix, small_params);

    // merge statistics
    AllRowsStatistics stats{
        .pairs_considered = big_stats.pairs_considered,
        .entries_considered = small_stats.entries_considered,

        .small_rows_duration = small_stats.duration,
        .big_rows_duration = big_stats.duration,
    };

    // merge results
    Result result = std::move(small_result);

    result.singular.append_range(big_result.singular);
    result.bipartite.append_range(big_result.bipartite);

    return {result, stats};
  }
};

}  // namespace seekers
