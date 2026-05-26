#include "Result.h"

#include <set>

seekers::SingularResult seekers::normalize_result(const Result& result) {
  std::set<std::pair<size_t, size_t>> pairs;

  const auto emplace_pair = [&](size_t i, size_t j) {
    if (i > j) {
      pairs.emplace(j, i);
    } else if (i < j) {
      pairs.emplace(i, j);
    }
  };

  for (auto [i, j] : result.singular) {
    emplace_pair(i, j);
  }

  for (const auto& [left, right] : result.bipartite) {
    for (size_t i : left) {
      for (size_t j : right) {
        emplace_pair(i, j);
      }
    }
  }

  return SingularResult{
      .singular = {pairs.begin(), pairs.end()},
  };
}

std::vector<size_t> seekers::only_one(size_t rows_count, const Result& result) {
  std::vector<size_t> chosen(rows_count, -1);

  for (const auto [i, j] : result.singular) {
    if (chosen[j] == -1) {
      chosen[j] = i;
    } else if (chosen[i] == -1) {
      chosen[i] = j;
    }
  }

  for (const auto& [left, right] : result.bipartite) {
    for (size_t i : left) {
      for (size_t j : right) {
        if (i == j) {
          continue;
        }

        if (chosen[j] == -1) {
          chosen[j] = i;
        } else if (chosen[i] == -1) {
          chosen[i] = j;
        }
      }
    }
  }

  return chosen;
}
