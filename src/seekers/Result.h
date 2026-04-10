#pragma once

#include <set>
#include <vector>

namespace seekers {

struct SingularResult {
  std::vector<std::pair<size_t, size_t>> singular;

  std::set<std::pair<size_t, size_t>> as_set() const {
    return {singular.begin(), singular.end()};
  }
};

struct Result {
  std::vector<std::pair<size_t, size_t>> singular;
  std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> bipartite;
};

SingularResult normalize_result(const Result& result);

std::vector<size_t> only_one(size_t rows_count, const Result& result);

}  // namespace seekers
