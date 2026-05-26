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

  // adds rows pair (i, j) to the singular part of the result
  void add(size_t i, size_t j);

  // removes all duplicates from the singular part
  void purge_singular();

  std::set<std::pair<size_t, size_t>> as_set() const;

  SingularResult to_singular() const;
};

std::vector<size_t> only_one(size_t rows_count, const Result& result);

}  // namespace seekers
