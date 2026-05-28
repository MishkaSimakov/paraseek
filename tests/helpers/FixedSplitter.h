#pragma once

#include <cassert>
#include <vector>

#include "matrix/CSCMatrix.h"

class FixedSplitter {
  const std::vector<std::vector<size_t>> groups_;

 public:
  explicit FixedSplitter(std::vector<std::vector<size_t>> groups)
      : groups_(std::move(groups)) {}

  template <typename Field>
  std::vector<std::vector<size_t>> split(const CSCMatrix<Field>& /* matrix */,
                                         size_t groups_count) {
    assert(groups_.size() == groups_count);

    return groups_;
  }
};
