#pragma once

#include <array>
#include <string>
#include <string_view>

#include "utils/String.h"

namespace mps {

enum class RowType { LESS_THAN, GREATER_THAN, EQUAL, OBJECTIVE };

struct ParsedRow {
  RowType type;
  std::string name;
};

class RowsParser {
  static RowType decode_row_type(const std::string& type) {
    auto trimmed = str::trim(type);

    switch (trimmed[0]) {
      case 'E':
        return RowType::EQUAL;
      case 'L':
        return RowType::LESS_THAN;
      case 'G':
        return RowType::GREATER_THAN;
      case 'N':
        return RowType::OBJECTIVE;
    }

    throw std::runtime_error("Unknown row type.");
  }

 public:
  ParsedRow parse(const std::array<std::string, 6>& parts) {
    return {decode_row_type(parts[0]), parts[1]};
  }
};

}  // namespace mps
