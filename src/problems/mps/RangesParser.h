#pragma once

#include <array>
#include <string>

#include "Common.h"
#include "Format.h"

namespace mps {

template <typename Field>
struct ParsedRanges {
  std::string name;

  std::array<std::pair<std::string, Field>, 2> rows;
};

template <typename Field>
class RangesParser {
 public:
  ParsedRanges<Field> parse(const std::array<std::string, 6>& parts,
                            Format format) {
    size_t start = 1;

    if (!parts[0].empty()) {
      if (format == Format::FIXED) {
        throw std::runtime_error(
            "In fixed format MPS in RANGES section names must start from "
            "column 2.");
      }

      start = 0;
    }

    ParsedRanges<Field> result;
    result.name = parts[start];

    for (size_t i = 0; i < 2; ++i) {
      const auto& row_name = parts[start + 1 + 2 * i];

      if (row_name.empty()) {
        break;
      }

      result.rows[i] = {row_name, parse_field<Field>(parts[start + 2 + 2 * i])};
    }

    return result;
  }
};

}  // namespace mps
