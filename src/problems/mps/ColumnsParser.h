#pragma once

#include <array>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <variant>

#include "Common.h"
#include "Format.h"

namespace mps {

enum class MarkerType { INTORG, INTEND };

template <typename Field>
struct ParsedMarkerColumn {
  MarkerType type;
};

template <typename Field>
struct ParsedVariableColumn {
  std::string name;

  std::array<std::pair<std::string, Field>, 2> rows;
};

template <typename Field>
using ParsedColumn =
    std::variant<ParsedMarkerColumn<Field>, ParsedVariableColumn<Field>>;

template <typename Field>
class ColumnsParser {
  static std::optional<MarkerType> get_marker_type(
      const std::array<std::string, 6>& parts, Format format) {
    size_t start = 2;

    if (!parts[0].empty()) {
      if (format == Format::FIXED) {
        throw std::runtime_error(
            "In fixed format MPS in RHS section names must start from "
            "column 2.");
      }

      start = 1;
    }

    if (parts[start] != "'MARKER'") {
      return std::nullopt;
    }

    if (parts[start + 1] == "'INTORG'") {
      return MarkerType::INTORG;
    }
    if (parts[start + 1] == "'INTEND'") {
      return MarkerType::INTEND;
    }

    throw std::runtime_error(
        std::format("Unknown marker type: {}", parts[start + 1]));
  }

  static ParsedVariableColumn<Field> parse_variable_column(
      const std::array<std::string, 6>& parts, Format format) {
    size_t start = 1;

    if (!parts[0].empty()) {
      if (format == Format::FIXED) {
        throw std::runtime_error(
            "In fixed format MPS in COLUMNS section names must start from "
            "column 2.");
      }

      start = 0;
    }

    ParsedVariableColumn<Field> result;
    result.name = parts[start];

    for (size_t i = 0; i < 2; ++i) {
      const auto& row_name = parts[start + 1 + i * 2];

      if (row_name.empty()) {
        break;
      }

      result.rows[i] = {row_name, parse_field<Field>(parts[start + 2 + i * 2])};
    }

    return result;
  }

 public:
  ParsedColumn<Field> parse(const std::array<std::string, 6>& parts,
                            Format format) {
    if (const auto marker = get_marker_type(parts, format)) {
      return ParsedMarkerColumn<Field>{*marker};
    }

    return parse_variable_column(parts, format);
  }
};

}  // namespace mps
