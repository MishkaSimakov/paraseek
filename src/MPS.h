#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <sstream>
#include <unordered_map>

#include "matrix/CSCMatrix.h"
#include "utils/String.h"

enum class MPSFieldsMode { FIXED_WIDTH, SPACE_SEPARATED };

template <typename Field>
class MPSReader {
  enum class SectionType {
    NAME,
    ROWS,
    COLUMNS,
    RHS,
    BOUNDS,
    RANGES,
    OBJECT,
    ENDATA
  };

  struct Column {
    std::vector<std::pair<size_t, Field>> rows;
  };

  enum class RowType { LESS_THAN, GREATER_THAN, EQUAL, OBJECTIVE };

  MPSFieldsMode mode_;
  std::unordered_map<std::string, Column> columns_;
  std::unordered_map<std::string, size_t> rows_enumeration_;

  std::string objective_row_;

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
      default:
        throw std::runtime_error("Unknown row type.");
    }
  }

  std::optional<std::string> get_marker_type(
      const std::array<std::string, 6>& parts) {
    if (parts[2] != "'MARKER'") {
      return std::nullopt;
    }

    return parts[3];
  }

  std::array<std::string, 6> get_parts(const std::string& str) const {
    constexpr size_t kNumFields = 6;
    constexpr size_t kFieldStartPos[kNumFields] = {1, 4, 14, 24, 39, 49};
    constexpr size_t kFieldLength[kNumFields] = {2, 8, 8, 12, 8, 12};

    std::array<std::string, 6> result;

    if (mode_ == MPSFieldsMode::FIXED_WIDTH) {
      for (size_t i = 0; i < kNumFields; ++i) {
        size_t start = kFieldStartPos[i];
        size_t length = kFieldLength[i];

        if (start >= str.size()) {
          break;
        }
        if (start + length > str.size()) {
          length = std::string::npos;
        }

        result[i] = str::rtrim(str.substr(start, length));
      }
    } else if (mode_ == MPSFieldsMode::SPACE_SEPARATED) {
      result[0] = str::rtrim(str.substr(kFieldStartPos[0], kFieldLength[0]));

      size_t current_index = 1;
      for (size_t i = kFieldStartPos[1]; i < str.size(); ++i) {
        if (i > 0 && std::isspace(str[i]) != 0 &&
            std::isspace(str[i - 1]) == 0) {
          ++current_index;
        } else if (std::isspace(str[i]) == 0) {
          result[current_index] += str[i];
        }
      }
    } else {
      throw std::runtime_error("Unknown field mode in MPSReader.");
    }

    return result;
  }

  static std::optional<SectionType> read_header_card(const std::string& line) {
    std::vector headers = {
        std::pair{"NAME", SectionType::NAME},
        std::pair{"ROWS", SectionType::ROWS},
        std::pair{"COLUMNS", SectionType::COLUMNS},
        std::pair{"RHS", SectionType::RHS},
        std::pair{"BOUNDS", SectionType::BOUNDS},
        std::pair{"RANGES", SectionType::RANGES},
        std::pair{"OBJECT", SectionType::OBJECT},
        std::pair{"ENDATA", SectionType::ENDATA},
    };

    for (const auto& [name, value] : headers) {
      if (line.starts_with(name)) {
        return value;
      }
    }

    throw std::runtime_error("Unknown header card.");
  }

  static Field parse_field(const std::string& str) {
    Field result;
    std::stringstream iss(str);
    iss >> result;

    return result;
  }

  bool should_skip_line(const std::string& line) {
    if (line.empty() || line[0] == '*') {
      return true;
    }

    if (str::ltrim(line).empty()) {
      return true;
    }

    return false;
  }

 public:
  explicit MPSReader(MPSFieldsMode mode) : mode_(mode) {}

  void read(const std::filesystem::path& filepath) {
    std::ifstream is(filepath);

    if (!is) {
      throw std::runtime_error("Failed to open file in MPS reader.");
    }

    std::string line;

    std::optional<SectionType> current_section = std::nullopt;

    while (std::getline(is, line)) {
      if (should_skip_line(line)) {
        continue;
      }

      // header card
      if (line[0] != ' ') {
        current_section = read_header_card(line);

        if (current_section == SectionType::ENDATA) {
          break;
        }

        continue;
      }

      // otherwise proceed with parsing of the current section
      if (!current_section.has_value()) {
        throw std::runtime_error("Section data must be inside of section.");
      }

      auto parts = get_parts(line);

      if (current_section == SectionType::ROWS) {
        if (decode_row_type(parts[0]) == RowType::OBJECTIVE) {
          objective_row_ = parts[1];
        }
      } else if (current_section == SectionType::COLUMNS) {
        auto marker = get_marker_type(parts);

        // skip marker rows
        if (marker.has_value()) {
          continue;
        }

        std::string variable_name = parts[1];

        for (size_t i = 2; i < parts.size(); i += 2) {
          std::string row_name = parts[i];

          if (row_name.empty() || row_name == objective_row_) {
            break;
          }

          auto [itr, _] =
              rows_enumeration_.emplace(row_name, rows_enumeration_.size());

          columns_[variable_name].rows.emplace_back(itr->second,
                                                    parse_field(parts[i + 1]));
        }
      }
    }
  }

  CSCMatrix<Field> get_A() const {
    CSCMatrix<Field> result(rows_enumeration_.size());

    for (const auto& [_, column] : columns_) {
      result.add_column();

      for (auto [index, coef] : column.rows) {
        result.push_to_last_column(index, coef);
      }
    }

    return result;
  }
};
