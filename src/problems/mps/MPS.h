#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <sstream>
#include <unordered_map>

#include "BoundsParser.h"
#include "ColumnsParser.h"
#include "Common.h"
#include "Format.h"
#include "RangesParser.h"
#include "RowsParser.h"
#include "matrix/Matrix.h"
#include "problems/Bound.h"
#include "problems/Problem.h"
#include "utils/String.h"
#include "utils/Variant.h"

namespace mps {

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

  enum class ObjectiveType { MINIMIZE, MAXIMIZE };

  struct Row {
    RowType type;
    // std::vector<std::pair<std::string, Field>> variables;
    Field rhs{0};
    std::optional<Field> range = std::nullopt;

    explicit Row(RowType type) : type(type) {}
  };

  struct VariableInfo {
    Bound<Field> bound = Bound<Field>(0, std::nullopt);
    bool is_integer = false;
    std::vector<std::pair<std::string, Field>> rows;
  };

  Format mode_;

  ObjectiveType objective_ = ObjectiveType::MINIMIZE;
  std::unordered_map<std::string, Row> rows_;
  std::unordered_map<std::string, VariableInfo> variables_;

  bool is_integer_section_ = false;

  std::array<std::string, 6> get_parts(const std::string& str) const {
    constexpr size_t kNumFields = 6;
    constexpr size_t kFieldStartPos[kNumFields] = {1, 4, 14, 24, 39, 49};
    constexpr size_t kFieldLength[kNumFields] = {2, 8, 8, 12, 8, 12};

    std::array<std::string, 6> result;

    if (mode_ == Format::FIXED) {
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
    } else if (mode_ == Format::FREE) {
      size_t current_index = 0;

      for (size_t i = 0; i < str.size(); ++i) {
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
  explicit MPSReader(Format mode) : mode_(mode) {}

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
        const auto parsed = RowsParser().parse(parts);

        rows_.emplace(parsed.name, Row(parsed.type));
      } else if (current_section == SectionType::COLUMNS) {
        const auto parsed = ColumnsParser<Field>().parse(parts, mode_);

        std::visit(
            Overload{
                [this](const ParsedMarkerColumn<Field>& marker) {
                  switch (marker.type) {
                    case MarkerType::INTORG:
                      is_integer_section_ = true;
                      break;
                    case MarkerType::INTEND:
                      is_integer_section_ = false;
                      break;
                    default:
                      throw std::runtime_error("Unsupported marker type.");
                  }
                },
                [this](const ParsedVariableColumn<Field>& variable) {
                  auto [itr, _] =
                      variables_.emplace(variable.name, VariableInfo{});
                  itr->second.is_integer = is_integer_section_;

                  for (auto [name, value] : variable.rows) {
                    if (name.empty()) {
                      break;
                    }

                    itr->second.rows.emplace_back(name, value);
                  }
                },
            },
            parsed);
      } else if (current_section == SectionType::RHS) {
        // rhs vector name is ignored
        // in free MPS format names can start from column 1
        size_t start = 2;

        if (!parts[0].empty()) {
          if (mode_ == Format::FIXED) {
            throw std::runtime_error(
                "In fixed format MPS in RHS section names must start from "
                "column 2.");
          }

          start = 1;
        }

        for (size_t i = start; i < parts.size(); i += 2) {
          std::string row_name = parts[i];
          if (row_name.empty()) {
            continue;
          }

          rows_.at(row_name).rhs = parse_field<Field>(parts[i + 1]);
        }
      } else if (current_section == SectionType::BOUNDS) {
        auto parsed = BoundsParser<Field>().parse(parts, mode_);
        auto& variable_info = variables_.at(parsed.variable_name);

        variable_info.bound = parsed.bound;
        variable_info.is_integer = parsed.is_integer;
      } else if (current_section == SectionType::RANGES) {
        auto parsed = RangesParser<Field>().parse(parts, mode_);

        for (const auto& [name, value] : parsed.rows) {
          if (name.empty()) {
            break;
          }

          rows_.at(name).range = value;
        }
      } else if (current_section == SectionType::OBJECT) {
        // skip for now
        continue;
      } else {
        throw std::runtime_error("Unknown section type.");
      }
    }
  }

  // TODO: support for ranges
  Problem<Field> get_problem(bool replace_inequalities) const {
    std::unordered_map<std::string, size_t> rows_enumeration;
    rows_enumeration.reserve(rows_.size());

    for (const auto& [name, row] : rows_) {
      if (row.type != RowType::OBJECTIVE) {
        rows_enumeration.emplace(name, rows_enumeration.size());
      }
    }

    CSCMatrix<Field> A(rows_enumeration.size());
    std::vector<Field> c(variables_.size(), 0);
    std::vector<Bound<Field>> bounds(variables_.size());

    size_t i = 0;
    for (const auto& [_, info] : variables_) {
      A.add_column();

      for (const auto& [row_name, coef] : info.rows) {
        const auto& row = rows_.at(row_name);

        if (row.type == RowType::OBJECTIVE) {
          c[i] = coef;
          continue;
        }

        A.push_to_last_column(rows_enumeration.at(row_name), coef);
      }

      ++i;
    }

    if (replace_inequalities) {
      for (const auto& [name, row] : rows_) {
        if (row.type == RowType::OBJECTIVE || row.type == RowType::EQUAL) {
          continue;
        }

        const size_t index = rows_enumeration.at(name);

        c.push_back(0);
        bounds.push_back(Bound<Field>(0, std::nullopt));

        if (row.type == RowType::LESS_THAN) {
          A.add_column();
          A.push_to_last_column(index, 1);
        } else if (row.type == RowType::GREATER_THAN) {
          A.add_column();
          A.push_to_last_column(index, -1);
        }
      }
    }

    // rhs column
    std::vector<Field> b(rows_enumeration.size());

    for (const auto& [name, row] : rows_) {
      if (row.type == RowType::OBJECTIVE) {
        continue;
      }

      const size_t index = rows_enumeration.at(name);
      b[index] = row.rhs;
    }

    return {A, b, c, bounds};
  }

  // generates a problem suitable for simplex method:
  // c x -> max, s.t. Ax = b, l <= x <= u
  // MILPProblem<Field> get_canonical_representation() {
  //   MILPProblem<Field> result;
  //   std::unordered_map<std::string, Variable<Field>> variables;
  //
  //   for (const Row& row : rows_ | std::views::values) {
  //     for (const std::string& name : row.variables | std::views::keys) {
  //       if (!variables.contains(name)) {
  //         const VariableInfo& variable = variables_.at(name);
  //
  //         auto variable_type =
  //             variable.is_integer ? VariableType::INTEGER :
  //             VariableType::REAL;
  //
  //         auto var = result.new_variable(
  //             name, variable_type, variable.bound.lower,
  //             variable.bound.upper);
  //         variables.emplace(name, var);
  //       }
  //     }
  //   }
  //
  //   auto objective_row = std::ranges::find_if(
  //       rows_, [](auto p) { return p.second.type == RowType::OBJECTIVE; });
  //
  //   if (objective_row == rows_.end()) {
  //     throw std::runtime_error("No objective row found.");
  //   }
  //
  //   Expression<Field> objective;
  //   for (auto [name, coef] : objective_row->second.variables) {
  //     if (objective_ == ObjectiveType::MINIMIZE) {
  //       coef = -coef;
  //     }
  //
  //     objective += variables.at(name) * coef;
  //   }
  //
  //   result.set_objective(objective);
  //
  //   // process constraints
  //   for (const Row& row : rows_ | std::views::values) {
  //     if (row.type == RowType::OBJECTIVE) {
  //       continue;
  //     }
  //
  //     Expression<Field> lhs;
  //     for (auto [name, coef] : row.variables) {
  //       lhs += variables.at(name) * coef;
  //     }
  //
  //     if (!row.range) {
  //       if (row.type == RowType::LESS_THAN) {
  //         result.add_constraint(lhs <= Expression{row.rhs});
  //       } else if (row.type == RowType::GREATER_THAN) {
  //         result.add_constraint(lhs >= Expression{row.rhs});
  //       } else {
  //         result.add_constraint(lhs == Expression{row.rhs});
  //       }
  //     } else {
  //       Field upper;
  //       Field lower;
  //
  //       if (row.type == RowType::LESS_THAN) {
  //         upper = row.rhs;
  //         lower = row.rhs - *row.range;
  //       } else if (row.type == RowType::GREATER_THAN) {
  //         upper = row.rhs + *row.range;
  //         lower = row.rhs;
  //       } else {
  //         upper = std::max(row.rhs, row.rhs + *row.range);
  //         lower = std::min(row.rhs, row.rhs + *row.range);
  //       }
  //
  //       result.add_constraint(lhs <= Expression{upper});
  //       result.add_constraint(lhs >= Expression{lower});
  //     }
  //   }
  //
  //   return result;
  // }
};

}  // namespace mps
