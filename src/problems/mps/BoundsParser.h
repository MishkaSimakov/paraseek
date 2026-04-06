#pragma once

#include <array>
#include <string>

#include "Common.h"
#include "Format.h"
#include "problems/Bound.h"

namespace mps {

template <typename Field>
struct ParsedBounds {
  std::string name;

  std::string variable_name;
  Bound<Field> bound;
  bool is_integer;
};

template <typename Field>
class BoundsParser {
  std::pair<Bound<Field>, bool> parse_bound_type(const std::string& type,
                                                 std::optional<Field> value) {
    bool is_integer = false;
    Bound<Field> bound;

    if (type == "LO") {
      bound = {value.value(), std::nullopt};
    } else if (type == "UP") {
      bound = {std::nullopt, value.value()};
    } else if (type == "FX") {
      bound = {value.value(), value.value()};
    } else if (type == "FR") {
      bound = {std::nullopt, std::nullopt};
    } else if (type == "MI") {
      bound = {std::nullopt, 0};
    } else if (type == "PL") {
      bound = {0, std::nullopt};
    } else if (type == "LI") {
      bound = {value.value(), std::nullopt};
      is_integer = true;
    } else if (type == "UI") {
      bound = {std::nullopt, value.value()};
      is_integer = true;
    } else if (type == "BV") {
      bound = {0, 1};
      is_integer = true;
    } else {
      throw std::runtime_error(std::format("Unsupported bound type: {}", type));
    }

    return {bound, is_integer};
  }

 public:
  ParsedBounds<Field> parse(const std::array<std::string, 6>& parts,
                            Format format) {
    ParsedBounds<Field> result;

    result.name = parts[1];
    result.variable_name = parts[2];

    auto value = parts[3].empty() ? std::nullopt
                                  : std::optional{parse_field<Field>(parts[3])};
    auto parsed_type = parse_bound_type(parts[0], value);

    result.bound = parsed_type.first;
    result.is_integer = parsed_type.second;

    return result;
  }
};

}  // namespace mps
