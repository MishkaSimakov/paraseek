#pragma once

#include <sstream>
#include <string>

namespace mps {

template <typename Field>
Field parse_field(const std::string& str) {
  Field result;
  std::stringstream iss(str);
  iss >> result;

  return result;
}

template <>
inline double parse_field(const std::string& str) {
  return std::stod(str);
}

}  // namespace mps
