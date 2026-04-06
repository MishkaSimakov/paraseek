#pragma once

#include <sstream>
#include <string>

namespace mps {

template <typename Field>
static Field parse_field(const std::string& str) {
  Field result;
  std::stringstream iss(str);
  iss >> result;

  return result;
}



}  // namespace mps
