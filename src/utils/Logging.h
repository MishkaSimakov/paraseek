#pragma once

#include <fstream>
#include <string_view>

#include "Paths.h"

namespace logging {

inline std::ofstream get_log_fstream(std::string_view filename) {
  auto path = paths::log(filename);

  std::filesystem::create_directories(path.parent_path());

  std::ofstream os(path, std::ofstream::app);

  if (!os) {
    throw std::runtime_error("Failed to write to log file.");
  }

  return os;
}

template <typename T>
void log_value(const T& value, std::string_view filename) {
  get_log_fstream(filename) << value << "\n";
}

void log(std::string_view text, std::string_view filename) {
  auto os = get_log_fstream(filename);
  os << text;
}

}  // namespace logging
