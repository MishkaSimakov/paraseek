#pragma once

#include <ostream>
#include <tuple>

struct CsvWriter {
  std::ostream& os;

  template <typename... Args>
  CsvWriter& operator<<(const std::tuple<Args...>& row) {
    std::apply(
        [&](const auto&... values) {
          const char* separator = "";
          ((os << separator << values, separator = ","), ...);
        },
        row);
    os << '\n';
    return *this;
  }

  void flush() { os.flush(); }
};
