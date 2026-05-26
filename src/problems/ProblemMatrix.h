#pragma once

#include <format>

#include "HighsMPS.h"
#include "matrix/CSCMatrix.h"
#include "mps/MPS.h"
#include "utils/Paths.h"

// template <typename Field = double>
// Problem<Field> get_problem(std::string name,
//                            bool replace_inequalities = false) {
//   auto archive_path = paths::problem(name + ".mps.gz");
//   auto path = paths::problem(name + ".mps");
//
//   if (!std::filesystem::exists(path)) {
//     std::string command = std::format("gzip -dk {}", archive_path.string());
//     int rs = std::system(command.c_str());
//
//     if (rs != 0) {
//       throw std::runtime_error(std::format(
//           "Failed to decompress problem file: \"{}\"",
//           archive_path.string()));
//     }
//   }
//
//   auto reader = mps::MPSReader<Field>(mps::Format::FREE);
//   reader.read(path);
//
//   return reader.get_problem(replace_inequalities);
// }

namespace detail {

template <typename Field>
void add_slack_variables(Problem<Field>& problem) {
  const auto [n, d] = problem.A.shape();

  for (size_t i = 0; i < n; ++i) {
    if (problem.rhs_bounds[i].is_fixed()) {
      continue;
    }

    problem.c.push_back(0);
    problem.bounds.push_back(-problem.rhs_bounds[i]);

    problem.rhs_bounds[i] = Bound<Field>{0, 0};

    problem.A.add_column();
    problem.A.push_to_last_column(i, 1);
  }
}

}  // namespace detail

template <typename Field = double>
Problem<Field> get_problem(const std::string& name,
                           bool replace_inequalities = false) {
  auto archive_path = paths::problem(name + ".mps.gz");
  auto path = paths::problem(name + ".mps");

  if (!std::filesystem::exists(path)) {
    std::string command = std::format("gzip -dk {}", archive_path.string());
    int rs = std::system(command.c_str());

    if (rs != 0) {
      throw std::runtime_error(std::format(
          "Failed to decompress problem file: \"{}\"", archive_path.string()));
    }
  }

  auto problem = mps::read_via_highs<Field>(path);

  if (replace_inequalities) {
    detail::add_slack_variables(problem);
  }

  return problem;
}
