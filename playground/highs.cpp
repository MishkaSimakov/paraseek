#include <Highs.h>

#include "problems/ArchivedProblem.h"
#include "problems/ProblemsNames.h"
#include "utils/Paths.h"

void highs_presolve(const std::filesystem::path& input_path,
                    const std::filesystem::path& output_path) {
  Highs highs;
  HighsStatus status;

  status = highs.readModel(input_path);
  if (status != HighsStatus::kOk && status != HighsStatus::kWarning) {
    throw std::runtime_error("Failed to read MPS: " +
                             highsStatusToString(status));
  }

  status = highs.presolve();
  if (status != HighsStatus::kOk && status != HighsStatus::kWarning) {
    throw std::runtime_error("Presolve failed: " + highsStatusToString(status));
  }

  if (highs.getModelPresolveStatus() != HighsPresolveStatus::kReduced &&
      highs.getModelPresolveStatus() != HighsPresolveStatus::kNotReduced) {
    throw std::runtime_error("Unexpected presolve result.");
  }

  status = highs.writePresolvedModel(output_path);
  if (status != HighsStatus::kOk) {
    throw std::runtime_error("Failed to write presolved model: " +
                             highsStatusToString(status));
  }
}

int main() {
  for (size_t problem_index = 0; problem_index < benchmark_set.size();
       ++problem_index) {
    try {
      const auto& problem_name = benchmark_set[problem_index];

      std::println("{}/{}: {}", problem_index + 1, benchmark_set.size(),
                   problem_name);

      auto input_name = std::format("{}.mps", problem_name);
      auto output_name = std::format("presolved_{}.mps", problem_name);

      // check that problem is unarchived
      unarchive(problem_name);

      highs_presolve(paths::problem(input_name), paths::problem(output_name));
    } catch (const std::exception& exception) {
      std::cerr << exception.what() << std::endl;
    }
  }
}
