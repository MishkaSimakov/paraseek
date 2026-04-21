#include "scip/scip.h"

#include "problems/ProblemsNames.h"
#include "scip/scipdefplugins.h"
#include "utils/Paths.h"

int main() {
  for (size_t problem_index = 0; problem_index < collection_set.size();
       ++problem_index) {
    const auto& problem_name = collection_set[problem_index];
    std::println("{}/{}: {}", problem_index + 1, collection_set.size(),
                 problem_name);

    SCIP* scip = NULL;

    SCIPcreate(&scip);
    SCIPincludeDefaultPlugins(scip);

    // Read MPS
    auto path = paths::problem(problem_name + ".mps");
    SCIPreadProb(scip, path.c_str(), NULL);

    // Presolve
    SCIPpresolve(scip);

    // Write transformed problem
    auto output_path =
        paths::problem(std::format("presolved_{}.mps", problem_name));
    SCIPwriteTransProblem(scip, output_path.c_str(), "mps", FALSE);

    // Cleanup
    SCIPfree(&scip);
  }

  return 0;
}
