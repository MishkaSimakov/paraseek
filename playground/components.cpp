#include <cassert>
#include <chrono>
#include <map>
#include <print>

#include "ExpressionDisjointSet.h"
#include "Reducer.h"
#include "problems/ProblemMatrix.h"
#include "problems/ProblemsNames.h"
#include "seekers/BruteForce.h"
#include "seekers/Tables.h"
#include "seekers/TablesSimpleHashing.h"
#include "utils/Printing.h"

using namespace std::chrono_literals;

void dfs(size_t row, const std::vector<SparseVector<double>>& rows,
         const CSCMatrix<double>& matrix, std::vector<size_t>& components,
         size_t component) {
  std::vector<size_t> stack;
  stack.push_back(row);

  while (!stack.empty()) {
    size_t curr = stack.back();
    components[row] = component;

    for (const size_t col : rows[row] | std::views::keys) {
      for (const size_t next : matrix.get_column(col) | std::views::keys) {
        if (next != row && components[next] == -1) {
          dfs(next, rows, matrix, components, component);
        }
      }
    }
  }
}

int main() {
  std::ofstream os(paths::log("components_count.csv"));
  std::println(os, "problem_name,count");

  constexpr size_t max_size = 10;

  for (size_t problem_index = 0; problem_index < problems_names.size();
       ++problem_index) {
    std::println("{}/{}: {}", problem_index + 1, problems_names.size(),
                 problems_names[problem_index]);

    auto problem = get_problem_matrices(problems_names[problem_index], true);

    const auto [n, d] = problem.A.shape();

    std::println("  size: {} x {} (nz = {})", n, d, problem.A.nonzero_count());

    auto transposed = problem.A.get_transposed();

    std::vector<size_t> component(n, -1);
    size_t components_cnt = 0;

    for (size_t row = 0; row < n; ++row) {
      if (transposed[row].size() > max_size) {
        continue;
      }

      if (component[row] != -1) {
        continue;
      }

      dfs(row, transposed, problem.A, component, components_cnt);
      ++components_cnt;
    }

    std::println("  components: {}", components_cnt);
    std::println(os, "{},{}", problems_names[problem_index], components_cnt);
  }
}
