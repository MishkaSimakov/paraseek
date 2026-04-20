#include <Highs.h>

#include "Reducer.h"
#include "problems/Problem.h"
#include "seekers/Tables.h"
#include "utils/Paths.h"

Problem<double> load_and_presolve(const std::string& filename) {
  Highs highs;

  // 1. Read MPS
  if (highs.readModel(filename) != HighsStatus::kOk) {
    throw std::runtime_error("Failed to read MPS");
  }

  // 2. Run presolve
  if (highs.presolve() != HighsStatus::kOk) {
    throw std::runtime_error("Presolve failed");
  }

  // 3. Check infeasibility
  auto model_status = highs.getModelStatus();
  if (model_status == HighsModelStatus::kInfeasible) {
    Problem<double> prob = Problem<double>::with_size(0, 0);
    prob.proven_unfeasible = true;
    return prob;
  }

  // 4. Extract presolved LP
  if (!highs.getPresolvedModel().isMip()) {
    throw std::runtime_error("Only LP models are supported.");
  }

  HighsLp lp = highs.getPresolvedModel().lp_;

  const size_t n = lp.num_row_;
  const size_t d = lp.num_col_;

  // 5. Convert matrix (CSC format)
  // HiGHS already stores matrix in column-wise format
  CSCMatrix<double> A(n);
  for (int j = 0; j < d; ++j) {
    A.add_column();

    int col_start = lp.a_matrix_.start_[j];
    int col_end = lp.a_matrix_.start_[j + 1];

    for (int idx = col_start; idx < col_end; ++idx) {
      int row = lp.a_matrix_.index_[idx];
      double val = lp.a_matrix_.value_[idx];

      A.push_to_last_column(row, val);
    }
  }

  // 6. Row bounds
  std::vector<Bound<double>> rhs_bounds(n);

  for (int i = 0; i < n; ++i) {
    const auto lb = lp.row_lower_[i];
    const auto ub = lp.row_upper_[i];

    rhs_bounds[i].lower = lb <= -kHighsInf ? std::nullopt : std::optional(lb);
    rhs_bounds[i].upper = ub >= kHighsInf ? std::nullopt : std::optional(ub);
  }

  // 7. Column bounds
  std::vector<Bound<double>> bounds(d);

  for (int j = 0; j < d; ++j) {
    const auto lb = lp.col_lower_[j];
    const auto ub = lp.col_upper_[j];

    bounds[j].lower = lb <= -kHighsInf ? std::nullopt : std::optional(lb);
    bounds[j].upper = ub >= kHighsInf ? std::nullopt : std::optional(ub);
  }

  // 8. Objective
  std::vector<double> c = lp.col_cost_;

  double shift = lp.offset_;

  // 9. Construct problem
  return Problem<double>(std::move(A), std::move(rhs_bounds), std::move(c),
                         std::move(bounds), shift);
}

int main() {
  auto problem = load_and_presolve(paths::problem("square47.mps"));

  std::println("problem: {} x {}", problem.A.shape().first,
               problem.A.shape().second);

  seekers::TablesParameters params{
      .groups_count = 6,
      .max_small_row_size = 8,
  };

  auto seeker = seekers::Tables<double, seekers::DoubleHasher>(3, params);
  auto result = seeker.seek(problem.A);

  std::println("  size = {}", seekers::normalize_result(result).singular.size());
}
